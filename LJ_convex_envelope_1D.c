#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "LJ_convex_envelope_1D.h"

double calculate_potential(double squared_distance){
    /*
    Calculates the LJ potential in terms of the squared interatomic distance, x^(-6) - 2*(x^(-3))
    */
    double second_term = pow(squared_distance,-3);
    return pow(second_term,2) - 2*second_term;
}

double get_1D_slope_difference(double x, double upper_bound){
    /*
    Calculates the difference between the slope of the 1D L-J potential x^(-6) - 2*(x^(-3))
    at a squared distance x, and the average slope between x and upper_bound. 
    Used in calculating the critical distance.
    Input:
        double x: point at which to evaluate the exact 1D L-J derivative (in terms of squared distance)
        double upper_bound: right endpoint of interval over which to evalue average derivative of 1D L-J
    Output:
        (double) difference in slopes
    */
    double potential_at_upper_bound = calculate_potential(upper_bound);
    double potential_at_x = calculate_potential(x);
    return 6 * (pow(x,-4) - pow(x,-7)) - (potential_at_upper_bound - potential_at_x)/(upper_bound - x);
}

double get_1D_derivative(double x, double upper_bound){
    /*
    Calculates the derivative of the difference between the slope of the 1D L-J potential x^(-6) - 2*(x^(-3))
    at a squared distance x, and the average slope between x and upper_bound. 
    Used in calculating the critical distance.
    Input:
        double x: point at which to evaluate the exact 1D L-J derivative (in terms of squared distance)
        double upper_bound: right endpoint of interval over which to evalue average derivative of 1D L-J
    Output:
        (double) derivative of the difference in slopes
    */
    double potential_at_upper_bound = calculate_potential(upper_bound);
    double potential_at_x = calculate_potential(x);
    double sq_distance_difference = upper_bound - x;
    return 42*pow(x,-8) - 24*pow(x,-5) - 6*(pow(x,-7) - pow(x,-4))/sq_distance_difference
        - (potential_at_upper_bound - potential_at_x)/(pow(sq_distance_difference,2));
}

double calculate_critical_distance(double upper_bound){
    /*
    Implements Newton's method with an intelligent initial guess
    Empirically, 3 iterations is usually enough to get within 1e-16 of the solution obtained by bracketing
    Input: upper bound for 1D envelope (in terms of squared distance)
    Output: point in [1,(7/4)^(1/3)] for tangent hyperplane
    */
    double newton_tolerance = 1e-16;
    int max_newton_iterations = 10;
    int current_iteration = 0;

    // Initial guess (empirical)
    double current_root;
    if (upper_bound < 2){
        current_root = 1.0 + 0.08956*pow(upper_bound - 0.6825,-1.277);
    }
    else{
        current_root = 1.0 + 0.06411*pow(upper_bound - 0.9904,-1.043);
    }

    while ((fabs(get_1D_slope_difference(current_root,upper_bound)) > newton_tolerance) && (current_iteration < max_newton_iterations)){
        current_root = current_root - 
            get_1D_slope_difference(current_root,upper_bound)/get_1D_derivative(current_root,upper_bound);
        current_iteration += 1;
    }

    return current_root;
}

void calculate_supporting_hyperplane_1D(double x_lower, double x, double x_upper, double *x_coeff_return, double *rhs_return){
    /*
    Input: bounds and value of squared distance, pointers to write results
        The potential in terms of the squared distance is x^(-6) - 2*(x^(-3))
    Output: coefficients of supporting hyperplane to the convex envelope of the LJ potential, in the form
        a0*d + a1*f <= rhs.
        Note: output coefficients are in the space of the epigraph; i.e. the last coefficient is that of the variable representing the objective function
        We don't worry about a sparse representation of the coefficients, since there are only two of them and the second one is always -1
        The only point at which the first coefficient is zero is when the squared distance is exactly 1
    */
    double x_coefficient; // order: x, f
    double current_underestimator_value;
    if ((x_upper <= inflection_point) || (x <= 1)){ // tangent at x
        current_underestimator_value = calculate_potential(x);
        x_coefficient = 6*(pow(x,-4) - pow(x,-7));
    }
    else if (x_lower >= inflection_point){ // secant between xL, xU
        double distance_fraction = (x - x_lower) / (x_upper - x_lower);
        double potential_upper = calculate_potential(x_upper);
        double potential_lower = calculate_potential(x_lower);
        current_underestimator_value = potential_lower + distance_fraction * (potential_upper - potential_lower);
        x_coefficient = (potential_upper - potential_lower) / (x_upper - x_lower);
    }
    else{
        double critical_distance = calculate_critical_distance(x_upper);
        if (x <= critical_distance){ // tangent at x
            current_underestimator_value = calculate_potential(x);
            x_coefficient = 6*(pow(x,-4) - pow(x,-7));
        }
        else if (critical_distance <= x_lower){  // secant between xL, xU
            double distance_fraction = (x - x_lower) / (x_upper - x_lower);
            double potential_upper = calculate_potential(x_upper);
            double potential_lower = calculate_potential(x_lower);
            current_underestimator_value = potential_lower + distance_fraction * (potential_upper - potential_lower);
            x_coefficient = (potential_upper - potential_lower) / (x_upper - x_lower);
        }
        else{  // tangent at xC, equivalent to secant between xC, xU
            double distance_fraction = (x - critical_distance) / (x_upper - critical_distance);
            double potential_upper = calculate_potential(x_upper);
            double potential_critical = calculate_potential(critical_distance);
            current_underestimator_value = potential_critical + distance_fraction * (potential_upper - potential_critical);
            x_coefficient = (potential_upper - potential_critical) / (x_upper - critical_distance);
        }
    }
    double right_hand_side = x_coefficient * x - current_underestimator_value;
    
    *x_coeff_return = x_coefficient;
    *rhs_return = right_hand_side;
}