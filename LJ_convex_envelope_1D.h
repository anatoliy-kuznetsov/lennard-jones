#ifndef LJ_CONVEX_ENVELOPE_1D
#define LJ_CONVEX_ENVELOPE_1D
double calculate_potential(double squared_distance);
void calculate_supporting_hyperplane_1D(double x_lower, double x, double x_upper, double *x_coeff_return, double *rhs_return);
#endif

#ifndef inflection_point
#define inflection_point 1.205071132087615
#endif