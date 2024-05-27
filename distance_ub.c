#include "putative_minimum_energies.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

double calculate_energy(int slice_populations[], int slice_count, double slice_difference_energies[]) {
    double energy = 0;
    for (int i = 0; i < slice_count; i++) {
        int interactions_in_slice = slice_populations[i] * (slice_populations[i] - 1) / 2;
        double within_slice_lower_bound = -1.0 * interactions_in_slice;
        energy += within_slice_lower_bound;
        for (int j = i + 1; j < slice_count; j++) {
            double interaction_energy; // between slices i and j
            int interaction_count = slice_populations[i] * slice_populations[j];
            interaction_energy = slice_difference_energies[j - i] * interaction_count;
            energy += interaction_energy;
        }
    }
    return energy;
}

int main(void) {
    int max_atom_count = 150;
    double slice_difference_energies[max_atom_count];
    for (int i = 0; i < max_atom_count; i++) {
        if (i < 3) {
            slice_difference_energies[i] = -1;
            continue;
        }
        slice_difference_energies[i] = pow((double)i, -12) - 2 * pow((double)i, -6);
    }

    printf("N,Distance upper bound\n");
    for (int atom_count = 5; atom_count <= max_atom_count; atom_count++) {
        int slice_count = atom_count - 1;
        bool bound_reached = false;
        while (!bound_reached) {
            int free_atom_count = atom_count - slice_count;
            int slice_populations[slice_count];
            for (int i = 0; i < slice_count; i++) {
                slice_populations[i] = 1;
            }
            slice_populations[slice_count / 2] += free_atom_count;
            double energy = calculate_energy(slice_populations, slice_count, slice_difference_energies);
            if (energy < putative_minimum_energies[atom_count]) {
                bound_reached = true;
                printf("%d,%d\n", atom_count, slice_count);
            }
            slice_count--;
        }
    }
    return 0;
}