/**
 * @file        partition_combinatorics_random.c
 * @authors     Jamie Oaks
 * @package     ABACUS (Approximate BAyesian C UtilitieS)
 * @brief       A collection of functions for partition combinatorics.
 * @copyright   Copyright (C) 2013 Jamie Oaks.
 *   This file is part of ABACUS.  ABACUS is free software; you can
 *   redistribute it and/or modify it under the terms of the GNU General Public
 *   License as published by the Free Software Foundation; either version 2 of
 *   the License, or (at your option) any later version.
 * 
 *   ABACUS is distributed in the hope that it will be useful, but WITHOUT ANY
 *   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *   details.
 * 
 *   You should have received a copy of the GNU General Public License along
 *   with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "partition_combinatorics_random.h"

int draw_int_partition_category(const gsl_rng * rng, int n) {
    assert(n > 0);
    d_array * cumulative_probs;
    cumulative_probs = init_d_array(n);
    double total_prob = cumulative_frequency_of_int_partitions_by_k(n,
            cumulative_probs);
    assert(almost_equal(total_prob, 1.0, 0.000001) != 0);
    double r = gsl_rng_uniform(rng);
    int i;
    for (i = 0; i < n; i++) {
        if (r < get_d_array(cumulative_probs, i)) {
            break;
        }
    }
    free_d_array(cumulative_probs);
    return (i+1);
}

/** 
 * A function for generating a random draw from a Dirichlet process.
 */
int dirichlet_process_draw(const gsl_rng * rng, int n, double alpha,
        i_array * elements) {
    assert(n > 0);
    int num_subsets;
    double subset_prob, new_subset_prob, u;
    i_array * subset_counts;
    (*elements).length = 0;
    append_i_array(elements, 0);
    subset_counts = init_i_array(n);
    append_i_array(subset_counts, 1);
    num_subsets = 1;
    int i, j;
    for (i = 1; i < n; i++) {
        new_subset_prob = (alpha / (alpha + (double)i));
        u = gsl_rng_uniform(rng);
        u -= new_subset_prob;
        if (u < 0.0) {
            append_i_array(elements, num_subsets);
            append_i_array(subset_counts, 1);
            num_subsets += 1;
            continue;
        }
        for (j = 0; j < num_subsets; j++) {
            subset_prob = ((double)subset_counts->a[j] / (alpha + (double)i));
            u -= subset_prob;
            if (u < 0.0) {
                append_i_array(elements, j);
                subset_counts->a[j] += 1;
                break;
            }
        }
        if (u > 0.0) {
            append_i_array(elements, (num_subsets-1));
            subset_counts->a[num_subsets-1] += 1;
        }
    }
    free_i_array(subset_counts);
    return num_subsets;
}

