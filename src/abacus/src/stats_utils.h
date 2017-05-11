/**
 * @file        stats_utils.h
 * @authors     Jamie Oaks
 * @package     ABACUS (Approximate BAyesian C UtilitieS)
 * @brief       A collection of types and functions for basic statistics.
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

#ifndef STATS_UTILS_H
#define STATS_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "array_utils.h"

typedef struct sample_sum_ {
    int n;
    double sum;
    double sum_of_squares;
} sample_sum;

typedef struct sample_sum_array_ {
    sample_sum ** a;
    int length;
} sample_sum_array;

sample_sum * init_sample_sum();
sample_sum_array * init_sample_sum_array(int length);
void free_sample_sum_array(sample_sum_array * v);
void update_sample_sum(sample_sum * s, double x);
double get_mean(const sample_sum * s);
double get_sample_variance(const sample_sum * s);
double get_std_dev(const sample_sum * s);
void update_sample_sum_array(sample_sum_array * s, const d_array * x);
void get_mean_array(const sample_sum_array * s, d_array * means);
void get_sample_variance_array(const sample_sum_array * s,
        d_array * v);
void get_std_dev_array(const sample_sum_array * s, d_array * std_devs);
void standardize_vector(d_array * v, const d_array * means,
        const d_array * std_devs);

#endif /* STATS_UTILS_H */

