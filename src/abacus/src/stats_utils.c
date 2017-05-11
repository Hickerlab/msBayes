/**
 * @file        stats_utils.c
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

#include "stats_utils.h"

sample_sum * init_sample_sum() {
    sample_sum * ss;
    ss = (typeof(*ss) *) malloc(sizeof(*ss));
    ss->n = 0;
    ss->sum = 0.0;
    ss->sum_of_squares = 0.0;
    return ss;
}

void free_sample_sum(sample_sum * s) {
    free(s);
    s = NULL;
}

sample_sum_array * init_sample_sum_array(int length) {
    assert(length > 0);
    sample_sum_array * v;
    v = (typeof(*v) *) malloc(sizeof(*v));
    v->length = length;
    if ((v->a = (typeof(*v->a) *) calloc(v->length,
            sizeof(*v->a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    int i;
    for (i = 0; i < v->length; i++) {
        v->a[i] = init_sample_sum();
    }
    return v;
}

void free_sample_sum_array(sample_sum_array * v) {
    int i;
    for (i = 0; i < v->length; i++) {
        free_sample_sum(v->a[i]);
    }
    free(v->a);
    free(v);
    v = NULL;
}

void update_sample_sum(sample_sum * s, double x) {
    (*s).n += 1;
    (*s).sum += x;
    (*s).sum_of_squares += pow(x, 2);
}

double get_mean(const sample_sum * s) {
    assert((*s).n > 0);
    double mn = ((*s).sum / (*s).n);
    return mn;
}

double get_sample_variance(const sample_sum * s) {
    assert((*s).n > 1);
    double mn = get_mean(s);
    double v = (((*s).sum_of_squares - (mn * (*s).sum)) / ((*s).n - 1));
    return v;
}

double get_std_dev(const sample_sum * s) {
    assert((*s).n > 1);
    double var = get_sample_variance(s);
    double sd = (sqrt(var));
    return sd;
}

void update_sample_sum_array(sample_sum_array * s, const d_array * x) {
    assert((*s).length == (*x).length);
    int i;
    for (i = 0; i < (*s).length; i++) {
        update_sample_sum((*s).a[i], (*x).a[i]);
    }
}

void get_mean_array(const sample_sum_array * s, d_array * means) {
    int i;
    (*means).length = 0;
    for (i = 0; i < (*s).length; i++) {
        append_d_array(means, get_mean((*s).a[i]));
    }
}

void get_sample_variance_array(const sample_sum_array * s,
        d_array * v) {
    int i;
    (*v).length = 0;
    for (i = 0; i < (*s).length; i++) {
        append_d_array(v, get_sample_variance((*s).a[i]));
    }
}
        
void get_std_dev_array(const sample_sum_array * s, d_array * std_devs) {
    int i;
    (*std_devs).length = 0;
    for (i = 0; i < (*s).length; i++) {
        append_d_array(std_devs, get_std_dev((*s).a[i]));
    }
}

void standardize_vector(d_array * v, const d_array * means,
        const d_array * std_devs) {
    assert (((*v).length == (*means).length) &&
            ((*v).length == (*std_devs).length));
    int i;
    for (i = 0; i < (*v).length; i++) {
        (*v).a[i] = ((*v).a[i] - (*means).a[i]) / (*std_devs).a[i];
    }
}

