/**
 * @file        partition_combinatorics.c
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

#include "partition_combinatorics.h"

int cumulative_number_of_int_partitions_by_k(int n,
        i_array * dest) {
    assert(n > 0);
    int i, j;
    int table[n+1][n+1];
    // initialize table with base cases
    for (i = 0; i <= n; i++)
        table[i][0] = 0;
    for (i = 1; i <= n; i++)
        table[0][i] = 1;
    // populate table
    for (i = 1 ; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            if (i - j < 0){
                table[i][j] = table[i][j-1];
                continue;
            }
            table[i][j] = table[i][j-1] + table[i-j][j];
        }
    }
    (*dest).length = 0;
    for (i = 0; i < n; i++) {
        append_i_array(dest, table[n][i+1]);
    }
    return table[n][n];
}

int number_of_int_partitions_by_k(int n, i_array * dest) {
    assert(n > 0);
    int i;
    i_array * v;
    v = init_i_array(n);
    int ip = cumulative_number_of_int_partitions_by_k(n, v);
    (*dest).length = 0;
    append_i_array(dest, 1);
    for (i = 1; i < n; i++) {
        append_i_array(dest, (v->a[i] - v->a[i-1]));
    }
    free_i_array(v);
    return ip;
}
    
double cumulative_frequency_of_int_partitions_by_k(int n,
        d_array * probs) {
    assert(n > 0);
    int i;
    i_array * v;
    v = init_i_array(n);
    int ip = cumulative_number_of_int_partitions_by_k(n, v);
    (*probs).length = 0;
    for (i = 0; i < n; i++) {
        append_d_array(probs, (v->a[i] / ((double) ip)));
    }
    free_i_array(v);
    return (get_d_array(probs, (n-1)));
}

double frequency_of_int_partitions_by_k(int n, d_array * probs) {
    assert(n > 0);
    int i;
    i_array * v;
    v = init_i_array(n);
    int ip = number_of_int_partitions_by_k(n, v);
    (*probs).length = 0;
    double sum = 0.0;
    for (i = 0; i < n; i++) {
        append_d_array(probs, (v->a[i] / ((double) ip)));
        sum += get_d_array(probs, i);
    }
    free_i_array(v);
    return sum;
}

int number_of_int_partitions(int n) {
    if (n < 0) return 0;
    if (n == 0) return 1;
    i_array * v;
    v = init_i_array(n);
    int ip = cumulative_number_of_int_partitions_by_k(n, v);
    free_i_array(v);
    return ip;
}

/** 
 * A function for generating all partitions of an integer.
 */ 
i_array_2d * generate_int_partitions(int n) {
    assert(n > 0);
    int i, ip;
    i_array_2d * partitions;
    ip = number_of_int_partitions(n);
    partitions = init_i_array_2d(ip, n);
    partitions->length = ip;
    int x[n];
    for (i = 0; i < n; i++) {
        x[i] = 1;
    }
    x[0] = n;
    int m = 0;
    int h = 0;
    int r, t, sum;
    append_el_i_array_2d(partitions, 0, x[0]);
    int part_index = 0;
    while (x[0] != 1) {
        part_index += 1;
        if (x[h] == 2) {
            m += 1;
            x[h] = 1;
            h -= 1;
        }
        else {
            r = x[h] - 1;
            t = m - h + 1;
            x[h] = r;
            while (t >= r) {
                h += 1;
                x[h] = r;
                t -= r;
            }
            if (t == 0) {
                m = h;
            }
            else {
                m = h + 1;
                if (t > 1) {
                    h += 1;
                    x[h] = t;
                }
            }
        }
        sum = 0;
        for (i = 0; i < n; i++) {
            append_el_i_array_2d(partitions, part_index, x[i]);
            sum += x[i];
            if (sum >= n) {
                break;
            }
        }
    }
    return partitions;
}

