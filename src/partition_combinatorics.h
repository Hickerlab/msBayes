/**
 * @file        partition_combinatorics.h
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

#ifndef PARTITION_COMBINATORICS_H
#define PARTITION_COMBINATORICS_H

#include "array_utils.h"

int cumulative_number_of_int_partitions_by_k(int n, i_array * dest);
int number_of_int_partitions_by_k(int n, i_array * dest);
double cumulative_frequency_of_int_partitions_by_k(int n, d_array * probs);
double frequency_of_int_partitions_by_k(int n, d_array * probs);
int number_of_int_partitions(int n);

/** 
 * A function for generating all partitions of an integer.
 *
 * This function generates all of the integer partitions of an integer in
 * anti-lexicographic order. The code is based on Algorithm ZS1 from:
 *
 * Zoghbi, A., and I. Stojmenovic. 1998. Fast algorithms for generating
 *     generating integer partitions. Intern. J. Computer Math.
 *     70:319--332.
 *
 * @param n
 *   The integer to partition.
 * @return
 *   The function returns a pointer to an `i_array_2d` struct. This is an array
 *   of `i_array`s. Each `i_array` contains a distinct integer partition of
 *   `n`; there will be `number_of_int_partitions(n)` of these.
 *   @see number_of_int_partitions()
 */ 
i_array_2d * generate_int_partitions(int n);

#endif /* PARTITION_COMBINATORICS_H */

