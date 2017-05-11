/**
 * @file        math_utils.c
 * @authors     Jamie Oaks
 * @package     ABACUS (Approximate BAyesian C UtilitieS)
 * @brief       A collection of math types and functions.
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

#include "math_utils.h"

double get_euclidean_distance(const d_array * v1, const d_array * v2) {
    assert((*v1).length == (*v2).length);
    double sum_of_squared_diffs = 0.0;
    int i;
    for (i = 0; i < (*v1).length; i++) {
        sum_of_squared_diffs += pow(((*v1).a[i] - (*v2).a[i]), 2);
    }
    return sqrt(sum_of_squared_diffs);
}

