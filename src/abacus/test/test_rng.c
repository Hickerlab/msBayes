/**
 * @file        test_rng.c
 * @authors     Jamie Oaks
 * @package     ABACUS
 * @brief       Random number generators for testing.
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

#include "test_rng.h"

gsl_rng * get_rng(int seed) {
    const gsl_rng_type * mt = gsl_rng_mt19937;
    gsl_rng * rng = gsl_rng_alloc(mt);
    if (seed < 1) {
        srand(time(NULL));
        gsl_rng_set(rng, rand());
        return rng;
    }
    gsl_rng_set(rng, rand());
    return rng;
}

void free_rng(gsl_rng * rng) {
    gsl_rng_free(rng);
}

