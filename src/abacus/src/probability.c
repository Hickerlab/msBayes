/**
 * @file        probability.c
 * @authors     Jamie Oaks
 * @package     ABACUS (Approximate BAyesian C UtilitieS)
 * @brief       A collection of functions for random variate generation.
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


#include "probability.h"

double draw_gamma_or_uniform(const gsl_rng * rng, double shape, double scale) {
    double draw;
    if ((shape > 0.0) && (scale > 0.0)) {
        draw = gsl_ran_gamma(rng, shape,
                scale);
    }
    else {
        double tau_a, tau_b;
        tau_a = fabs(shape);
        tau_b = fabs(scale);
        if (tau_a < tau_b) {
            draw = gsl_ran_flat(rng, tau_a, tau_b);
        }
        else {
            draw = gsl_ran_flat(rng, tau_b, tau_a);
        }
    }
    return draw;
}

double get_gamma_or_uniform_mean(double shape, double scale) {
    double mn;
    if ((shape > 0.0) && (scale > 0.0)) {
        mn = shape * scale;
    }
    else {
        double tau_a, tau_b;
        tau_a = fabs(shape);
        tau_b = fabs(scale);
        mn = ((tau_a + tau_b) / 2);
    }
    return mn;
}

