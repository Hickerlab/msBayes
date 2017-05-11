/**
 * @file        abacus.c
 * @authors     Jamie Oaks
 * @package     ABACUS (Approximate BAyesian C UtilitieS)
 * @brief       Package info.
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

#include "abacus.h"

char * abacus_preamble() {
    char * s;
    if ((s = (typeof(*s) *) calloc(512, sizeof(*s))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    sprintf(s, "%s %s (%s)\n\n%s\n\n%s\n", ABACUS_NAME, ABACUS_VERSION,
            ABACUS_DESCRIPTION, ABACUS_COPYRIGHT, ABACUS_LICENSE);
    return s;
}

