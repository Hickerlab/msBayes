/**
 * @file        intpart.h
 * @authors     Jamie Oaks
 * @package     ABACUS (Approximate BAyesian C UtilitieS)
 * @brief       Integer partitioning.
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

#ifndef INTPART_H
#define INTPART_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> // for getopt
#include <string.h>
#include <ctype.h>

#include "array_utils.h"
#include "partition_combinatorics.h"
#include "abacus.h"

#define INTPART_VERSION "0.1.0"

typedef struct config_ {
    int num_elements;
    int summary;
} config;

config * init_config();
void free_config(config * c);

void intpart_preamble();
void help();
void parse_args(config * conf, int argc, char ** argv);

int intpart_main(int argc, char ** argv);

#endif /* INTPART_H */

