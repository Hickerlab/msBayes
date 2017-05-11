/**
 * @file        parsing.h
 * @authors     Jamie Oaks
 * @package     ABACUS (Approximate BAyesian C UtilitieS)
 * @brief       A collection of parsing functions.
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

#ifndef PARSING_H
#define PARSING_H

#include <stdio.h>
#include <string.h>

#include "array_utils.h"


void parse_header(const char * path, c_array * line_buffer, s_array * header);
void parse_observed_stats_file(const char * path, c_array * line_buffer,
        s_array * header, d_array * stats);
void parse_summary_file(const char * path,
        c_array * line_buffer,
        s_array * header,
        d_array * means,
        d_array * std_devs,
        i_array * sample_sizes);
int strcmp_i(const char * a, const char * b);
char * strip(const char * s);

#endif /* PARSING_H */

