/**
 * @file        array_utils.h
 * @authors     Jamie Oaks
 * @package     ABACUS (Approximate BAyesian C UtilitieS)
 * @brief       A collection of types and functions for array manipulations.
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

#ifndef ARRAY_UTILS_H
#define ARRAY_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

typedef struct d_array_ {
    double * a;
    int length;
    int capacity;
} d_array;

typedef struct i_array_ {
    int * a;
    int length;
    int capacity;
} i_array;

typedef struct c_array_ {
    char * a;
    int capacity;
} c_array;

typedef struct s_array_ {
    c_array ** a;
    int length;
    int capacity;
} s_array;

typedef struct i_array_2d_ {
    i_array ** a;
    int length;
    int capacity;
    int initial_element_capacity;
} i_array_2d;

d_array * init_d_array(int length);
void expand_d_array(d_array * v);
void append_d_array(d_array * v, double x);
void extend_d_array(d_array * dest, const d_array * to_add);
double get_d_array(const d_array * v, int index);
void set_d_array(d_array * v, int index, double x);
void write_d_array(FILE * stream, const d_array * v, const char * sep);
void free_d_array(d_array * v);

c_array * init_c_array(int length);
void expand_c_array(c_array * v);
void assign_c_array(c_array * v, const char * s);
void free_c_array(c_array * v);
char * get_c_array(const c_array * v);

i_array * init_i_array(int length);
void expand_i_array(i_array * v);
void append_i_array(i_array * v, int x);
void extend_i_array(i_array * dest, const i_array * to_add);
int get_i_array(const i_array * v, int index);
void set_i_array(i_array * v, int index, int x);
void write_i_array(FILE * stream, const i_array * v, const char * sep);
void free_i_array(i_array * v);

s_array * init_s_array(int length);
void expand_s_array(s_array * v);
void append_s_array(s_array * v, const char * x);
void extend_s_array(s_array * dest, const s_array * to_add);
char * get_s_array(const s_array * v, int index);
void set_s_array(s_array * v, int index, const char * s);
void write_s_array(FILE * stream, const s_array * v, const char * sep);
void free_s_array(s_array * v);

int almost_equal(const double x, const double y, const double error);
int d_arrays_equal(const d_array * v1, const d_array * v2, double error);
int i_arrays_equal(const i_array * v1, const i_array * v2);
int s_arrays_equal(const s_array * v1, const s_array * v2);

int split_str(char * string, s_array * words, int expected_num);
int split_str_d(char * string, d_array * v, int expected_num);
int split_str_i(char * string, i_array * v, int expected_num);
void get_matching_indices(const s_array * search_strings,
        const s_array * target_strings,
        i_array * indices);
int get_doubles(const s_array * strings, const i_array * indices,
        d_array * doubles_dest);

i_array_2d * init_i_array_2d(int capacity, int initial_element_capacity);
void set_i_array_2d(i_array_2d * v, int index, const i_array * x);
void set_el_i_array_2d(i_array_2d * v, int i_array_index, int el_index, int x);
void expand_i_array_2d(i_array_2d * v);
void append_i_array_2d(i_array_2d * v, const i_array * x);
void append_el_i_array_2d(i_array_2d * v, int i_array_index, int x);
void extend_i_array_2d(i_array_2d * dest, const i_array_2d * to_add);
i_array * get_i_array_2d(const i_array_2d * v, int index);
int get_el_i_array_2d(const i_array_2d * v, int i_array_index, int el_index);
void free_i_array_2d(i_array_2d * v);

#endif /* ARRAY_UTILS_H */

