/**
 * @file        array_utils.c
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

#include "array_utils.h"

d_array * init_d_array(int capacity) {
    assert(capacity > 0);
    d_array * v;
    v = (typeof(*v) *) malloc(sizeof(*v));
    v->capacity = capacity;
    if ((v->a = (typeof(*v->a) *) calloc(v->capacity,
            sizeof(*v->a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    v->length = 0;
    return v;
}

c_array * init_c_array(int capacity) {
    assert(capacity > 0);
    c_array * v;
    v = (typeof(*v) *) malloc(sizeof(*v));
    v->capacity = capacity;
    if ((v->a = (typeof(*v->a) *) calloc((v->capacity + 1),
            sizeof(*v->a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    v->a[v->capacity] = '\0';
    return v;
}

i_array * init_i_array(int capacity) {
    assert(capacity > 0);
    i_array * v;
    v = (typeof(*v) *) malloc(sizeof(*v));
    v->capacity = capacity;
    if ((v->a = (typeof(*v->a) *) calloc(v->capacity,
            sizeof(*v->a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    v->length = 0;
    return v;
}

s_array * init_s_array(int capacity) {
    assert(capacity > 0);
    int i;
    s_array * v;
    v = (typeof(*v) *) malloc(sizeof(*v));
    v->capacity = capacity;
    if ((v->a = (typeof(*v->a) *) calloc(v->capacity,
            sizeof(*v->a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    for (i = 0; i < v->capacity; i++) {
        v->a[i] = init_c_array(64);
    }
    v->length = 0;
    return v;
}

void set_d_array(d_array * v, int index, double x) {
    assert((index >= 0) && (index < v->length));
    v->a[index] = x;
}

void assign_c_array(c_array * v, const char * s) {
    while (strlen(s) > v->capacity) {
        expand_c_array(v);
    }
    strncpy(v->a, s, v->capacity);
}

void set_i_array(i_array * v, int index, int x) {
    assert((index >= 0) && (index < v->length));
    v->a[index] = x;
}

void set_s_array(s_array * v, int index, const char * s) {
    assert((index >= 0) && (index < v->length));
    assign_c_array(v->a[index], s);
}

void expand_d_array(d_array * v) {
    v->capacity *= 2;
    if ((v->a = (typeof(*v->a) *) realloc(v->a ,
            (v->capacity * sizeof(*v->a)))) == NULL) {
        perror("out of memory");
        exit(1);
    }
}

void expand_c_array(c_array * v) {
    v->capacity = ((v->capacity + 1) * 2) - 1;
    if ((v->a = (typeof(*v->a) *) realloc(v->a ,
            ((v->capacity + 1) * sizeof(*v->a)))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    v->a[v->capacity] = '\0';
}

void expand_i_array(i_array * v) {
    v->capacity *= 2;
    if ((v->a = (typeof(*v->a) *) realloc(v->a ,
            (v->capacity * sizeof(*v->a)))) == NULL) {
        perror("out of memory");
        exit(1);
    }
}

void expand_s_array(s_array * v) {
    int i;
    if ((v->a = (typeof(*v->a) *) realloc(v->a ,
            ((v->capacity * 2) * sizeof(*v->a)))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    for (i = v->capacity; i < (v->capacity * 2); i++) {
        v->a[i] = init_c_array(64);
    }
    v->capacity *= 2;
}

void append_d_array(d_array * v, double x) {
    if (v->length >= v->capacity) {
        expand_d_array(v);
    }
    v->length++;
    set_d_array(v, (v->length - 1), x);
}

void append_i_array(i_array * v, int x) {
    if (v->length >= v->capacity) {
        expand_i_array(v);
    }
    v->length++;
    set_i_array(v, (v->length - 1), x);
}

void append_s_array(s_array * v, const char * x) {
    if (v->length >= v->capacity) {
        expand_s_array(v);
    }
    v->length++;
    assign_c_array(v->a[(v->length - 1)], x);
}

void extend_d_array(d_array * dest, const d_array * to_add) {
    int i;
    for (i = 0; i < to_add->length; i++) {
        append_d_array(dest, to_add->a[i]);
    }
}

void extend_i_array(i_array * dest, const i_array * to_add) {
    int i;
    for (i = 0; i < (*to_add).length; i++) {
        append_i_array(dest, (*to_add).a[i]);
    }
}

void extend_s_array(s_array * dest, const s_array * to_add) {
    int i;
    for (i = 0; i < to_add->length; i++) {
        append_s_array(dest, get_s_array(to_add, i));
    }
}

double get_d_array(const d_array * v, int index) {
    assert((index >= 0) && (index < v->length));
    return (v->a[index]);
}

int get_i_array(const i_array * v, int index) {
    assert((index >= 0) && (index < v->length));
    return (v->a[index]);
}

char * get_s_array(const s_array * v, int index) {
    assert((index >= 0) && (index < v->length));
    return (v->a[index]->a);
}

char * get_c_array(const c_array * v) {
    return (v->a);
}

void write_d_array(FILE * stream, const d_array * v, const char * sep) {
    int i;
    for (i = 0; i < (v->length - 1); i++) {
        fprintf(stream, "%.12lf%s", v->a[i], sep);
    }
    fprintf(stream, "%.12lf\n", v->a[(v->length - 1)]);
}

void write_i_array(FILE * stream, const i_array * v, const char * sep) {
    int i;
    for (i = 0; i < (v->length - 1); i++) {
        fprintf(stream, "%d%s", v->a[i], sep);
    }
    fprintf(stream, "%d\n", v->a[(v->length - 1)]);
}

void write_s_array(FILE * stream, const s_array * v, const char * sep) {
    int i;
    for (i = 0; i < (v->length - 1); i++) {
        fprintf(stream, "%s%s", get_s_array(v, i), sep);
    }
    fprintf(stream, "%s\n", get_s_array(v, (v->length - 1)));
}

void free_d_array(d_array * v) {
    free(v->a);
    free(v);
    v = NULL;
}

void free_c_array(c_array * v) {
    free(v->a);
    free(v);
    v = NULL;
}

void free_i_array(i_array * v) {
    free(v->a);
    free(v);
    v = NULL;
}

void free_s_array(s_array * v) {
    int i;
    for (i = 0; i < v->capacity; i++) {
        free_c_array(v->a[i]);
    }
    free(v->a);
    free(v);
    v = NULL;
}

int almost_equal(const double x, const double y, const double error) {
    return ((fabs(x - y) < error));
}

int d_arrays_equal(const d_array * v1, const d_array * v2, double error) {
    int i;
    if (v1->length != v2->length) {
        return 0;
    }
    for (i = 0; i < v1->length; i++) {
        if (almost_equal(get_d_array(v1, i), get_d_array(v2, i), error) == 0) {
            return 0;
        }
    }
    return 1;
}

int i_arrays_equal(const i_array * v1, const i_array * v2) {
    int i;
    if (v1->length != v2->length) {
        return 0;
    }
    for (i = 0; i < v1->length; i++) {
        if (get_i_array(v1, i) != get_i_array(v2, i)) {
            return 0;
        }
    }
    return 1;
}

int s_arrays_equal(const s_array * v1, const s_array * v2) {
    int i;
    if (v1->length != v2->length) {
        return 0;
    }
    for (i = 0; i < v1->length; i++) {
        if (strcmp(get_s_array(v1, i), get_s_array(v2, i)) != 0) {
            return 0;
        }
    }
    return 1;
}

int split_str(char * string, s_array * words, int expected_num) {
    int word_idx, n;
    char * ptr;
    c_array * match;
    match = init_c_array(63);
    word_idx = 0;
    ptr = string;
    words->length = 0;
    while(*ptr) {
        if ((sscanf(ptr, "%s%n", match->a, &n)) == 1) {
            ptr += n;
            append_s_array(words, match->a);
            word_idx++;
        }
        ++ptr;
    }
    free_c_array(match);
    if ((expected_num > 0) && (expected_num != word_idx)) {
        if (word_idx == 0) word_idx--;
        return word_idx;
    }
    return 0;
}

int split_str_d(char * string, d_array * v, int expected_num) {
    int word_idx, n;
    char * ptr;
    double match;
    word_idx = 0;
    ptr = string;
    v->length = 0;
    while(*ptr) {
        if ((sscanf(ptr, "%lf%n", &match, &n)) == 1) {
            ptr += n;
            append_d_array(v, match);
            word_idx++;
        }
        ++ptr;
    }
    if ((expected_num > 0) && (expected_num != word_idx)) {
        if (word_idx == 0) word_idx--;
        return word_idx;
    }
    return 0;
}

int split_str_i(char * string, i_array * v, int expected_num) {
    int word_idx, n;
    char * ptr;
    int match;
    word_idx = 0;
    ptr = string;
    v->length = 0;
    while(*ptr) {
        if ((sscanf(ptr, "%d%n", &match, &n)) == 1) {
            ptr += n;
            append_i_array(v, match);
            word_idx++;
        }
        ++ptr;
    }
    if ((expected_num > 0) && (expected_num != word_idx)) {
        if (word_idx == 0) word_idx--;
        return word_idx;
    }
    return 0;
}

void get_matching_indices(const s_array * search_strings,
        const s_array * target_strings,
        i_array * indices) {
    int i, j, found;
    indices->length = 0;
    for (i = 0; i < (*search_strings).length; i++) {
        found = 0;
        for (j = 0; j < (*target_strings).length; j++) {
            if (strcmp(get_s_array(search_strings, i),
                    get_s_array(target_strings, j)) == 0) {
                if (found != 0) {
                    fprintf(stderr, "ERROR: get_matching_indices: string %s "
                            "found more than once\n",
                            get_s_array(search_strings, i));
                    exit(1);
                }
                found = 1;
                append_i_array(indices, j);
            }
        }
        if (found == 0) {
            fprintf(stderr, "ERROR: get_matching_indices: string %s was not "
                    "found\n", get_s_array(search_strings, i));
            exit(1);
        }
    }
}

int get_doubles(const s_array * strings, const i_array * indices,
        d_array * doubles_dest) {
    int i, ret;
    char * end_ptr;
    char * end_ptr_orig;
    end_ptr = (typeof(*end_ptr) *) malloc(sizeof(end_ptr) * 64);
    end_ptr_orig = end_ptr;
    ret = 0;
    (*doubles_dest).length = 0;
    for (i = 0; i < (*indices).length; i++) {
        append_d_array(doubles_dest, strtod(
                get_s_array(strings, (*indices).a[i]),
                &end_ptr));
        if (end_ptr == get_s_array(strings, (*indices).a[i])) {
            fprintf(stderr, "ERROR: get_doubles : column %d is not a valid "
                    "number\n", ((*indices).a[i] + 1));
            ret++;
        }
    }
    free(end_ptr_orig);
    return ret;
}

i_array_2d * init_i_array_2d(int capacity, int initial_element_capacity) {
    assert(capacity > 0);
    int i;
    i_array_2d * v;
    v = (typeof(*v) *) malloc(sizeof(*v));
    v->capacity = capacity;
    v->initial_element_capacity = initial_element_capacity;
    if ((v->a = (typeof(*v->a) *) calloc(v->capacity,
            sizeof(*v->a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    for (i = 0; i < v->capacity; i++) {
        v->a[i] = init_i_array(initial_element_capacity);
    }
    v->length = 0;
    return v;
}

void set_i_array_2d(i_array_2d * v, int index, const i_array * x) {
    assert((index >= 0) && (index < v->length));
    v->a[index]->length = 0;
    extend_i_array(v->a[index], x);
}

void set_el_i_array_2d(i_array_2d * v, int i_array_index, int el_index, int x) {
    set_i_array(get_i_array_2d(v, i_array_index), el_index, x);
}

void expand_i_array_2d(i_array_2d * v) {
    int i;
    if ((v->a = (typeof(*v->a) *) realloc(v->a ,
            ((v->capacity * 2) * sizeof(*v->a)))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    for (i = v->capacity; i < (v->capacity * 2); i++) {
        v->a[i] = init_i_array(v->initial_element_capacity);
    }
    v->capacity *= 2;
}

void append_i_array_2d(i_array_2d * v, const i_array * x) {
    if (v->length >= v->capacity) {
        expand_i_array_2d(v);
    }
    v->length++;
    set_i_array_2d(v, (v->length - 1), x);
}

void append_el_i_array_2d(i_array_2d * v, int i_array_index, int x) {
    if (i_array_index >= v->length) {
        while (i_array_index >= v->capacity) {
            expand_i_array_2d(v);
        }
        v->length = i_array_index + 1;
    }
    append_i_array(get_i_array_2d(v, i_array_index), x);
}

void extend_i_array_2d(i_array_2d * dest, const i_array_2d * to_add) {
    int i;
    for (i = 0; i < to_add->length; i++) {
        append_i_array_2d(dest, get_i_array_2d(to_add, i));
    }
}

i_array * get_i_array_2d(const i_array_2d * v, int index) {
    assert((index >= 0) && (index < v->length));
    return (v->a[index]);
}

int get_el_i_array_2d(const i_array_2d * v, int i_array_index, int el_index) {
    return(get_i_array(get_i_array_2d(v, i_array_index), el_index));
}

void free_i_array_2d(i_array_2d * v) {
    int i;
    for (i = 0; i < v->capacity; i++) {
        free_i_array(v->a[i]);
    }
    free(v->a);
    free(v);
    v = NULL;
}

