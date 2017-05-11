/**
 * @file        parsing.c
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

#include "parsing.h"

void parse_header(const char * path, c_array * line_buffer, s_array * header) {
    FILE * f;
    if ((f = fopen(path, "r")) == NULL) {
        perror(path);
        exit(1);
    }
    if ((fgets((*line_buffer).a, (((*line_buffer).capacity) - 1), f)) == NULL) {
        fprintf(stderr, "ERROR: found no lines in %s", path);
        exit(1);
    }
    split_str((*line_buffer).a, header, 0);
    fclose(f);
}

void parse_observed_stats_file(const char * path, c_array * line_buffer,
        s_array * header, d_array * stats) {
    FILE * f;
    (*header).length = 0;
    (*stats).length = 0;
    if ((f = fopen(path, "r")) == NULL) {
        perror(path);
        exit(1);
    }
    // parse header
    if ((fgets((*line_buffer).a, (((*line_buffer).capacity) - 1), f)) == NULL) {
        fprintf(stderr, "ERROR: found no header in %s\n", path);
        exit(1);
    }
    split_str((*line_buffer).a, header, 0);
    // parse stats
    if ((fgets((*line_buffer).a, (((*line_buffer).capacity) - 1), f)) == NULL) {
        fprintf(stderr, "ERROR: found no stats in %s\n", path);
        exit(1);
    }
    split_str_d((*line_buffer).a, stats, 0);
    fclose(f);
    if ((*header).length != (*stats).length) {
        fprintf(stderr, "ERROR: found %d column headers, but %d stats in "
                "file %s\n", (*header).length, (*stats).length, path);
        exit(1);
    }
}

void parse_summary_file(const char * path,
        c_array * line_buffer,
        s_array * header,
        d_array * means,
        d_array * std_devs,
        i_array * sample_sizes) {
    FILE * f;
    (*header).length = 0;
    (*means).length = 0;
    (*std_devs).length = 0;
    (*sample_sizes).length = 0;
    if ((f = fopen(path, "r")) == NULL) {
        perror(path);
        exit(1);
    }
    // parse header
    if ((fgets((*line_buffer).a, (((*line_buffer).capacity) - 1), f)) == NULL) {
        fprintf(stderr, "ERROR: found no header in %s\n", path);
        exit(1);
    }
    split_str((*line_buffer).a, header, 0);
    // parse means
    if ((fgets((*line_buffer).a, (((*line_buffer).capacity) - 1), f)) == NULL) {
        fprintf(stderr, "ERROR: found no means in %s\n", path);
        exit(1);
    }
    split_str_d((*line_buffer).a, means, 0);
    if ((*header).length != (*means).length) {
        fprintf(stderr, "ERROR: found %d column headers, but %d means in "
                "file %s\n", (*header).length, (*means).length, path);
        exit(1);
    }
    // parse std deviations
    if ((fgets((*line_buffer).a, (((*line_buffer).capacity) - 1), f)) == NULL) {
        fprintf(stderr, "ERROR: found no std devs in %s\n", path);
        exit(1);
    }
    split_str_d((*line_buffer).a, std_devs, 0);
    if ((*header).length != (*std_devs).length) {
        fprintf(stderr, "ERROR: found %d column headers, but %d std devs in "
                "file %s\n", (*header).length, (*std_devs).length, path);
        exit(1);
    }
    // parse sample sizes
    if ((fgets((*line_buffer).a, (((*line_buffer).capacity) - 1), f)) == NULL) {
        fprintf(stderr, "ERROR: found no sample sizes in %s\n", path);
        exit(1);
    }
    split_str_i((*line_buffer).a, sample_sizes, 0);
    fclose(f);
    if ((*header).length != (*sample_sizes).length) {
        fprintf(stderr, "ERROR: found %d column headers, but %d sample sizes "
                "in file %s\n", (*header).length, (*sample_sizes).length, path);
        exit(1);
    }
}

int strcmp_i(const char * a, const char * b) {
    for (;; a++, b++) {
        int x = tolower(*a) - tolower(*b);
        if (x != 0 || ! *a) {
            return x;
        }
    }
}

char * strip(const char * str) {
    const char *end;
    size_t out_length;
    char *s;
    int len;
    len = strlen(str);
    if ((s = (typeof(*s) *) malloc(sizeof(s) * (len))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    end = str + len - 1;
    while (isspace(*str)) {
        str++;
    }
    if (*str == 0) {
        *s = 0;
        return s;
    }
    
    while ((end > str) && (isspace(*end))) {
        end--;
    }
    end++;
    out_length = (end - str) < len ? (end - str) : len;
    memcpy(s, str, out_length);
    s[out_length] = '\0';
    return s;
}


