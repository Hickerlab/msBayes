/**
 * @file        intpart.c
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

#include "intpart.h"

config * init_config() {
    config * c;
    c = (typeof(*c) *) malloc(sizeof(*c));
    c->summary = 1;
    return c;
}

void free_config(config * c) {
    free(c);
    c = NULL;
}

void intpart_preamble() {
    char * version = INTPART_VERSION;
    char * ab_preamble = abacus_preamble();
    fprintf(stderr, "%s\n", ab_preamble);
    fprintf(stderr, "IntPart Version %s\n\n", version);
    fprintf(stderr,
        "    Integer partitioning\n\n");
    free(ab_preamble);
}

void help() {
    intpart_preamble();
    printf("Usage:\n");
    printf("  intpart [ -a ] N\n\n");
    printf("Options:\n");
    printf(" -a  Report all partitions. The default is to show only a\n");
    printf("     summary of the number and relative frequency of partitions\n");
    printf("     per number of categories. NOTE: only use this for\n");
    printf("     relatively small numbers of elements (e.g., < 50).\n");
    printf(" -h  Display this help message and exit\n");
}

void parse_args(config * conf, int argc, char ** argv) {
    int i;
    /* opterr = 0; */
    while((i = getopt(argc, argv, "ah")) != -1) {
        switch(i) {
            case 'a':
                conf->summary = 0;
                break;
            case 'h':
                help();
                exit(0);
                break;
            case '?':
                fprintf(stderr, "ERROR: unknown option character `\\x%x'\n",
                        optopt);
                break;
            default:
                help();
                exit(1);
        }
    }
    if ((optind + 1) > argc) {
        fprintf(stderr, "ERROR: Please provide the number of elements "
                "as an argument\n");
        help();
        exit(1);
    }
    else if ((optind + 1) < argc) {
        fprintf(stderr, "ERROR: Too many positional arguments\n");
        help();
        exit(1);
    }
    else {
        conf->num_elements = atoi(argv[optind]);
    }
    // vetting
    if (conf->num_elements < 1) {
        fprintf(stderr, "ERROR: expecting a positive integer as an argument "
                "(the number of elements)\n");
        help();
        exit(1);
    }
}

int intpart_main(int argc, char ** argv) {
    int i, ip;
    double total_prob;
    config * conf;

    if (argc < 2) {
        help();
        exit(1);
    }
    conf = init_config();
    parse_args(conf, argc, argv);

    if (conf->summary!= 0) {
        i_array * header;
        i_array * counts;
        d_array * probs;
        header = init_i_array(conf->num_elements);
        counts = init_i_array(conf->num_elements);
        probs = init_d_array(conf->num_elements);
        for (i = 1; i <= conf->num_elements; i++) {
            append_i_array(header, i);
        }
        write_i_array(stdout, header, "\t");
        ip = number_of_int_partitions_by_k(conf->num_elements, counts);
        total_prob = frequency_of_int_partitions_by_k(conf->num_elements, probs);
        assert(almost_equal(total_prob, 1.0, 0.000001));
        write_i_array(stdout, counts, "\t");
        write_d_array(stdout, probs, "\t");
        free_i_array(header);
        free_i_array(counts);
        free_d_array(probs);
        free_config(conf);
        return 0;
    }

    i_array_2d * partitions;
    partitions = generate_int_partitions(conf->num_elements);
    for (i = 0; i < partitions->length; i++) {
        write_i_array(stdout, get_i_array_2d(partitions, i), "\t");
    }
    free_i_array_2d(partitions);
    free_config(conf);
    return(0);
}
