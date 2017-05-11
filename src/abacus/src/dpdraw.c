/**
 * @file        dpdraw.c
 * @authors     Jamie Oaks
 * @package     ABACUS (Approximate BAyesian C UtilitieS)
 * @brief       Random draws from a Dirichlet process.
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

#include "dpdraw.h"

config * init_config() {
    config * c;
    c = (typeof(*c) *) malloc(sizeof(*c));
    c->reps = 1;
    c->alpha = 1.0;
    srand(time(NULL));
    c->seed = rand();
    return c;
}

void free_config(config * c) {
    free(c);
    c = NULL;
}

void dpdraw_preamble() {
    char * version = DPDRAW_VERSION;
    char * ab_preamble = abacus_preamble();
    fprintf(stderr, "%s\n", ab_preamble);
    fprintf(stderr, "DPDraw Version %s\n\n", version);
    fprintf(stderr,
        "    Simulating draws from a Dirichlet process\n\n");
    free(ab_preamble);
}

void help() {
    dpdraw_preamble();
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  dpdraw [ -r REPS -a ALPHA -s SEED ] NUM_ELEMENTS\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, " -r  Number of replicates. Default 1.\n");
    fprintf(stderr,
        " -a  Alpha (concentration) parameter of the Dirichlet process.\n"
        "     Default 1.0.\n");
    fprintf(stderr, " -s  Random number seed.\n");
    fprintf(stderr, " -h  Display this help message and exit\n");
}

void parse_args(config * conf, int argc, char ** argv) {
    int i;
    char * end_ptr;
    char * end_ptr_orig;
    /* opterr = 0; */
    end_ptr = (typeof(*end_ptr) *) malloc(sizeof(end_ptr) * 64);
    end_ptr_orig = end_ptr;
    while((i = getopt(argc, argv, "r:a:s:h")) != -1) {
        switch(i) {
            case 'r':
                conf->reps = atoi(optarg);
                break;
            case 'a':
                conf->alpha = strtod(optarg, &end_ptr);
                if (end_ptr == optarg) {
                    fprintf(stderr, "ERROR: %s is not a valid number for "
                            "alpha\n", optarg);
                    exit(1);
                }
                break;
            case 's':
                conf->seed = atoi(optarg);
                break;
            case 'h':
                help();
                exit(0);
                break;
            case '?':
                if (optopt == 'r') {
                    fprintf(stderr, "ERROR: option `-%c' requires an "
                            "argument\n", optopt);
                }
                else if (optopt == 'n') {
                    fprintf(stderr, "ERROR: option `-%c' requires an "
                            "argument\n", optopt);
                }
                else if (optopt == 'a') {
                    fprintf(stderr, "ERROR: option `-%c' requires an "
                            "argument\n", optopt);
                }
                else {
                    fprintf(stderr, "ERROR: unknown option character `\\x%x'\n",
                            optopt);
                }
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
    if (conf->reps < 1) {
        fprintf(stderr, "ERROR: `-r` should be a positive integer\n");
        help();
        exit(1);
    }
    if (conf->alpha <= 0) {
        fprintf(stderr, "ERROR: `-a` should be a positive number\n");
        help();
        exit(1);
    }
    if (conf->seed < 1) {
        fprintf(stderr, "ERROR: `-s` should be a positive integer\n");
        help();
        exit(1);
    }
    free(end_ptr_orig);
}

int dpdraw_main(int argc, char ** argv) {
    int i, num_cats;
    config * conf;
    const gsl_rng_type * mt = gsl_rng_mt19937;
    i_array * elements;

    if (argc < 2) {
        help();
        exit(1);
    }
    conf = init_config();
    parse_args(conf, argc, argv);

    gsl_rng * rng = gsl_rng_alloc(mt);
    gsl_rng_set(rng, conf->seed);
    elements = init_i_array(conf->num_elements);

    for (i = 0; i < conf->reps; i++) {
        num_cats= dirichlet_process_draw(rng,
                conf->num_elements,
                conf->alpha,
                elements);
        assert((num_cats > 0) && (num_cats <= conf->num_elements));
        write_i_array(stdout, elements, "\t");
    }
    gsl_rng_free(rng);
    free_i_array(elements);
    free_config(conf);
    return 0;
}
