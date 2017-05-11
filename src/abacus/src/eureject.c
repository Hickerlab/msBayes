/**
 * @file        eureject.c
 * @authors     Jamie Oaks
 * @package     ABACUS (Approximate BAyesian C UtilitieS)
 * @brief       A collection of functions for Euclidean-distance-based
 *              rejection.
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

#include "eureject.h"

config * init_config() {
    config * c;
    c = (typeof(*c) *) malloc(sizeof(*c));
    c->num_retain = 1000;
    c->num_subsample = 10000;
    c->summary_provided = 0;
    c->means = init_d_array(1);
    c->std_devs = init_d_array(1);
    c->sim_paths = init_s_array(1);
    c->observed_path = init_c_array(63);
    c->summary_path = init_c_array(63);
    c->summary_out_path = init_s_array(1);
    c->include_distance = 0;
    return c;
}

void free_config(config * c) {
    free_d_array(c->means);
    free_d_array(c->std_devs);
    free_s_array(c->sim_paths);
    free_c_array(c->observed_path);
    free_c_array(c->summary_path);
    free_s_array(c->summary_out_path);
    free(c);
    c = NULL;
}

sample * init_sample(
        char * file_path,
        const int line_num,
        const s_array * line_array,
        const i_array * stat_indices,
        const d_array * std_observed_stats,
        const d_array * means,
        const d_array * std_devs) {
    int get_stats_return;
    d_array * stats;
    sample * s;
    s = (typeof(*s) *) malloc(sizeof(*s));
    s->file_path = init_c_array(63);
    assign_c_array(s->file_path, file_path);
    s->line_num = line_num;
    s->line_array = init_s_array(line_array->length);
    extend_s_array(s->line_array, line_array);
    stats = init_d_array(stat_indices->length);
    get_stats_return = get_doubles(line_array, stat_indices, stats);
    if (get_stats_return != 0) {
        fprintf(stderr, "ERROR: file %s line %d contains %d invalid stats "
                "columns\n",
                file_path, line_num, get_stats_return);
    }
    standardize_vector(stats, means, std_devs);
    s->distance = get_euclidean_distance(std_observed_stats, stats);
    free_d_array(stats);
    return s;
}

void write_sample(FILE * stream, const sample * s, const int include_distance) {
    if (include_distance != 0) {
        fprintf(stream, "%lf\t", s->distance);
    }
    write_s_array(stream, s->line_array, "\t");
}

void free_sample(sample * s) {
    free_s_array(s->line_array);
    free_c_array(s->file_path);
    free(s);
    s = NULL;
}
    
sample_array * init_sample_array(int capacity) {
    if (capacity < 1) {
        fprintf(stderr, "ERROR: init_d_array: capacity must be positive int "
                "greater than 0\n");
        exit(1);
    }
    sample_array * v;
    v = (typeof(*v) *) malloc(sizeof(*v));
    v->capacity = capacity;
    if ((v->a = (typeof(*v->a) *) calloc(v->capacity,
            sizeof(*v->a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    v->length = 0;
    v->num_processed = 0;
    v->header = init_s_array(64);
    v->paths_processed = init_s_array(1);
    return v;
}

int process_sample(sample_array * samples, sample * s) {
    samples->num_processed++;
    if (samples->length == 0) {
        samples->a[0] = s;
        samples->length++;
        return 0;
    }
    if (samples->a[(samples->length - 1)]->distance <= s->distance) {
        if (samples->length >= samples->capacity) {
            return -1;
        }
        else {
            samples->a[samples->length] = s;
            samples->length++;
            return (samples->length - 1);
        }
    }
    int i;
    for (i = 0; i < samples->length; i++) {
        if (samples->a[i]->distance > s->distance) {
            rshift_samples(samples, i);
            samples->a[i] = s;
            return i;
        }
    }
    return -1;
}

void rshift_samples(sample_array * s, int index) {
    int i, inc;
    inc = 0;
    if (s->length < s->capacity) {
        s->a[s->length] = s->a[(s->length - 1)];
        inc = 1;
    }
    else {
        free_sample(s->a[(s->length - 1)]);
    }
    for (i = (s->length - 1); i > index; i--) {
        s->a[i] = s->a[i - 1];
    }
    s->length += inc;
}

void write_sample_array(FILE * stream, const sample_array * s,
        const int include_distance) {
    int i;
    if (include_distance != 0) {
        fprintf(stream, "distance\t");
    }
    write_s_array(stream, s->header, "\t");
    for (i = 0; i < s->length; i++) {
        write_sample(stream, s->a[i], include_distance);
    }
}

void free_sample_array(sample_array * v) {
    int i;
    for (i = 0; i < v->length; i++) {
        free_sample(v->a[i]);
    }
    free(v->a);
    free_s_array(v->header);
    free_s_array(v->paths_processed);
    free(v);
    v = NULL;
}

void eureject_preamble() {
    char * version = EUREJECT_VERSION;
    char * ab_preamble = abacus_preamble();
    fprintf(stderr, "%s\n", ab_preamble);
    fprintf(stderr, "EuReject Version %s\n\n", version);
    fprintf(stderr,
        "    Euclidean-distance based rejection\n\n");
    free(ab_preamble);
}

void help() {
    eureject_preamble();
    fprintf(stderr, "Usage:\n");
    fprintf(stderr,
        "  eureject -f OBS-FILE [-k INT] [-n INT] [-e] [-s SUM-FILE] \\\n"
        "      [-o SUM-OUT-FILE] SIMS-FILE1 [ SIMS-FILE2 [...] ]\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr,
        " -f  Path to file containing observed summary statistics\n");
    fprintf(stderr,
        " -k  Number of samples to keep. Default: 1000.\n"
        "     If set to 0, only the stat means and standard deviations\n"
        "     are calculated and reported (i.e., no rejection is\n"
        "     performed).\n");
    fprintf(stderr,
        " -n  Number of samples to use for calculating stat means and\n"
        "     standard deviations for standardizing the statistics.\n"
        "     This option is ignored if `-s` is provided. Default: 10000\n");
    fprintf(stderr,
        " -s  Tab-delimited file containing the means and standard\n"
        "     deviations to use for standardizing statistics. The file\n"
        "     must contain the same header as the file with the observed\n"
        "     summary statistics (`-f`), with unique names identifying\n"
        "     the statistics in each column. This header must be\n"
        "     followed by a line containing the means for each statistic\n"
        "     and a third line containing the standard deviations.\n");
    fprintf(stderr,
        " -o  Output file path for means and standard deviations used for\n"
        "     standardizing statistics. The means and standard deviations are\n"
        "     always reported to standard error along with other run info.\n"
        "     It does not make much sense to use this option in combination\n"
        "     with `-s` (it will simply recreate the file specified by `-s`).\n"
        "     This option is nice if you plan to re-use the means/standard\n"
        "     deviations, or just want a permanent record of them.\n");
    fprintf(stderr,
        " -e  Report Euclidean distances of retained samples in the\n"
        "     first column of the output. Default is not to report\n"
        "     the distance column.\n");
    fprintf(stderr, " -h  Display this help message and exit\n");
}

void write_summary(FILE * stream,
        const s_array * sum_paths_processed,
        const int sum_sample_size,
        const s_array * reject_paths_processed,
        const int num_samples_processed,
        const int num_samples_retained) {
    fprintf(stream, "=======\nSUMMARY\n=======\n");
    fprintf(stream, "Files used for calculating means/std deviations: ");
    if (sum_paths_processed->length < 1) {
        fprintf(stream, "None\n");
    }
    else {
        write_s_array(stream, sum_paths_processed, ", ");
    }
    fprintf(stream, "Number of samples used for means/std deviations: %d\n",
            sum_sample_size);
    fprintf(stream, "Files processed for rejection: ");
    if (reject_paths_processed->length < 1) {
        fprintf(stream, "None\n");
    }
    else {
        write_s_array(stream, reject_paths_processed, ", ");
    }
    fprintf(stream, "Total number of samples processed during rejection: %d\n",
            num_samples_processed);
    fprintf(stream, "Number of samples retained: %d\n", num_samples_retained);
}

void write_config(FILE * stream, const config * c) {
    int i;
    fprintf(stream, "========\nSETTINGS\n========\n");
    fprintf(stream, "Number of samples to retain: %d\n", c->num_retain);
    fprintf(stream, "Number of samples to use for standardization: %d\n",
            c->num_subsample);
    fprintf(stream, "Observed stats path: %s\n", c->observed_path->a);
    fprintf(stream, "Path(s) to file(s) with simulated draws: ");
    for (i = 0; i < (c->sim_paths->length - 1); i++) {
        fprintf(stream, "%s, ", get_s_array(c->sim_paths, i));
    }
    fprintf(stream, "%s\n", get_s_array(c->sim_paths,
                (c->sim_paths->length-1)));
    fprintf(stream, "Means for standardization: ");
    if (c->summary_provided == 0) {
        fprintf(stream, "None\n");
    }
    else {
        write_d_array(stream, c->means, ", ");
    }
    fprintf(stream, "Standard deviations for standardization: ");
    if (c->summary_provided == 0) {
        fprintf(stream, "None\n");
    }
    else {
        write_d_array(stream, c->std_devs, ", ");
    }
}

void parse_args(config * conf, int argc, char ** argv) {
    int i, j;
    char * p;
    char * end_ptr;
    char * end_ptr_orig;
    /* opterr = 0; */
    end_ptr = (typeof(*end_ptr) *) malloc(sizeof(end_ptr) * 64);
    end_ptr_orig = end_ptr;
    conf->means->length = 0;
    conf->std_devs->length = 0;
    conf->sim_paths->length = 0;
    while((i = getopt(argc, argv, "f:k:n:s:o:eh")) != -1) {
        switch(i) {
            case 'f':
                assign_c_array(conf->observed_path, optarg);
                break;
            case 'k':
                conf->num_retain = atoi(optarg);
                break;
            case 'n':
                if (atoi(optarg) >= 0) {
                    conf->num_subsample = atoi(optarg);
                    break;
                }
                else {
                    fprintf(stderr, "ERROR: `-n' must be positive integer\n");
                    help();
                    exit(1);
                }
                break;
            case 's':
                conf->summary_provided = 1;
                assign_c_array(conf->summary_path, optarg);
                break;
            case 'o':
                append_s_array(conf->summary_out_path, optarg);
                break;
            case 'e':
                conf->include_distance = 1;
                break;
            case 'h':
                help();
                exit(0);
                break;
            case '?':
                if (optopt == 'o') {
                    fprintf(stderr, "ERROR: option `-%c' requires an "
                            "argument\n", optopt);
                }
                else if (optopt == 'k') {
                    fprintf(stderr, "ERROR: option `-%c' requires an "
                            "argument\n", optopt);
                }
                else if (optopt == 'n') {
                    fprintf(stderr, "ERROR: option `-%c' requires an "
                            "argument\n", optopt);
                }
                else if (optopt == 's') {
                    fprintf(stderr, "ERROR: option `-%c' requires an "
                            "argument\n", optopt);
                }
                else if (isprint(optopt)) {
                    fprintf(stderr, "ERROR: unknown option `-%c'\n", optopt);
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
    for (i = optind, j = 0; i < argc; i++, j++) {
        append_s_array(conf->sim_paths, argv[i]);
    }
    if ((conf->num_subsample < 1) && (conf->summary_provided != 1)) {
        fprintf(stderr, "ERROR: If `-n` is 0, a summary file must be provided "
                "via `-s`\n");
        help();
        exit(1);
    }
    if (conf->summary_provided == 1) {
        conf->num_subsample = 0;
    }
    // vetting
    if (conf->observed_path->a == NULL) {
        fprintf(stderr, "ERROR: Please provide path to observed stats\n");
        help();
        exit(1);
    }
    if ((conf->sim_paths->length < 1) ||
            (get_s_array(conf->sim_paths, 0) == NULL)) {
        fprintf(stderr, "ERROR: Please provide at least one simulation file\n");
        help();
        exit(1);
    }
    free(end_ptr_orig);
}

sample_array * reject(const s_array * paths,
        const c_array * line_buffer,
        const i_array * stat_indices,
        const d_array * std_observed_stats,
        const d_array * means,
        const d_array * std_devs,
        int num_retain,
        const s_array * header) {
    FILE * f;
    int i, line_num, ncols, sample_idx;
    s_array * line_array;
    sample_array * retained_samples;
    line_array = init_s_array((*header).length);
    retained_samples = init_sample_array(num_retain);
    extend_s_array(retained_samples->header, header);
    for (i = 0; i < (*paths).length; i++) {
        line_num = 0;
        if ((f = fopen(get_s_array(paths, i), "r")) == NULL) {
            perror(get_s_array(paths, i));
            exit(1);
        }
        append_s_array(retained_samples->paths_processed,
                get_s_array(paths, i));
        while (fgets((*line_buffer).a, (((*line_buffer).capacity) - 1),
                    f) != NULL) {
            line_num++;
            ncols = split_str((*line_buffer).a, line_array, (*header).length);
            if (ncols == -1) continue; //empty line
            if (ncols != 0) {
                fprintf(stderr, "ERROR: file %s line %d has %d columns "
                        "(expected %d)\n", get_s_array(paths, i), line_num,
                        ncols, (*header).length);
                exit(1);
            }
            if (line_num == 1) continue;
            sample * s;
            s = init_sample(get_s_array(paths, i), line_num, line_array,
                    stat_indices, std_observed_stats, means, std_devs);
            sample_idx = process_sample(retained_samples, s);
            if (sample_idx < 0) {
                free_sample(s);
            }
        }
        fclose(f);
    }
    free_s_array(line_array);
    return retained_samples;
}

void summarize_stat_samples(const s_array * paths,
        c_array * line_buffer,
        const i_array * stat_indices,
        sample_sum_array * ss_array,
        d_array * means,
        d_array * std_devs,
        int num_to_sample,
        int expected_num_columns,
        s_array * paths_processed) {
    assert(ss_array->length == stat_indices->length);
    FILE * f;
    int i, line_num, ncols, get_stats_return;
    s_array * line_array;
    d_array * stats;
    line_array = init_s_array(expected_num_columns);
    stats = init_d_array((*stat_indices).length);
    for (i = 0; i < (*paths).length; i++) {
        line_num = 0;
        if ((f = fopen(get_s_array(paths, i), "r")) == NULL) {
            perror(get_s_array(paths, i));
            exit(1);
        }
        append_s_array(paths_processed, get_s_array(paths, i));
        while (fgets((*line_buffer).a, (((*line_buffer).capacity) - 1),
                    f) != NULL) {
            line_num++;
            ncols = split_str((*line_buffer).a, line_array,
                    expected_num_columns);
            if (ncols == -1) continue; //empty line
            if (ncols != 0) {
                fprintf(stderr, "ERROR: file %s line %d has %d columns "
                        "(expected %d)\n", get_s_array(paths, i), line_num,
                        ncols, expected_num_columns);
                exit(1);
            }
            if (line_num == 1) continue;
            get_stats_return = get_doubles(line_array, stat_indices, stats);
            if (get_stats_return != 0) {
                fprintf(stderr, "ERROR: file %s line %d has %d invalid stats "
                        "columns\n", get_s_array(paths, i), line_num,
                        get_stats_return);
                exit(1);
            }
            update_sample_sum_array(ss_array, stats);
            if (ss_array->a[0]->n >= num_to_sample) break;
        }
        fclose(f);
        if (ss_array->a[0]->n >= num_to_sample) break;
    }
    get_mean_array(ss_array, means);
    get_std_dev_array(ss_array, std_devs);
    free_d_array(stats);
    free_s_array(line_array);
}

int eureject_main(int argc, char ** argv) {
    c_array * line_buffer;
    s_array * obs_header;
    s_array * summary_header;
    s_array * sim_header;
    s_array * sim_header_comp;
    d_array * obs_stats;
    int i, heads_match, sum_sample_size;
    i_array * indices;
    i_array * summary_sample_sizes;
    sample_sum_array * sample_sums;
    sample_array * retained_samples;
    s_array * sum_paths_used;
    config * conf;
    FILE * summary_out_stream;
    line_buffer = init_c_array(pow(2, 14));
    obs_header = init_s_array(1);
    obs_stats = init_d_array(1);
    sum_paths_used = init_s_array(1);
    retained_samples = init_sample_array(1);
    if (argc < 2) {
        help();
        exit(1);
    }
    conf = init_config();
    parse_args(conf, argc, argv);

    parse_observed_stats_file(conf->observed_path->a, line_buffer, obs_header,
            obs_stats);

    summary_sample_sizes = init_i_array(obs_header->length);

    // parse means and std devs from summary  file
    if (conf->summary_provided != 0) {
        summary_header = init_s_array(obs_header->length);
        parse_summary_file(conf->summary_path->a, line_buffer, summary_header,
                conf->means, conf->std_devs, summary_sample_sizes);
        heads_match = s_arrays_equal(obs_header, summary_header);
        if (heads_match == 0) {
            fprintf(stderr, "ERROR: Files %s and %s have different headers\n",
                    get_c_array(conf->observed_path),
                    get_c_array(conf->summary_path));
            help();
            exit(1);
        }
        free_s_array(summary_header);
    }

    sim_header = init_s_array(obs_header->length);
    parse_header(get_s_array(conf->sim_paths, 0), line_buffer, sim_header);
    // check all simulation file headers
    if (conf->sim_paths->length > 1) {
        sim_header_comp = init_s_array(sim_header->length);
        for (i = 1; i < conf->sim_paths->length; i++) {
            parse_header(get_s_array(conf->sim_paths, i), line_buffer,
                    sim_header_comp);
            heads_match = s_arrays_equal(sim_header, sim_header_comp);
            if (heads_match == 0) {
                fprintf(stderr, "ERROR: Files %s and %s have different "
                        "headers\n", get_s_array(conf->sim_paths, 0),
                        get_s_array(conf->sim_paths, i));
                help();
                exit(1);
            }
        }
        free_s_array(sim_header_comp);
    }
    indices = init_i_array(obs_header->length);
    get_matching_indices(obs_header, sim_header, indices);

    eureject_preamble();
    write_config(stderr, conf);

    // calc means and standard devs
    sum_sample_size = 0;
    if (conf->summary_provided == 0) {
        fprintf(stderr, "\nCalculating means and standard deviations... ");
        sample_sums = init_sample_sum_array(obs_header->length);
        summarize_stat_samples(conf->sim_paths, line_buffer,
                indices, sample_sums, conf->means, conf->std_devs,
                conf->num_subsample, sim_header->length, sum_paths_used);
        summary_sample_sizes->length = 0;
        for (i = 0; i < sample_sums->length; i++){
            append_i_array(summary_sample_sizes, sample_sums->a[i]->n);
        }
        sum_sample_size = get_i_array(summary_sample_sizes, 0);
        free_sample_sum_array(sample_sums);
        fprintf(stderr, "Done!\n");
    }

    // rejection
    if (conf->num_retain > 0) {
        free_sample_array(retained_samples);
        fprintf(stderr, "\nPerforming rejection... ");
        standardize_vector(obs_stats, conf->means, conf->std_devs);
        retained_samples = reject(conf->sim_paths, line_buffer, indices,
                obs_stats, conf->means, conf->std_devs, conf->num_retain,
                sim_header);
        fprintf(stderr, "Done!\n\n");
    }

    // write run stats
    write_summary(stderr,
        sum_paths_used,
        sum_sample_size,
        retained_samples->paths_processed,
        retained_samples->num_processed,
        retained_samples->length);

    // write means and standard devs
    if (conf->summary_out_path->length == 1) {
        if ((summary_out_stream = fopen(get_s_array(conf->summary_out_path, 0),
                "w")) == NULL) {
            fprintf(stderr,
                    "ERROR: Could not open %s for writing means and std devs;\n"
                    "using standard error instead\n",
                    get_s_array(conf->summary_out_path, 0));
        }
        else {
            write_s_array(summary_out_stream, obs_header, "\t");
            write_d_array(summary_out_stream, conf->means, "\t");
            write_d_array(summary_out_stream, conf->std_devs, "\t");
            write_i_array(summary_out_stream, summary_sample_sizes, "\t");
            fclose(summary_out_stream);
        }
    }
    fprintf(stderr,
                  "\n=====================================================\n");
    fprintf(stderr, "MEANS AND STD DEVIATIONS USED FOR STANDARDIZING STATS\n");
    fprintf(stderr, "=====================================================\n");
    write_s_array(stderr, obs_header, "\t");
    write_d_array(stderr, conf->means, "\t");
    write_d_array(stderr, conf->std_devs, "\t");
    write_i_array(stderr, summary_sample_sizes, "\t");
    fprintf(stderr, "\n");

    // write retained samples
    if (conf->num_retain > 0) {
        write_sample_array(stdout, retained_samples, conf->include_distance);
    }

    free_sample_array(retained_samples);
    free_i_array(indices);
    free_c_array(line_buffer);
    free_s_array(obs_header);
    free_s_array(sim_header);
    free_d_array(obs_stats);
    free_s_array(sum_paths_used);
    free_config(conf);
    free_i_array(summary_sample_sizes);
    return 0;
}

