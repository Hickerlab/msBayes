#include <stdlib.h>
#include <check.h>
#include <signal.h>
#include "../src/eureject.c"
#include "test_utils.h"


START_TEST (test_init_free_config) {
    config * c;
    c = init_config();
    ck_assert_int_eq(c->summary_provided, 0);
    free_config(c);
}
END_TEST

START_TEST (test_summarize_stat_samples_p1_n4) {
    int i;
    s_array * paths;
    s_array * paths_used;
    c_array * line_buffer;
    i_array * stat_indices;
    sample_sum_array * ss_array;
    d_array * means;
    d_array * std_devs;
    d_array * exp_means;
    d_array * exp_std_devs;
    int num_to_sample;
    int expected_num_cols;
    expected_num_cols = 6;
    num_to_sample = 4;
    paths = init_s_array(1);
    paths_used = init_s_array(1);
    line_buffer = init_c_array(1023);
    stat_indices = init_i_array(4);
    ss_array = init_sample_sum_array(4);
    means = init_d_array(1);
    std_devs = init_d_array(1);
    exp_means = init_d_array(1);
    exp_std_devs = init_d_array(1);
    append_s_array(paths, "data/test_parameter_stat_samples.txt");
    for (i = 2; i < expected_num_cols; i++) {
        append_i_array(stat_indices, i);
    }
    summarize_stat_samples(paths, line_buffer, stat_indices, ss_array, means,
            std_devs, num_to_sample, expected_num_cols, paths_used);
    ck_assert_int_eq(means->length, stat_indices->length);
    ck_assert_int_eq(std_devs->length, stat_indices->length);
    ck_assert_int_eq(ss_array->a[0]->n, num_to_sample);
    ck_assert_msg((s_arrays_equal(paths, paths_used) != 0),
            "paths used do not match input paths");
    append_d_array(exp_means, 0.25);
    append_d_array(exp_means, 0.2225);
    append_d_array(exp_means, 2.5);
    append_d_array(exp_means, 3.0);
    append_d_array(exp_std_devs, 0.1290994);
    append_d_array(exp_std_devs, 0.009574271);
    append_d_array(exp_std_devs, 1.290994);
    append_d_array(exp_std_devs, 1.154701);
    ck_assert_msg((d_arrays_equal(means, exp_means, 0.000001) != 0),
            "unexpected means");
    ck_assert_msg((d_arrays_equal(std_devs, exp_std_devs, 0.000001) != 0),
            "unexpected std deviations");
    free_s_array(paths);
    free_s_array(paths_used);
    free_c_array(line_buffer);
    free_i_array(stat_indices);
    free_sample_sum_array(ss_array);
    free_d_array(means);
    free_d_array(std_devs);
    free_d_array(exp_means);
    free_d_array(exp_std_devs);
}
END_TEST

START_TEST (test_summarize_stat_samples_p1_n3) {
    int i;
    s_array * paths;
    s_array * paths_used;
    c_array * line_buffer;
    i_array * stat_indices;
    sample_sum_array * ss_array;
    d_array * means;
    d_array * std_devs;
    d_array * exp_means;
    d_array * exp_std_devs;
    int num_to_sample;
    int expected_num_cols;
    expected_num_cols = 6;
    num_to_sample = 3;
    paths = init_s_array(1);
    paths_used = init_s_array(1);
    line_buffer = init_c_array(1023);
    stat_indices = init_i_array(4);
    ss_array = init_sample_sum_array(4);
    means = init_d_array(1);
    std_devs = init_d_array(1);
    exp_means = init_d_array(1);
    exp_std_devs = init_d_array(1);
    append_s_array(paths, "data/test_parameter_stat_samples.txt");
    for (i = 2; i < expected_num_cols; i++) {
        append_i_array(stat_indices, i);
    }
    summarize_stat_samples(paths, line_buffer, stat_indices, ss_array, means,
            std_devs, num_to_sample, expected_num_cols, paths_used);
    ck_assert_int_eq(means->length, stat_indices->length);
    ck_assert_int_eq(std_devs->length, stat_indices->length);
    ck_assert_int_eq(ss_array->a[0]->n, num_to_sample);
    ck_assert_msg((s_arrays_equal(paths, paths_used) != 0),
            "paths used do not match input paths");
    append_d_array(exp_means, 0.2);
    append_d_array(exp_means, 0.22);
    append_d_array(exp_means, 2.0);
    append_d_array(exp_means, 3.33333333);
    append_d_array(exp_std_devs, 0.1);
    append_d_array(exp_std_devs, 0.01);
    append_d_array(exp_std_devs, 1.0);
    append_d_array(exp_std_devs, 1.154701);
    ck_assert_msg((d_arrays_equal(means, exp_means, 0.000001) != 0),
            "unexpected means");
    ck_assert_msg((d_arrays_equal(std_devs, exp_std_devs, 0.000001) != 0),
            "unexpected std deviations");
    free_s_array(paths);
    free_s_array(paths_used);
    free_c_array(line_buffer);
    free_i_array(stat_indices);
    free_sample_sum_array(ss_array);
    free_d_array(means);
    free_d_array(std_devs);
    free_d_array(exp_means);
    free_d_array(exp_std_devs);
}
END_TEST

START_TEST (test_summarize_stat_samples_p1_n3_c2) {
    int i;
    s_array * paths;
    s_array * paths_used;
    c_array * line_buffer;
    i_array * stat_indices;
    sample_sum_array * ss_array;
    d_array * means;
    d_array * std_devs;
    d_array * exp_means;
    d_array * exp_std_devs;
    int num_to_sample;
    int expected_num_cols;
    expected_num_cols = 6;
    num_to_sample = 3;
    paths = init_s_array(1);
    paths_used = init_s_array(1);
    line_buffer = init_c_array(1023);
    stat_indices = init_i_array(4);
    ss_array = init_sample_sum_array(2);
    means = init_d_array(1);
    std_devs = init_d_array(1);
    exp_means = init_d_array(1);
    exp_std_devs = init_d_array(1);
    append_s_array(paths, "data/test_parameter_stat_samples.txt");
    for (i = 3; i < expected_num_cols; i += 2) {
        append_i_array(stat_indices, i);
    }
    summarize_stat_samples(paths, line_buffer, stat_indices, ss_array, means,
            std_devs, num_to_sample, expected_num_cols, paths_used);
    ck_assert_int_eq(means->length, stat_indices->length);
    ck_assert_int_eq(std_devs->length, stat_indices->length);
    ck_assert_int_eq(ss_array->a[0]->n, num_to_sample);
    ck_assert_msg((s_arrays_equal(paths, paths_used) != 0),
            "paths used do not match input paths");
    append_d_array(exp_means, 0.22);
    append_d_array(exp_means, 3.33333333);
    append_d_array(exp_std_devs, 0.01);
    append_d_array(exp_std_devs, 1.154701);
    ck_assert_msg((d_arrays_equal(means, exp_means, 0.000001) != 0),
            "unexpected means");
    ck_assert_msg((d_arrays_equal(std_devs, exp_std_devs, 0.000001) != 0),
            "unexpected std deviations");
    free_s_array(paths);
    free_s_array(paths_used);
    free_c_array(line_buffer);
    free_i_array(stat_indices);
    free_sample_sum_array(ss_array);
    free_d_array(means);
    free_d_array(std_devs);
    free_d_array(exp_means);
    free_d_array(exp_std_devs);
}
END_TEST

START_TEST (test_summarize_stat_samples_p2_n3) {
    int i;
    s_array * paths;
    s_array * paths_used;
    c_array * line_buffer;
    i_array * stat_indices;
    sample_sum_array * ss_array;
    d_array * means;
    d_array * std_devs;
    d_array * exp_means;
    d_array * exp_std_devs;
    int num_to_sample;
    int expected_num_cols;
    expected_num_cols = 6;
    num_to_sample = 3;
    paths = init_s_array(1);
    paths_used = init_s_array(1);
    line_buffer = init_c_array(1023);
    stat_indices = init_i_array(4);
    ss_array = init_sample_sum_array(4);
    means = init_d_array(1);
    std_devs = init_d_array(1);
    exp_means = init_d_array(1);
    exp_std_devs = init_d_array(1);
    append_s_array(paths, "data/test_parameter_stat_samples.txt");
    append_s_array(paths, "data/test_parameter_stat_samples2.txt");
    for (i = 2; i < expected_num_cols; i++) {
        append_i_array(stat_indices, i);
    }
    summarize_stat_samples(paths, line_buffer, stat_indices, ss_array, means,
            std_devs, num_to_sample, expected_num_cols, paths_used);
    ck_assert_int_eq(means->length, stat_indices->length);
    ck_assert_int_eq(std_devs->length, stat_indices->length);
    ck_assert_int_eq(ss_array->a[0]->n, num_to_sample);
    ck_assert_int_eq(paths_used->length, 1);
    ck_assert_msg((*get_s_array(paths, 0) == *get_s_array(paths_used, 0)),
            "paths used do not match input paths");
    append_d_array(exp_means, 0.2);
    append_d_array(exp_means, 0.22);
    append_d_array(exp_means, 2.0);
    append_d_array(exp_means, 3.33333333);
    append_d_array(exp_std_devs, 0.1);
    append_d_array(exp_std_devs, 0.01);
    append_d_array(exp_std_devs, 1.0);
    append_d_array(exp_std_devs, 1.154701);
    ck_assert_msg((d_arrays_equal(means, exp_means, 0.000001) != 0),
            "unexpected means");
    ck_assert_msg((d_arrays_equal(std_devs, exp_std_devs, 0.000001) != 0),
            "unexpected std deviations");
    free_s_array(paths);
    free_s_array(paths_used);
    free_c_array(line_buffer);
    free_i_array(stat_indices);
    free_sample_sum_array(ss_array);
    free_d_array(means);
    free_d_array(std_devs);
    free_d_array(exp_means);
    free_d_array(exp_std_devs);
}
END_TEST

START_TEST (test_summarize_stat_samples_p2_n4) {
    int i;
    s_array * paths;
    s_array * paths_used;
    c_array * line_buffer;
    i_array * stat_indices;
    sample_sum_array * ss_array;
    d_array * means;
    d_array * std_devs;
    d_array * exp_means;
    d_array * exp_std_devs;
    int num_to_sample;
    int expected_num_cols;
    expected_num_cols = 6;
    num_to_sample = 4;
    paths = init_s_array(1);
    paths_used = init_s_array(1);
    line_buffer = init_c_array(1023);
    stat_indices = init_i_array(4);
    ss_array = init_sample_sum_array(4);
    means = init_d_array(1);
    std_devs = init_d_array(1);
    exp_means = init_d_array(1);
    exp_std_devs = init_d_array(1);
    append_s_array(paths, "data/test_parameter_stat_samples.txt");
    append_s_array(paths, "data/test_parameter_stat_samples2.txt");
    for (i = 2; i < expected_num_cols; i++) {
        append_i_array(stat_indices, i);
    }
    summarize_stat_samples(paths, line_buffer, stat_indices, ss_array, means,
            std_devs, num_to_sample, expected_num_cols, paths_used);
    ck_assert_int_eq(means->length, stat_indices->length);
    ck_assert_int_eq(std_devs->length, stat_indices->length);
    ck_assert_int_eq(ss_array->a[0]->n, num_to_sample);
    ck_assert_int_eq(paths_used->length, 1);
    ck_assert_msg((*get_s_array(paths, 0) == *get_s_array(paths_used, 0)),
            "paths used do not match input paths");
    append_d_array(exp_means, 0.25);
    append_d_array(exp_means, 0.2225);
    append_d_array(exp_means, 2.5);
    append_d_array(exp_means, 3.0);
    append_d_array(exp_std_devs, 0.1290994);
    append_d_array(exp_std_devs, 0.009574271);
    append_d_array(exp_std_devs, 1.290994);
    append_d_array(exp_std_devs, 1.154701);
    ck_assert_msg((d_arrays_equal(means, exp_means, 0.000001) != 0),
            "unexpected means");
    ck_assert_msg((d_arrays_equal(std_devs, exp_std_devs, 0.000001) != 0),
            "unexpected std deviations");
    free_s_array(paths);
    free_s_array(paths_used);
    free_c_array(line_buffer);
    free_i_array(stat_indices);
    free_sample_sum_array(ss_array);
    free_d_array(means);
    free_d_array(std_devs);
    free_d_array(exp_means);
    free_d_array(exp_std_devs);
}
END_TEST

START_TEST (test_summarize_stat_samples_p2_n5) {
    int i;
    s_array * paths;
    s_array * paths_used;
    c_array * line_buffer;
    i_array * stat_indices;
    sample_sum_array * ss_array;
    d_array * means;
    d_array * std_devs;
    d_array * exp_means;
    d_array * exp_std_devs;
    int num_to_sample;
    int expected_num_cols;
    expected_num_cols = 6;
    num_to_sample = 5;
    paths = init_s_array(1);
    paths_used = init_s_array(1);
    line_buffer = init_c_array(1023);
    stat_indices = init_i_array(4);
    ss_array = init_sample_sum_array(4);
    means = init_d_array(1);
    std_devs = init_d_array(1);
    exp_means = init_d_array(1);
    exp_std_devs = init_d_array(1);
    append_s_array(paths, "data/test_parameter_stat_samples.txt");
    append_s_array(paths, "data/test_parameter_stat_samples2.txt");
    for (i = 2; i < expected_num_cols; i++) {
        append_i_array(stat_indices, i);
    }
    summarize_stat_samples(paths, line_buffer, stat_indices, ss_array, means,
            std_devs, num_to_sample, expected_num_cols, paths_used);
    ck_assert_int_eq(means->length, stat_indices->length);
    ck_assert_int_eq(std_devs->length, stat_indices->length);
    ck_assert_int_eq(ss_array->a[0]->n, num_to_sample);
    ck_assert_msg((s_arrays_equal(paths, paths_used) != 0),
            "paths used do not match input paths");
    append_d_array(exp_means, 0.222);
    append_d_array(exp_means, 0.218);
    append_d_array(exp_means, 4.2);
    append_d_array(exp_means, 3.6);
    append_d_array(exp_std_devs, 0.1281405);
    append_d_array(exp_std_devs, 0.0130384);
    append_d_array(exp_std_devs, 3.962323);
    append_d_array(exp_std_devs, 1.67332);
    ck_assert_msg((d_arrays_equal(means, exp_means, 0.000001) != 0),
            "unexpected means");
    ck_assert_msg((d_arrays_equal(std_devs, exp_std_devs, 0.000001) != 0),
            "unexpected std deviations");
    free_s_array(paths);
    free_s_array(paths_used);
    free_c_array(line_buffer);
    free_i_array(stat_indices);
    free_sample_sum_array(ss_array);
    free_d_array(means);
    free_d_array(std_devs);
    free_d_array(exp_means);
    free_d_array(exp_std_devs);
}
END_TEST

START_TEST (test_summarize_stat_samples_p1_n4_missing_cell) {
    int i;
    s_array * paths;
    s_array * paths_used;
    c_array * line_buffer;
    i_array * stat_indices;
    sample_sum_array * ss_array;
    d_array * means;
    d_array * std_devs;
    d_array * exp_means;
    d_array * exp_std_devs;
    int num_to_sample;
    int expected_num_cols;
    expected_num_cols = 6;
    num_to_sample = 4;
    paths = init_s_array(1);
    paths_used = init_s_array(1);
    line_buffer = init_c_array(1023);
    stat_indices = init_i_array(4);
    ss_array = init_sample_sum_array(4);
    means = init_d_array(1);
    std_devs = init_d_array(1);
    exp_means = init_d_array(1);
    exp_std_devs = init_d_array(1);
    append_s_array(paths, "data/test_parameter_stat_samples.missing_cell.txt");
    for (i = 2; i < expected_num_cols; i++) {
        append_i_array(stat_indices, i);
    }
    summarize_stat_samples(paths, line_buffer, stat_indices, ss_array, means,
            std_devs, num_to_sample, expected_num_cols, paths_used); // exit(1)
    free_s_array(paths);
    free_s_array(paths_used);
    free_c_array(line_buffer);
    free_i_array(stat_indices);
    free_sample_sum_array(ss_array);
    free_d_array(means);
    free_d_array(std_devs);
    free_d_array(exp_means);
    free_d_array(exp_std_devs);
}
END_TEST

START_TEST (test_reject_p1_n1_c4) {
    int i;
    s_array * paths;
    c_array * line_buffer;
    i_array * stat_indices;
    d_array * obs_stats;
    d_array * means;
    d_array * std_devs;
    s_array * header;
    sample_array * samples;
    int num_to_retain;
    num_to_retain = 1;
    paths = init_s_array(1);
    line_buffer = init_c_array(1023);
    stat_indices = init_i_array(4);
    means = init_d_array(1);
    std_devs = init_d_array(1);
    obs_stats = init_d_array(1);
    header = init_s_array(1);
    append_s_array(paths, "data/test_parameter_stat_samples.txt");
    parse_header(get_s_array(paths, 0), line_buffer, header);
    for (i = 2; i < header->length; i++) {
        append_i_array(stat_indices, i);
    }
    append_d_array(obs_stats, 0.1);
    append_d_array(obs_stats, 0.21);
    append_d_array(obs_stats, 1.0);
    append_d_array(obs_stats, 2.0);
    append_d_array(means, 0.25);
    append_d_array(means, 0.2225);
    append_d_array(means, 2.5);
    append_d_array(means, 3.0);
    append_d_array(std_devs, 0.1290994);
    append_d_array(std_devs, 0.009574271);
    append_d_array(std_devs, 1.290994);
    append_d_array(std_devs, 1.154701);
    standardize_vector(obs_stats, means, std_devs);
    samples = reject(paths, line_buffer, stat_indices, obs_stats, means,
            std_devs, num_to_retain, header);
    ck_assert_int_eq(samples->length, num_to_retain);
    ck_assert_int_eq(samples->capacity, num_to_retain);
    ck_assert_msg((almost_equal(samples->a[0]->distance, 0.0, 0.000001)),
            "euclidean distance was %lf, expected %lf",
            samples->a[0]->distance, 0.0);
    free_s_array(paths);
    free_c_array(line_buffer);
    free_i_array(stat_indices);
    free_d_array(obs_stats);
    free_sample_array(samples);
    free_d_array(means);
    free_d_array(std_devs);
    free_s_array(header);
}
END_TEST

START_TEST (test_reject_p1_n1_c2) {
    int i;
    s_array * paths;
    c_array * line_buffer;
    i_array * stat_indices;
    d_array * obs_stats;
    d_array * means;
    d_array * std_devs;
    s_array * header;
    sample_array * samples;
    int num_to_retain;
    num_to_retain = 1;
    paths = init_s_array(1);
    line_buffer = init_c_array(1023);
    stat_indices = init_i_array(4);
    means = init_d_array(1);
    std_devs = init_d_array(1);
    obs_stats = init_d_array(1);
    header = init_s_array(1);
    append_s_array(paths, "data/test_parameter_stat_samples.txt");
    parse_header(get_s_array(paths, 0), line_buffer, header);
    for (i = 3; i < header->length; i += 2) {
        append_i_array(stat_indices, i);
    }
    append_d_array(obs_stats, 0.21);
    append_d_array(obs_stats, 2.0);
    append_d_array(means, 0.2225);
    append_d_array(means, 3.0);
    append_d_array(std_devs, 0.009574271);
    append_d_array(std_devs, 1.154701);
    standardize_vector(obs_stats, means, std_devs);
    samples = reject(paths, line_buffer, stat_indices, obs_stats, means,
            std_devs, num_to_retain, header);
    ck_assert_int_eq(samples->length, num_to_retain);
    ck_assert_int_eq(samples->capacity, num_to_retain);
    ck_assert_msg((almost_equal(samples->a[0]->distance, 0.0, 0.000001)),
            "euclidean distance was %lf, expected %lf",
            samples->a[0]->distance, 0.0);
    free_s_array(paths);
    free_c_array(line_buffer);
    free_i_array(stat_indices);
    free_d_array(obs_stats);
    free_sample_array(samples);
    free_d_array(means);
    free_d_array(std_devs);
    free_s_array(header);
}
END_TEST

START_TEST (test_reject_p2_n1_c4) {
    int i;
    s_array * paths;
    c_array * line_buffer;
    i_array * stat_indices;
    d_array * obs_stats;
    d_array * means;
    d_array * std_devs;
    s_array * header;
    sample_array * samples;
    int num_to_retain;
    num_to_retain = 1;
    paths = init_s_array(1);
    line_buffer = init_c_array(1023);
    stat_indices = init_i_array(4);
    means = init_d_array(1);
    std_devs = init_d_array(1);
    obs_stats = init_d_array(1);
    header = init_s_array(1);
    append_s_array(paths, "data/test_parameter_stat_samples.txt");
    append_s_array(paths, "data/test_parameter_stat_samples2.txt");
    parse_header(get_s_array(paths, 0), line_buffer, header);
    for (i = 2; i < header->length; i++) {
        append_i_array(stat_indices, i);
    }
    append_d_array(obs_stats, 0.1);
    append_d_array(obs_stats, 0.21);
    append_d_array(obs_stats, 1.0);
    append_d_array(obs_stats, 2.0);
    append_d_array(means, 0.25);
    append_d_array(means, 0.2225);
    append_d_array(means, 2.5);
    append_d_array(means, 3.0);
    append_d_array(std_devs, 0.1290994);
    append_d_array(std_devs, 0.009574271);
    append_d_array(std_devs, 1.290994);
    append_d_array(std_devs, 1.154701);
    standardize_vector(obs_stats, means, std_devs);
    samples = reject(paths, line_buffer, stat_indices, obs_stats, means,
            std_devs, num_to_retain, header);
    ck_assert_int_eq(samples->length, num_to_retain);
    ck_assert_int_eq(samples->capacity, num_to_retain);
    ck_assert_msg((almost_equal(samples->a[0]->distance, 0.0, 0.000001)),
            "euclidean distance was %lf, expected %lf",
            samples->a[0]->distance, 0.0);
    free_s_array(paths);
    free_c_array(line_buffer);
    free_i_array(stat_indices);
    free_d_array(obs_stats);
    free_sample_array(samples);
    free_d_array(means);
    free_d_array(std_devs);
    free_s_array(header);
}
END_TEST

START_TEST (test_reject_p1_n2_c4) {
    int i;
    s_array * paths;
    c_array * line_buffer;
    i_array * stat_indices;
    d_array * obs_stats;
    d_array * means;
    d_array * std_devs;
    s_array * header;
    sample_array * samples;
    int num_to_retain;
    num_to_retain = 2;
    paths = init_s_array(1);
    line_buffer = init_c_array(1023);
    stat_indices = init_i_array(4);
    means = init_d_array(1);
    std_devs = init_d_array(1);
    obs_stats = init_d_array(1);
    header = init_s_array(1);
    append_s_array(paths, "data/test_parameter_stat_samples3.txt");
    parse_header(get_s_array(paths, 0), line_buffer, header);
    for (i = 2; i < header->length; i++) {
        append_i_array(stat_indices, i);
    }
    append_d_array(obs_stats, 0.1);
    append_d_array(obs_stats, 0.21);
    append_d_array(obs_stats, 1.0);
    append_d_array(obs_stats, 2.0);
    append_d_array(means, 0.25);
    append_d_array(means, 0.2225);
    append_d_array(means, 2.5);
    append_d_array(means, 3.0);
    append_d_array(std_devs, 0.1290994);
    append_d_array(std_devs, 0.009574271);
    append_d_array(std_devs, 1.290994);
    append_d_array(std_devs, 1.154701);
    standardize_vector(obs_stats, means, std_devs);
    samples = reject(paths, line_buffer, stat_indices, obs_stats, means,
            std_devs, num_to_retain, header);
    ck_assert_int_eq(samples->length, num_to_retain);
    ck_assert_int_eq(samples->capacity, num_to_retain);
    ck_assert_msg((almost_equal(samples->a[0]->distance, 0.0, 0.000001)),
            "euclidean distance was %lf, expected %lf",
            samples->a[0]->distance, 0.0);
    ck_assert_msg((almost_equal(samples->a[1]->distance, 0.130035, 0.000001)),
            "euclidean distance was %lf, expected %lf",
            samples->a[1]->distance, 0.130035);
    free_s_array(paths);
    free_c_array(line_buffer);
    free_i_array(stat_indices);
    free_d_array(obs_stats);
    free_sample_array(samples);
    free_d_array(means);
    free_d_array(std_devs);
    free_s_array(header);
}
END_TEST

START_TEST (test_reject_p2_n2_c4) {
    int i;
    s_array * paths;
    c_array * line_buffer;
    i_array * stat_indices;
    d_array * obs_stats;
    d_array * means;
    d_array * std_devs;
    s_array * header;
    sample_array * samples;
    int num_to_retain;
    num_to_retain = 2;
    paths = init_s_array(1);
    line_buffer = init_c_array(1023);
    stat_indices = init_i_array(4);
    means = init_d_array(1);
    std_devs = init_d_array(1);
    obs_stats = init_d_array(1);
    header = init_s_array(1);
    append_s_array(paths, "data/test_parameter_stat_samples4.txt");
    append_s_array(paths, "data/test_parameter_stat_samples.txt");
    parse_header(get_s_array(paths, 0), line_buffer, header);
    for (i = 2; i < header->length; i++) {
        append_i_array(stat_indices, i);
    }
    append_d_array(obs_stats, 0.1);
    append_d_array(obs_stats, 0.21);
    append_d_array(obs_stats, 1.0);
    append_d_array(obs_stats, 2.0);
    append_d_array(means, 0.25);
    append_d_array(means, 0.2225);
    append_d_array(means, 2.5);
    append_d_array(means, 3.0);
    append_d_array(std_devs, 0.1290994);
    append_d_array(std_devs, 0.009574271);
    append_d_array(std_devs, 1.290994);
    append_d_array(std_devs, 1.154701);
    standardize_vector(obs_stats, means, std_devs);
    samples = reject(paths, line_buffer, stat_indices, obs_stats, means,
            std_devs, num_to_retain, header);
    ck_assert_int_eq(samples->length, num_to_retain);
    ck_assert_int_eq(samples->capacity, num_to_retain);
    ck_assert_msg((almost_equal(samples->a[0]->distance, 0.0, 0.000001)),
            "euclidean distance was %lf, expected %lf",
            samples->a[0]->distance, 0.0);
    ck_assert_msg((almost_equal(samples->a[1]->distance, 0.130035, 0.000001)),
            "euclidean distance was %lf, expected %lf",
            samples->a[1]->distance, 0.130035);
    free_s_array(paths);
    free_c_array(line_buffer);
    free_i_array(stat_indices);
    free_d_array(obs_stats);
    free_sample_array(samples);
    free_d_array(means);
    free_d_array(std_devs);
    free_s_array(header);
}
END_TEST

START_TEST (test_reject_p2_n3_c4) {
    int i;
    s_array * paths;
    c_array * line_buffer;
    i_array * stat_indices;
    d_array * obs_stats;
    d_array * means;
    d_array * std_devs;
    s_array * header;
    sample_array * samples;
    int num_to_retain;
    num_to_retain = 3;
    paths = init_s_array(1);
    line_buffer = init_c_array(1023);
    stat_indices = init_i_array(4);
    means = init_d_array(1);
    std_devs = init_d_array(1);
    obs_stats = init_d_array(1);
    header = init_s_array(1);
    append_s_array(paths, "data/test_parameter_stat_samples4.txt");
    append_s_array(paths, "data/test_parameter_stat_samples.txt");
    parse_header(get_s_array(paths, 0), line_buffer, header);
    for (i = 2; i < header->length; i++) {
        append_i_array(stat_indices, i);
    }
    append_d_array(obs_stats, 0.1);
    append_d_array(obs_stats, 0.21);
    append_d_array(obs_stats, 1.0);
    append_d_array(obs_stats, 2.0);
    append_d_array(means, 0.25);
    append_d_array(means, 0.2225);
    append_d_array(means, 2.5);
    append_d_array(means, 3.0);
    append_d_array(std_devs, 0.1290994);
    append_d_array(std_devs, 0.009574271);
    append_d_array(std_devs, 1.290994);
    append_d_array(std_devs, 1.154701);
    standardize_vector(obs_stats, means, std_devs);
    samples = reject(paths, line_buffer, stat_indices, obs_stats, means,
            std_devs, num_to_retain, header);
    ck_assert_int_eq(samples->length, num_to_retain);
    ck_assert_int_eq(samples->capacity, num_to_retain);
    ck_assert_msg((almost_equal(samples->a[0]->distance, 0.0, 0.000001)),
            "euclidean distance was %lf, expected %lf",
            samples->a[0]->distance, 0.0);
    ck_assert_msg((almost_equal(samples->a[1]->distance, 0.130035, 0.000001)),
            "euclidean distance was %lf, expected %lf",
            samples->a[1]->distance, 0.130035);
    ck_assert_msg((almost_equal(samples->a[2]->distance, 0.1433004, 0.000001)),
            "euclidean distance was %lf, expected %lf",
            samples->a[2]->distance, 0.1433004);
    free_s_array(paths);
    free_c_array(line_buffer);
    free_i_array(stat_indices);
    free_d_array(obs_stats);
    free_sample_array(samples);
    free_d_array(means);
    free_d_array(std_devs);
    free_s_array(header);
}
END_TEST

Suite * eureject_suite(void) {
    Suite * s = suite_create("eureject");

    TCase * tc_eureject_config = tcase_create("config_test_case");
    tcase_add_test(tc_eureject_config, test_init_free_config);
    suite_add_tcase(s, tc_eureject_config);

    TCase * tc_summarize_stat_samples = tcase_create(
            "summarize_stat_samples_test_case");
    tcase_add_test(tc_summarize_stat_samples,
            test_summarize_stat_samples_p1_n4);
    tcase_add_test(tc_summarize_stat_samples,
            test_summarize_stat_samples_p1_n3);
    tcase_add_test(tc_summarize_stat_samples,
            test_summarize_stat_samples_p1_n3_c2);
    tcase_add_test(tc_summarize_stat_samples,
            test_summarize_stat_samples_p2_n3);
    tcase_add_test(tc_summarize_stat_samples,
            test_summarize_stat_samples_p2_n4);
    tcase_add_test(tc_summarize_stat_samples,
            test_summarize_stat_samples_p2_n5);
    tcase_add_exit_test(tc_summarize_stat_samples,
            test_summarize_stat_samples_p1_n4_missing_cell, 1);
    suite_add_tcase(s, tc_summarize_stat_samples);

    TCase * tc_reject = tcase_create("reject_test_case");
    tcase_add_test(tc_reject, test_reject_p1_n1_c4);
    tcase_add_test(tc_reject, test_reject_p1_n1_c2);
    tcase_add_test(tc_reject, test_reject_p2_n1_c4);
    tcase_add_test(tc_reject, test_reject_p1_n2_c4);
    tcase_add_test(tc_reject, test_reject_p2_n2_c4);
    tcase_add_test(tc_reject, test_reject_p2_n3_c4);
    suite_add_tcase(s, tc_reject);

    return s;
}

int main(void) {
    int number_failed;
    Suite * s = eureject_suite();
    SRunner * sr = srunner_create(s);
    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

