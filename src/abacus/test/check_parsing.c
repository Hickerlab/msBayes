#include <stdlib.h>
#include <check.h>
#include <signal.h>
#include "../src/parsing.c"
#include "test_utils.h"

START_TEST (test_parse_header) {
    char * path = "data/observed_stats.txt";
    int ret;
    c_array * line_buffer;
    s_array * header;
    s_array * expected_header;
    line_buffer = init_c_array(1023);
    header = init_s_array(1);
    expected_header = init_s_array(1);
    parse_header(path, line_buffer, header);
    ck_assert_msg((header->length == 4), "header length is %d, "
            "expecting %d", header->length, 4);
    append_s_array(expected_header, "stat.1");
    append_s_array(expected_header, "stat.2");
    append_s_array(expected_header, "stat.3");
    append_s_array(expected_header, "stat.4");
    ret = s_arrays_equal(header, expected_header);
    ck_assert_msg((ret != 0), "parsed header is incorrect");
    free_c_array(line_buffer);
    free_s_array(header);
    free_s_array(expected_header);
}
END_TEST

START_TEST (test_parse_observed_stats_file) {
    char * path = "data/observed_stats.txt";
    int ret;
    double e = 0.000001;
    c_array * line_buffer;
    s_array * header;
    s_array * expected_header;
    d_array * stats;
    d_array * expected_stats;
    line_buffer = init_c_array(1023);
    header = init_s_array(1);
    expected_header = init_s_array(1);
    append_s_array(expected_header, "stat.1");
    append_s_array(expected_header, "stat.2");
    append_s_array(expected_header, "stat.3");
    append_s_array(expected_header, "stat.4");
    stats = init_d_array(1);
    expected_stats = init_d_array(1);
    append_d_array(expected_stats, 0.4);
    append_d_array(expected_stats, 0.3);
    append_d_array(expected_stats, 0.3);
    append_d_array(expected_stats, 0.4);

    parse_observed_stats_file(path, line_buffer, header, stats);
    ck_assert_msg((header->length == 4), "header length is %d, "
            "expecting %d", header->length, 4);
    ck_assert_msg((stats->length == 4), "stats length is %d, "
            "expecting %d", header->length, 4);
    ret = s_arrays_equal(header, expected_header);
    ck_assert_msg((ret != 0), "parsed header is incorrect");
    ret = d_arrays_equal(stats, expected_stats, e);
    ck_assert_msg((ret != 0), "parsed stats are incorrect");
    free_c_array(line_buffer);
    free_s_array(header);
    free_s_array(expected_header);
    free_d_array(stats);
    free_d_array(expected_stats);
}
END_TEST

START_TEST (test_parse_observed_stats_file_extra_line) {
    char * path = "data/observed_stats_extra_line.txt";
    int ret;
    double e = 0.000001;
    c_array * line_buffer;
    s_array * header;
    s_array * expected_header;
    d_array * stats;
    d_array * expected_stats;
    line_buffer = init_c_array(1023);
    header = init_s_array(1);
    expected_header = init_s_array(1);
    append_s_array(expected_header, "stat.1");
    append_s_array(expected_header, "stat.2");
    append_s_array(expected_header, "stat.3");
    append_s_array(expected_header, "stat.4");
    stats = init_d_array(1);
    expected_stats = init_d_array(1);
    append_d_array(expected_stats, 0.1);
    append_d_array(expected_stats, 0.2);
    append_d_array(expected_stats, 0.3);
    append_d_array(expected_stats, 0.4);

    parse_observed_stats_file(path, line_buffer, header, stats);
    ck_assert_msg((header->length == 4), "header length is %d, "
            "expecting %d", header->length, 4);
    ck_assert_msg((stats->length == 4), "stats length is %d, "
            "expecting %d", stats->length, 4);
    ret = s_arrays_equal(header, expected_header);
    ck_assert_msg((ret != 0), "parsed header is incorrect");
    ret = d_arrays_equal(stats, expected_stats, e);
    ck_assert_msg((ret != 0), "parsed stats are incorrect");
    free_c_array(line_buffer);
    free_s_array(header);
    free_s_array(expected_header);
    free_d_array(stats);
    free_d_array(expected_stats);
}
END_TEST

START_TEST (test_parse_summary_file) {
    char * path = "data/observed_stats_extra_line.txt";
    int ret;
    double e = 0.000001;
    c_array * line_buffer;
    s_array * header;
    s_array * expected_header;
    d_array * means;
    d_array * expected_means;
    d_array * std_devs;
    d_array * expected_std_devs;
    i_array * sample_sizes;
    i_array * expected_sample_sizes;
    line_buffer = init_c_array(1023);
    header = init_s_array(1);
    expected_header = init_s_array(1);
    append_s_array(expected_header, "stat.1");
    append_s_array(expected_header, "stat.2");
    append_s_array(expected_header, "stat.3");
    append_s_array(expected_header, "stat.4");
    means = init_d_array(1);
    expected_means = init_d_array(1);
    append_d_array(expected_means, 0.1);
    append_d_array(expected_means, 0.2);
    append_d_array(expected_means, 0.3);
    append_d_array(expected_means, 0.4);
    std_devs = init_d_array(1);
    expected_std_devs = init_d_array(1);
    append_d_array(expected_std_devs, 0.4);
    append_d_array(expected_std_devs, 0.3);
    append_d_array(expected_std_devs, 0.3);
    append_d_array(expected_std_devs, 0.4);
    sample_sizes = init_i_array(1);
    expected_sample_sizes = init_i_array(1);
    append_i_array(expected_sample_sizes, 10);
    append_i_array(expected_sample_sizes, 10);
    append_i_array(expected_sample_sizes, 10);
    append_i_array(expected_sample_sizes, 10);

    parse_summary_file(path, line_buffer, header, means, std_devs,
            sample_sizes);
    ck_assert_msg((header->length == 4), "header length is %d, "
            "expecting %d", header->length, 4);
    ck_assert_msg((means->length == 4), "means length is %d, "
            "expecting %d", means->length, 4);
    ck_assert_msg((std_devs->length == 4), "std devs length is %d, "
            "expecting %d", std_devs->length, 4);
    ck_assert_msg((sample_sizes->length == 4), "sample sizes length is %d, "
            "expecting %d", sample_sizes->length, 4);
    ret = s_arrays_equal(header, expected_header);
    ck_assert_msg((ret != 0), "parsed header is incorrect");
    ret = d_arrays_equal(means, expected_means, e);
    ck_assert_msg((ret != 0), "parsed means are incorrect");
    ret = d_arrays_equal(std_devs, expected_std_devs, e);
    ck_assert_msg((ret != 0), "parsed std_devs are incorrect");
    ret = i_arrays_equal(sample_sizes, expected_sample_sizes);
    ck_assert_msg((ret != 0), "parsed sample_sizes are incorrect");
    free_c_array(line_buffer);
    free_s_array(header);
    free_s_array(expected_header);
    free_d_array(means);
    free_d_array(expected_means);
    free_d_array(std_devs);
    free_d_array(expected_std_devs);
    free_i_array(sample_sizes);
    free_i_array(expected_sample_sizes);
}
END_TEST

START_TEST (test_strcmp_i) {
    int r;
    char * a;
    char * b;
    a = "test";
    b = "test";
    r = strcmp_i(a, b);
    ck_assert_msg((r == 0), "strcmp_i(%s, %s) returned %d", a, b, r);

    a = "TeSt";
    b = "tEsT";
    r = strcmp_i(a, b);
    ck_assert_msg((r == 0), "strcmp_i(%s, %s) returned %d", a, b, r);

    a = "testt";
    b = "test";
    r = strcmp_i(a, b);
    ck_assert_msg((r != 0), "strcmp_i(%s, %s) returned %d", a, b, r);

    a = "test";
    b = "testt";
    r = strcmp_i(a, b);
    ck_assert_msg((r != 0), "strcmp_i(%s, %s) returned %d", a, b, r);

    a = "ttest";
    b = "test";
    r = strcmp_i(a, b);
    ck_assert_msg((r != 0), "strcmp_i(%s, %s) returned %d", a, b, r);

    a = "test";
    b = "ttest";
    r = strcmp_i(a, b);
    ck_assert_msg((r != 0), "strcmp_i(%s, %s) returned %d", a, b, r);
}
END_TEST

START_TEST (test_strip) {
    char * a = " testing 1 2 3 ";
    char * exp = "testing 1 2 3";
    char * b = strip(a);
    printf("<%s>\n", a);
    printf("<%s>\n", b);
    ck_assert_msg((strncmp(b, exp, 16) == 0),
            "strip(<%s>) returned <%s>", a, b);
    free(b);

    a = "testing 1 2 3        ";
    exp = "testing 1 2 3";
    b = strip(a);
    printf("<%s>\n", a);
    printf("<%s>\n", b);
    ck_assert_msg((strncmp(b, exp, 16) == 0),
            "strip(<%s>) returned <%s>", a, b);
    free(b);

    a = "      testing 1 2 3";
    exp = "testing 1 2 3";
    b = strip(a);
    printf("<%s>\n", a);
    printf("<%s>\n", b);
    ck_assert_msg((strncmp(b, exp, 16) == 0),
            "strip(<%s>) returned <%s>", a, b);
    free(b);

    a = "     testing 1 2 3       ";
    exp = "testing 1 2 3";
    b = strip(a);
    printf("<%s>\n", a);
    printf("<%s>\n", b);
    ck_assert_msg((strncmp(b, exp, 16) == 0),
            "strip(<%s>) returned <%s>", a, b);
    free(b);

    a = "testing 1 2 3";
    exp = "testing 1 2 3";
    b = strip(a);
    printf("<%s>\n", a);
    printf("<%s>\n", b);
    ck_assert_msg((strncmp(b, exp, 16) == 0),
            "strip(<%s>) returned <%s>", a, b);
    free(b);

    a = "test";
    exp = "test";
    b = strip(a);
    printf("<%s>\n", a);
    printf("<%s>\n", b);
    ck_assert_msg((strncmp(b, exp, 16) == 0),
            "strip(<%s>) returned <%s>", a, b);
    free(b);

    a = "";
    exp = "";
    b = strip(a);
    printf("<%s>\n", a);
    printf("<%s>\n", b);
    ck_assert_msg((strncmp(b, exp, 16) == 0),
            "strip(<%s>) returned <%s>", a, b);
    free(b);

    a = "       ";
    exp = "";
    b = strip(a);
    printf("<%s>\n", a);
    printf("<%s>\n", b);
    ck_assert_msg((strncmp(b, exp, 16) == 0),
            "strip(<%s>) returned <%s>", a, b);
    free(b);

    a = "testing 1 2 3\n";
    exp = "testing 1 2 3";
    b = strip(a);
    printf("<%s>\n", a);
    printf("<%s>\n", b);
    ck_assert_msg((strncmp(b, exp, 16) == 0),
            "strip(<%s>) returned <%s>", a, b);
    free(b);

    a = "  testing 1 2 3  \n";
    exp = "testing 1 2 3";
    b = strip(a);
    printf("<%s>\n", a);
    printf("<%s>\n", b);
    ck_assert_msg((strncmp(b, exp, 16) == 0),
            "strip(<%s>) returned <%s>", a, b);
    free(b);

    a = "  testing 1 2 3  \n\n";
    exp = "testing 1 2 3";
    b = strip(a);
    printf("<%s>\n", a);
    printf("<%s>\n", b);
    ck_assert_msg((strncmp(b, exp, 16) == 0),
            "strip(<%s>) returned <%s>", a, b);
    free(b);

    a = "  testing 1 2 3  \n  \n  ";
    exp = "testing 1 2 3";
    b = strip(a);
    printf("<%s>\n", a);
    printf("<%s>\n", b);
    ck_assert_msg((strncmp(b, exp, 16) == 0),
            "strip(<%s>) returned <%s>", a, b);
    free(b);
}
END_TEST

Suite * parsing_suite(void) {
    Suite * s = suite_create("parsing");

    TCase * tc_parse_header = tcase_create("parse_header_test_case");
    tcase_add_test(tc_parse_header, test_parse_header);
    suite_add_tcase(s, tc_parse_header);

    TCase * tc_parse_observed_stats_file = tcase_create(
            "parse_observed_stats_file_test_case");
    tcase_add_test(tc_parse_observed_stats_file,
            test_parse_observed_stats_file);
    tcase_add_test(tc_parse_observed_stats_file,
            test_parse_observed_stats_file_extra_line);
    suite_add_tcase(s, tc_parse_observed_stats_file);

    TCase * tc_parse_summary_file = tcase_create("parse_summary_file");
    tcase_add_test(tc_parse_summary_file, test_parse_summary_file);
    suite_add_tcase(s, tc_parse_summary_file);

    TCase * tc_strcmp_i = tcase_create("strcmp_i_test_case");
    tcase_add_test(tc_strcmp_i, test_strcmp_i);
    suite_add_tcase(s, tc_strcmp_i);

    TCase * tc_strip = tcase_create("strip_test_case");
    tcase_add_test(tc_strip, test_strip);
    suite_add_tcase(s, tc_strip);

    return s;
}

int main(void) {
    int number_failed;
    Suite * s = parsing_suite();
    SRunner * sr = srunner_create(s);
    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

