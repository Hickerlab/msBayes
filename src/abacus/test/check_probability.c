#include <stdlib.h>
#include <check.h>
#include <signal.h>
#include "../src/probability.c"
#include "test_rng.h"

START_TEST (test_draw_gamma_or_uniform_gamma) {
    double e, shape, scale, d, sum, mn;
    int i, reps;
    gsl_rng * rng;
    rng = get_rng(0);
    reps = 100000;
    shape = 1.0;
    scale = 10.0;
    for (i = 0; i < reps; i++) {
        d = draw_gamma_or_uniform(rng, shape, scale);
        ck_assert_msg((d > 0.0), "gamma draw %lf is less than zero", d);
        sum += d;
    }
    e = 0.1;
    mn = sum / reps;
    ck_assert_msg(almost_equal(mn, 10.0, e), "mean of gamma draws "
            "is %lf, expecting 10.0", mn);
    free_rng(rng);
}
END_TEST

START_TEST (test_draw_gamma_or_uniform_zero) {
    double e, shape, scale, d, sum, mn;
    int i, reps;
    gsl_rng * rng;
    rng = get_rng(0);
    reps = 100000;
    shape = 0.0;
    scale = 10.0;
    for (i = 0; i < reps; i++) {
        d = draw_gamma_or_uniform(rng, shape, scale);
        ck_assert_msg((d > 0.0), "uniform draw %lf is less than zero", d);
        sum += d;
    }
    e = 0.1;
    mn = sum / reps;
    ck_assert_msg(almost_equal(mn, 5.0, e), "mean of uniform draws "
            "is %lf, expecting 5.0", mn);
    free_rng(rng);
}
END_TEST

START_TEST (test_draw_gamma_or_uniform_zero_neg) {
    double e, shape, scale, d, sum, mn;
    int i, reps;
    gsl_rng * rng;
    rng = get_rng(0);
    reps = 100000;
    shape = 0.0;
    scale = -10.0;
    for (i = 0; i < reps; i++) {
        d = draw_gamma_or_uniform(rng, shape, scale);
        ck_assert_msg((d > 0.0), "uniform draw %lf is less than zero", d);
        sum += d;
    }
    e = 0.1;
    mn = sum / reps;
    ck_assert_msg(almost_equal(mn, 5.0, e), "mean of uniform draws "
            "is %lf, expecting 5.0", mn);
    free_rng(rng);
}
END_TEST

START_TEST (test_draw_gamma_or_uniform_pos_neg) {
    double e, shape, scale, d, sum, mn;
    int i, reps;
    gsl_rng * rng;
    rng = get_rng(0);
    reps = 100000;
    shape = 20.0;
    scale = -10.0;
    for (i = 0; i < reps; i++) {
        d = draw_gamma_or_uniform(rng, shape, scale);
        ck_assert_msg((d > 0.0), "uniform draw %lf is less than zero", d);
        sum += d;
    }
    e = 0.1;
    mn = sum / reps;
    ck_assert_msg(almost_equal(mn, 15.0, e), "mean of uniform draws "
            "is %lf, expecting 15.0", mn);
    free_rng(rng);
}
END_TEST

START_TEST (test_draw_gamma_or_uniform_neg_neg) {
    double e, shape, scale, d, sum, mn;
    int i, reps;
    gsl_rng * rng;
    rng = get_rng(0);
    reps = 100000;
    shape = -2.0;
    scale = -3.0;
    for (i = 0; i < reps; i++) {
        d = draw_gamma_or_uniform(rng, shape, scale);
        ck_assert_msg((d > 0.0), "uniform draw %lf is less than zero", d);
        sum += d;
    }
    e = 0.1;
    mn = sum / reps;
    ck_assert_msg(almost_equal(mn, 2.5, e), "mean of uniform draws "
            "is %lf, expecting 2.5", mn);
    free_rng(rng);
}
END_TEST

START_TEST (test_get_gamma_or_uniform_mean_pos_pos) {
    double e, shape, scale, mn;
    e = 0.0000001;
    shape = 4.0;
    scale = 5.0;
    mn = get_gamma_or_uniform_mean(shape, scale);
    ck_assert_msg(almost_equal(mn, (shape * scale), e), "mean is %lf, expecting "
            "%lf", mn, (shape * scale));
}
END_TEST

START_TEST (test_get_gamma_or_uniform_mean_zero_pos) {
    double e, shape, scale, mn, exp;
    e = 0.0000001;
    shape = 4.0;
    scale = 0.0;
    mn = get_gamma_or_uniform_mean(shape, scale);
    exp = (shape + scale) / 2;
    ck_assert_msg(almost_equal(mn, exp, e), "mean is %lf, expecting "
            "%lf", mn, exp);
}
END_TEST

START_TEST (test_get_gamma_or_uniform_mean_neg_pos) {
    double e, shape, scale, mn, exp;
    e = 0.0000001;
    shape = 5.0;
    scale = -4.0;
    mn = get_gamma_or_uniform_mean(shape, scale);
    exp = 4.5;
    ck_assert_msg(almost_equal(mn, exp, e), "mean is %lf, expecting "
            "%lf", mn, exp);
}
END_TEST

START_TEST (test_get_gamma_or_uniform_mean_neg_neg) {
    double e, shape, scale, mn, exp;
    e = 0.0000001;
    shape = -4.0;
    scale = -5.0;
    mn = get_gamma_or_uniform_mean(shape, scale);
    exp = 4.5;
    ck_assert_msg(almost_equal(mn, exp, e), "mean is %lf, expecting "
            "%lf", mn, exp);
}
END_TEST
    

Suite * probability_suite(void) {
    Suite * s = suite_create("probability");

    TCase * tc_draw_gamma_or_uniform = tcase_create(
            "draw_gamma_or_uniform_test_case");
    tcase_add_test(tc_draw_gamma_or_uniform, test_draw_gamma_or_uniform_gamma);
    tcase_add_test(tc_draw_gamma_or_uniform, test_draw_gamma_or_uniform_zero);
    tcase_add_test(tc_draw_gamma_or_uniform,
            test_draw_gamma_or_uniform_zero_neg);
    tcase_add_test(tc_draw_gamma_or_uniform,
            test_draw_gamma_or_uniform_pos_neg);
    tcase_add_test(tc_draw_gamma_or_uniform,
            test_draw_gamma_or_uniform_neg_neg);
    suite_add_tcase(s, tc_draw_gamma_or_uniform);

    TCase * tc_get_gamma_or_uniform_mean = tcase_create(
            "get_gamma_or_uniform_mean_test_case");
    tcase_add_test(tc_get_gamma_or_uniform_mean,
            test_get_gamma_or_uniform_mean_pos_pos);
    tcase_add_test(tc_get_gamma_or_uniform_mean,
            test_get_gamma_or_uniform_mean_zero_pos);
    tcase_add_test(tc_get_gamma_or_uniform_mean,
            test_get_gamma_or_uniform_mean_neg_pos);
    tcase_add_test(tc_get_gamma_or_uniform_mean,
            test_get_gamma_or_uniform_mean_neg_neg);
    suite_add_tcase(s, tc_get_gamma_or_uniform_mean);

    return s;
}

int main(void) {
    int number_failed;
    Suite * s = probability_suite();
    SRunner * sr = srunner_create(s);
    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

