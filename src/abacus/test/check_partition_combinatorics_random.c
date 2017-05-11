#include <stdlib.h>
#include <check.h>
#include <signal.h>
#include "../src/partition_combinatorics_random.c"
#include "test_rng.h"


START_TEST (test_draw_int_partition_category_n1) {
    gsl_rng * rng;
    rng = get_rng(0);
    int i, n, ret, target, sum, reps;
    double e, p, ep;
    e = 0.000001;
    n = 1;
    target = 1;
    reps = 100;
    sum = 0;
    for (i = 0; i < reps; i++) {
        ret = draw_int_partition_category(rng, n);
        if (ret == target) sum++;
    }
    ep = 1.0;
    p = sum / (double)reps;
    ck_assert_msg((almost_equal(p, ep, e) != 0),
            "target freq was %lf, expecting %lf", p, ep);
    free_rng(rng);
}
END_TEST

START_TEST (test_draw_int_partition_category_n2) {
    gsl_rng * rng;
    rng = get_rng(0);
    int i, n, ret, target, sum, reps;
    double e, p, ep;
    e = 0.01;
    n = 2;
    target = 1;
    reps = 100000;
    sum = 0;
    for (i = 0; i < reps; i++) {
        ret = draw_int_partition_category(rng, n);
        if (ret == target) sum++;
    }
    ep = 0.5;
    p = sum / (double)reps;
    ck_assert_msg((almost_equal(p, ep, e) != 0),
            "target freq was %lf, expecting %lf", p, ep);
    free_rng(rng);
}
END_TEST

START_TEST (test_draw_int_partition_category_n3) {
    gsl_rng * rng;
    rng = get_rng(0);
    int i, n, ret, target, sum, reps;
    double e, p, ep;
    e = 0.01;
    n = 3;
    target = 1;
    reps = 100000;
    sum = 0;
    for (i = 0; i < reps; i++) {
        ret = draw_int_partition_category(rng, n);
        if (ret == target) sum++;
    }
    ep = 0.3333333;
    p = sum / (double)reps;
    ck_assert_msg((almost_equal(p, ep, e) != 0),
            "target freq was %lf, expecting %lf", p, ep);
    free_rng(rng);
}
END_TEST

START_TEST (test_draw_int_partition_category_n7) {
    gsl_rng * rng;
    rng = get_rng(0);
    int i, n, ret, target, sum, reps;
    double e, p, ep;
    e = 0.01;
    n = 7;
    target = 3;
    reps = 100000;
    sum = 0;
    for (i = 0; i < reps; i++) {
        ret = draw_int_partition_category(rng, n);
        if (ret == target) sum++;
    }
    ep = 4 / (double)15;
    p = sum / (double)reps;
    ck_assert_msg((almost_equal(p, ep, e) != 0),
            "target freq was %lf, expecting %lf", p, ep);
    free_rng(rng);
}
END_TEST

START_TEST (test_draw_int_partition_category_n22) {
    gsl_rng * rng;
    rng = get_rng(0);
    int i, n, ret, target, sum, reps;
    double e, p, ep;
    e = 0.01;
    n = 22;
    target = 6;
    reps = 100000;
    sum = 0;
    for (i = 0; i < reps; i++) {
        ret = draw_int_partition_category(rng, n);
        if (ret == target) sum++;
    }
    ep = 136 / (double)1002;
    p = sum / (double)reps;
    ck_assert_msg((almost_equal(p, ep, e) != 0),
            "target freq was %lf, expecting %lf", p, ep);
    free_rng(rng);
}
END_TEST


START_TEST (test_dirichlet_process_draw_n10_a1) {
    gsl_rng * rng;
    rng = get_rng(0);
    i_array * elements;
    int i, n, idx1, idx2, sum, reps, ret;
    double e, p, ep, alpha;
    e = 0.01;
    n = 10;
    alpha = 1.0;
    reps = 100000;
    sum = 0;
    elements = init_i_array(n);
    for (i = 0; i < reps; i++) {
        ret = dirichlet_process_draw(rng, n, alpha, elements);
        ck_assert_msg(((ret > 0) && (ret <= n)), "invalid number of "
                "categories: %d", ret);
        idx1 = gsl_rng_uniform_int(rng, n);
        idx2 = gsl_rng_uniform_int(rng, n);
        while(idx1 == idx2) {
            idx2 = gsl_rng_uniform_int(rng, n);
        }
        if (get_i_array(elements, idx1) == get_i_array(elements, idx2)) {
            sum++;
        }
    }
    ep = 1 / (1.0 + alpha);
    p = sum / (double)reps;
    ck_assert_msg((almost_equal(p, ep, e) != 0),
            "target freq was %lf, expecting %lf", p, ep);
    free_i_array(elements);
    free_rng(rng);
}
END_TEST

START_TEST (test_dirichlet_process_draw_n10_a3) {
    gsl_rng * rng;
    rng = get_rng(0);
    i_array * elements;
    int i, n, idx1, idx2, sum, reps, ret;
    double e, p, ep, alpha;
    e = 0.01;
    n = 10;
    alpha = 3.0;
    reps = 100000;
    sum = 0;
    elements = init_i_array(n);
    for (i = 0; i < reps; i++) {
        ret = dirichlet_process_draw(rng, n, alpha, elements);
        ck_assert_msg(((ret > 0) && (ret <= n)), "invalid number of "
                "categories: %d", ret);
        idx1 = gsl_rng_uniform_int(rng, n);
        idx2 = gsl_rng_uniform_int(rng, n);
        while(idx1 == idx2) {
            idx2 = gsl_rng_uniform_int(rng, n);
        }
        if (get_i_array(elements, idx1) == get_i_array(elements, idx2)) {
            sum++;
        }
    }
    ep = 1 / (1.0 + alpha);
    p = sum / (double)reps;
    ck_assert_msg((almost_equal(p, ep, e) != 0),
            "target freq was %lf, expecting %lf", p, ep);
    free_i_array(elements);
    free_rng(rng);
}
END_TEST

START_TEST (test_dirichlet_process_draw_n5_a3) {
    gsl_rng * rng;
    rng = get_rng(0);
    i_array * elements;
    int i, n, idx1, idx2, sum, reps, ret;
    double e, p, ep, alpha;
    e = 0.01;
    n = 5;
    alpha = 3.0;
    reps = 100000;
    sum = 0;
    elements = init_i_array(n);
    for (i = 0; i < reps; i++) {
        ret = dirichlet_process_draw(rng, n, alpha, elements);
        ck_assert_msg(((ret > 0) && (ret <= n)), "invalid number of "
                "categories: %d", ret);
        idx1 = gsl_rng_uniform_int(rng, n);
        idx2 = gsl_rng_uniform_int(rng, n);
        while(idx1 == idx2) {
            idx2 = gsl_rng_uniform_int(rng, n);
        }
        if (get_i_array(elements, idx1) == get_i_array(elements, idx2)) {
            sum++;
        }
    }
    ep = 1 / (1.0 + alpha);
    p = sum / (double)reps;
    ck_assert_msg((almost_equal(p, ep, e) != 0),
            "target freq was %lf, expecting %lf", p, ep);
    free_i_array(elements);
    free_rng(rng);
}
END_TEST

START_TEST (test_dirichlet_process_draw_n10_a9) {
    gsl_rng * rng;
    rng = get_rng(0);
    i_array * elements;
    int i, n, idx1, idx2, sum, reps, ret;
    double e, p, ep, alpha;
    e = 0.01;
    n = 10;
    alpha = 9.0;
    reps = 100000;
    sum = 0;
    elements = init_i_array(n);
    for (i = 0; i < reps; i++) {
        ret = dirichlet_process_draw(rng, n, alpha, elements);
        ck_assert_msg(((ret > 0) && (ret <= n)), "invalid number of "
                "categories: %d", ret);
        idx1 = gsl_rng_uniform_int(rng, n);
        idx2 = gsl_rng_uniform_int(rng, n);
        while(idx1 == idx2) {
            idx2 = gsl_rng_uniform_int(rng, n);
        }
        if (get_i_array(elements, idx1) == get_i_array(elements, idx2)) {
            sum++;
        }
    }
    ep = 1 / (1.0 + alpha);
    p = sum / (double)reps;
    ck_assert_msg((almost_equal(p, ep, e) != 0),
            "target freq was %lf, expecting %lf", p, ep);
    free_i_array(elements);
    free_rng(rng);
}
END_TEST

Suite * partition_combinatorics_random_suite(void) {
    Suite * s = suite_create("partition_combinatorics_random");

    TCase * tc_draw_cat = tcase_create("draw_int_partition_category");
    tcase_add_test(tc_draw_cat,
            test_draw_int_partition_category_n1);
    tcase_add_test(tc_draw_cat,
            test_draw_int_partition_category_n2);
    tcase_add_test(tc_draw_cat,
            test_draw_int_partition_category_n3);
    tcase_add_test(tc_draw_cat,
            test_draw_int_partition_category_n7);
    tcase_add_test(tc_draw_cat,
            test_draw_int_partition_category_n22);
    suite_add_tcase(s, tc_draw_cat);

    TCase * tc_dirichlet_process_draw = tcase_create("dirichlet_process_draw");
    tcase_add_test(tc_dirichlet_process_draw,
            test_dirichlet_process_draw_n10_a1);
    tcase_add_test(tc_dirichlet_process_draw,
            test_dirichlet_process_draw_n10_a3);
    tcase_add_test(tc_dirichlet_process_draw,
            test_dirichlet_process_draw_n5_a3);
    tcase_add_test(tc_dirichlet_process_draw,
            test_dirichlet_process_draw_n10_a9);
    suite_add_tcase(s, tc_dirichlet_process_draw);

    return s;
}

int main(void) {
    int number_failed;
    Suite * s = partition_combinatorics_random_suite();
    SRunner * sr = srunner_create(s);
    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

