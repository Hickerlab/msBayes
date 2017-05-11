#include <stdlib.h>
#include <check.h>
#include <signal.h>
#include "../src/stats_utils.c"
#include "test_utils.h"

START_TEST (test_init_free_sample_sum) {
    sample_sum * ss;
    ss = init_sample_sum();
    ck_assert_msg((sizeof(ss->n) == sizeof(int)), "`n` "
            "attribute is not an int");
    ck_assert_msg((sizeof(ss->sum) == sizeof(double)), "`sum` "
            "attribute is not an double");
    ck_assert_msg((sizeof(ss->sum_of_squares) == sizeof(double)), 
            "`sum_of_squares` attribute is not an double");
    ck_assert_msg((ss->n == 0), "`sample_sum.n` did not initiate to 0");
    ck_assert_msg((ss->sum == 0.0), "`sample_sum.sum` did not initiate to 0");
    ck_assert_msg((ss->sum_of_squares == 0.0), "`sample_sum.sum_of_squares` "
            "did not initiate to 0");
    free_sample_sum(ss);
}
END_TEST

START_TEST (test_update_sample_sum) {
    double e = 0.000001;
    sample_sum * ss;
    ss = init_sample_sum();
    update_sample_sum(ss, 1.0);
    ck_assert_msg((ss->n == 1), "`sample_sum.n` is %d, expecting %d", ss->n, 1);
    ck_assert_msg(almost_equal(ss->sum, 1.0, e), "`sample_sum.sum` is %lf, "
            "expecting %lf", ss->sum, 1.0);
    ck_assert_msg(almost_equal(ss->sum_of_squares, 1.0, e),
            "`sample_sum.sum_of_squares is %lf, expecting %lf", ss->sum, 1.0);

    update_sample_sum(ss, 2.0);
    ck_assert_msg((ss->n == 2), "`sample_sum.n` is %d, expecting %d", ss->n, 2);
    ck_assert_msg(almost_equal(ss->sum, 3.0, e), "`sample_sum.sum` is %lf, "
            "expecting %lf", ss->sum, 3.0);
    ck_assert_msg(almost_equal(ss->sum_of_squares, 5.0, e),
            "`sample_sum.sum_of_squares is %lf, expecting %lf", ss->sum, 5.0);

    update_sample_sum(ss, 3.0);
    ck_assert_msg((ss->n == 3), "`sample_sum.n` is %d, expecting %d", ss->n, 2);
    ck_assert_msg(almost_equal(ss->sum, 6.0, e), "`sample_sum.sum` is %lf, "
            "expecting %lf", ss->sum, 6.0);
    ck_assert_msg(almost_equal(ss->sum_of_squares, 14.0, e),
            "`sample_sum.sum_of_squares is %lf, expecting %lf", ss->sum, 14.0);
    free_sample_sum(ss);
}
END_TEST

START_TEST (test_get_mean_n0) {
    sample_sum * ss;
    ss = init_sample_sum();
    get_mean(ss);
}
END_TEST

START_TEST (test_get_mean) {
    double e = 0.000001;
    sample_sum * ss;
    ss = init_sample_sum();
    update_sample_sum(ss, 100.0);
    ck_assert_msg(almost_equal(get_mean(ss), 100.0, e), "mean is %lf, "
            "expecting %lf", get_mean(ss), 100.0);
    update_sample_sum(ss, 101.0);
    ck_assert_msg(almost_equal(get_mean(ss), 100.5, e), "mean is %lf, "
            "expecting %lf", get_mean(ss), 100.5);
    update_sample_sum(ss, 117.0);
    ck_assert_msg(almost_equal(get_mean(ss), 106.0, e), "mean is %lf, "
            "expecting %lf", get_mean(ss), 106.0);
    update_sample_sum(ss, -500.0);
    ck_assert_msg(almost_equal(get_mean(ss), -45.5, e), "mean is %lf, "
            "expecting %lf", get_mean(ss), -45.5);
    free_sample_sum(ss);
}
END_TEST

START_TEST (test_get_sample_variance_n0) {
    sample_sum * ss;
    ss = init_sample_sum();
    get_sample_variance(ss);
}
END_TEST

START_TEST (test_get_sample_variance_n1) {
    sample_sum * ss;
    ss = init_sample_sum();
    update_sample_sum(ss, 100.0);
    get_sample_variance(ss);
}
END_TEST

START_TEST (test_get_sample_variance) {
    double e = 0.000001;
    sample_sum * ss;
    ss = init_sample_sum();
    update_sample_sum(ss, 100.0);
    update_sample_sum(ss, 101.0);
    ck_assert_msg(almost_equal(get_sample_variance(ss), 0.5, e),
            "variance is %lf, expecting %lf", get_sample_variance(ss), 0.5);
    update_sample_sum(ss, 117.0);
    ck_assert_msg(almost_equal(get_sample_variance(ss), 91.0, e),
            "variance is %lf, expecting %lf", get_sample_variance(ss), 91.0);
    update_sample_sum(ss, -500.0);
    ck_assert_msg(almost_equal(get_sample_variance(ss), 91869.66666667, e),
            "variance is %lf, expecting %lf", get_sample_variance(ss),
            91869.666666667);
    free_sample_sum(ss);
}
END_TEST

START_TEST (test_get_std_dev_n0) {
    sample_sum * ss;
    ss = init_sample_sum();
    get_std_dev(ss);
}
END_TEST

START_TEST (test_get_std_dev_n1) {
    sample_sum * ss;
    ss = init_sample_sum();
    update_sample_sum(ss, 100.0);
    get_std_dev(ss);
}
END_TEST

START_TEST (test_get_std_dev) {
    double e = 0.00001;
    sample_sum * ss;
    ss = init_sample_sum();
    update_sample_sum(ss, 100.0);
    update_sample_sum(ss, 101.0);
    ck_assert_msg(almost_equal(get_std_dev(ss), 0.7071068, e),
            "std deviation is %lf, expecting %lf", get_std_dev(ss), 0.7071068);
    update_sample_sum(ss, 117.0);
    ck_assert_msg(almost_equal(get_std_dev(ss), 9.539392, e),
            "std deviation is %lf, expecting %lf", get_std_dev(ss), 9.539392);
    update_sample_sum(ss, -500.0);
    ck_assert_msg(almost_equal(get_std_dev(ss), 303.1001, e),
            "std deviation is %lf, expecting %lf", get_std_dev(ss),
            303.1001);
    free_sample_sum(ss);
}
END_TEST

START_TEST (test_init_sample_sum_array_fail) {
    sample_sum_array * v;
    v = init_sample_sum_array(0); // SIGABRT
}
END_TEST

START_TEST (test_init_free_sample_sum_array) {
    int i;
    sample_sum_array * v;
    int size = 5;
    v = init_sample_sum_array(size);
    ck_assert_msg((sizeof(v->length) == sizeof(int)), "`length` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->a) == sizeof(long)), "Array attribute "
            "of `sample_sum_array` is not a pointer");
    ck_assert_msg((sizeof(*v->a) == sizeof(sample_sum *)), "Pointer to array"
            "of `sample_sum_array` is not pointing to array of `sample_sum`s");
    ck_assert_int_eq(v->length, size);
    for (i = 0; i < v->length; i++) {
        ck_assert_msg((v->a[i]->n == 0),
                "`n` of element %d did not initiate to 0", i);
        ck_assert_msg((v->a[i]->sum == 0.0),
                "`sum` of element %d did not initiate to 0", i);
        ck_assert_msg((v->a[i]->sum_of_squares == 0.0),
                "`sum_of_squares` of element %d "
                "did not initiate to 0", i);
    }
    free_sample_sum_array(v);
}
END_TEST

START_TEST (test_update_sample_sum_array) {
    int i;
    sample_sum_array * v;
    d_array * x;
    int size = 3;
    v = init_sample_sum_array(size);
    x = init_d_array(size);
    for (i = 0; i < size; i++) {
        append_d_array(x, ((i+1)/((double)100)));
    }
    update_sample_sum_array(v, x);
    for (i = 0; i < size; i++) {
        ck_assert_msg((v->a[i]->n == 1), "element %d has `n` of %d, should "
                "be %d", i, v->a[i]->n, 1);
        ck_assert_msg((v->a[i]->sum == get_d_array(x, i)), "element %d has "
                "`sum` of %lf, should be %lf",
                i, v->a[i]->sum, get_d_array(x, i));
        ck_assert_msg((v->a[i]->sum_of_squares == pow(get_d_array(x, i), 2)),
                "element %d has `sum_of_squares` of %lf, should be %lf",
                i, v->a[i]->sum_of_squares, pow(get_d_array(x, i), 2));
    }
    update_sample_sum_array(v, x);
    for (i = 0; i < size; i++) {
        ck_assert_msg((v->a[i]->n == 2), "element %d has `n` of %d, should "
                "be %d", i, v->a[i]->n, 2);
        ck_assert_msg((v->a[i]->sum == (get_d_array(x, i)*2)), "element %d has "
                "`sum` of %lf, should be %lf",
                i, v->a[i]->sum, (get_d_array(x, i)*2));
        ck_assert_msg((v->a[i]->sum_of_squares == (pow(get_d_array(x, i), 2)*2)),
                "element %d has `sum_of_squares` of %lf, should be %lf",
                i, v->a[i]->sum_of_squares, (pow(get_d_array(x, i), 2)*2));
    }
    free_sample_sum_array(v);
    free_d_array(x);
}
END_TEST

START_TEST (test_get_mean_array) {
    int i, j, ret;
    sample_sum_array * v;
    d_array * x;
    d_array * means;
    int size = 3;
    double e = 0.0000001;
    v = init_sample_sum_array(size);
    x = init_d_array(size);
    means = init_d_array(size);
    for (i = 0; i < size; i++) {
        append_d_array(x, ((i+1)/((double)100)));
    }
    update_sample_sum_array(v, x);
    get_mean_array(v, means);
    ret = d_arrays_equal(means, x, e);
    ck_assert_msg((ret != 0), "did not get expected means");

    update_sample_sum_array(v, x);
    get_mean_array(v, means);
    ret = d_arrays_equal(means, x, e);
    ck_assert_msg((ret != 0), "did not get expected means");

    free_sample_sum_array(v);
    free_d_array(x);
    free_d_array(means);

    v = init_sample_sum_array(size);
    x = init_d_array(size);
    means = init_d_array(size);

    for (i = 0; i < 3; i++) {
        x->length = 0;
        for (j = 0; j < size; j++) {
            append_d_array(x, ((double)(i+1)));
        }
        update_sample_sum_array(v, x);
    }
    get_mean_array(v, means);
    for (i = 0; i < size; i++) {
        ret = almost_equal(get_d_array(means, i), 2.0, e);
        ck_assert_msg((ret != 0), "mean %d is %lf, expecting %lf",
                i, get_d_array(means, i), 2.0);
    }

    free_sample_sum_array(v);
    free_d_array(x);
    free_d_array(means);
}
END_TEST

START_TEST (test_get_sample_variance_array) {
    int i, j, ret;
    sample_sum_array * v;
    d_array * x;
    d_array * vars;
    int size = 3;
    double e = 0.0000001;

    v = init_sample_sum_array(size);
    x = init_d_array(size);
    vars = init_d_array(size);
    
    for (i = 0; i < 3; i++) {
        x->length = 0;
        for (j = 0; j < size; j++) {
            append_d_array(x, ((double)(i+1) * 2));
        }
        update_sample_sum_array(v, x);
    }
    get_sample_variance_array(v, vars);
    for (i = 0; i < size; i++) {
        ret = almost_equal(get_d_array(vars, i), 4.0, e);
        ck_assert_msg((ret != 0), "variance %d is %lf, expecting %lf",
                i, get_d_array(vars, i), 4.0);
    }

    free_sample_sum_array(v);
    free_d_array(x);
    free_d_array(vars);
}
END_TEST

START_TEST (test_get_std_dev_array) {
    int i, j, ret;
    sample_sum_array * v;
    d_array * x;
    d_array * std_devs;
    int size = 3;
    double e = 0.0000001;

    v = init_sample_sum_array(size);
    x = init_d_array(size);
    std_devs = init_d_array(size);
    
    for (i = 0; i < 3; i++) {
        x->length = 0;
        for (j = 0; j < size; j++) {
            append_d_array(x, ((double)(i+1) * 2));
        }
        update_sample_sum_array(v, x);
    }
    get_std_dev_array(v, std_devs);
    for (i = 0; i < size; i++) {
        ret = almost_equal(get_d_array(std_devs, i), 2.0, e);
        ck_assert_msg((ret != 0), "std deviation %d is %lf, expecting %lf",
                i, get_d_array(std_devs, i), 2.0);
    }

    free_sample_sum_array(v);
    free_d_array(x);
    free_d_array(std_devs);
}
END_TEST

START_TEST (test_standardize_vector) {
    int i, size, ret;
    double e = 0.0000001;
    double x, mn, sd, ex;
    size = 5;
    d_array * expected;
    d_array * v;
    d_array * means;
    d_array * std_devs;
    expected = init_d_array(size);
    v = init_d_array(size);
    means = init_d_array(size);
    std_devs = init_d_array(size);
    for (i = 0; i < size; i++) {
        x = (double)(i+1);
        mn = (double)(size-i);
        sd = 2.0;
        append_d_array(expected, ((x-mn)/sd));
        append_d_array(v, x);
        append_d_array(means, mn);
        append_d_array(std_devs, sd);
    }
    standardize_vector(v, means, std_devs);
    ret = d_arrays_equal(v, expected, e);
    ck_assert_msg((ret != 0), "standardized vector is incorrect");

    free_d_array(expected);
    free_d_array(v);
    free_d_array(means);
    free_d_array(std_devs);

    expected = init_d_array(size);
    v = init_d_array(size);
    means = init_d_array(size);
    std_devs = init_d_array(size);
    for (i = 0; i < size; i++) {
        x = (double)(-(i+1));
        mn = (double)(-(size-i));
        sd = 2.0;
        append_d_array(expected, ((x-mn)/sd));
        append_d_array(v, x);
        append_d_array(means, mn);
        append_d_array(std_devs, sd);
    }
    standardize_vector(v, means, std_devs);
    ret = d_arrays_equal(v, expected, e);
    ck_assert_msg((ret != 0), "standardized vector is incorrect");

    free_d_array(expected);
    free_d_array(v);
    free_d_array(means);
    free_d_array(std_devs);
}
END_TEST;


Suite * stats_utils_suite(void) {
    Suite * s = suite_create("stats_utils");

    TCase * tc_sample_sum = tcase_create("sample_sum_test_case");
    tcase_add_test(tc_sample_sum, test_init_free_sample_sum);
    tcase_add_test(tc_sample_sum, test_update_sample_sum);
    tcase_add_test_raise_signal(tc_sample_sum, test_get_mean_n0, SIGABRT);
    tcase_add_test(tc_sample_sum, test_get_mean);
    tcase_add_test_raise_signal(tc_sample_sum, test_get_sample_variance_n0,
            SIGABRT);
    tcase_add_test_raise_signal(tc_sample_sum, test_get_sample_variance_n1,
            SIGABRT);
    tcase_add_test(tc_sample_sum, test_get_sample_variance);
    tcase_add_test_raise_signal(tc_sample_sum, test_get_std_dev_n0,
            SIGABRT);
    tcase_add_test_raise_signal(tc_sample_sum, test_get_std_dev_n1,
            SIGABRT);
    tcase_add_test(tc_sample_sum, test_get_std_dev);
    suite_add_tcase(s, tc_sample_sum);

    TCase * tc_sample_sum_array = tcase_create("sample_sum_array_test_case");
    tcase_add_test_raise_signal(tc_sample_sum_array,
            test_init_sample_sum_array_fail,
            SIGABRT);
    tcase_add_test(tc_sample_sum_array, test_init_free_sample_sum_array);
    tcase_add_test(tc_sample_sum_array, test_update_sample_sum_array);
    tcase_add_test(tc_sample_sum_array, test_get_mean_array);
    tcase_add_test(tc_sample_sum_array, test_get_sample_variance_array);
    tcase_add_test(tc_sample_sum_array, test_get_std_dev_array);
    suite_add_tcase(s, tc_sample_sum_array);

    TCase * tc_standardize_vector = tcase_create(
            "standardize_vector_test_case");
    tcase_add_test(tc_standardize_vector, test_standardize_vector);
    suite_add_tcase(s, tc_standardize_vector);

    return s;
}

int main(void) {
    int number_failed;
    Suite * s = stats_utils_suite();
    SRunner * sr = srunner_create(s);
    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

