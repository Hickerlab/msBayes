#include <stdlib.h>
#include <check.h>
#include <signal.h>
#include "../src/math_utils.c"
#include "test_utils.h"


START_TEST(test_get_euclidean_distance) {
    int i;
    double e = 0.0000001;
    int size = 5;
    d_array * v1 = init_d_array(size);
    d_array * v2 = init_d_array(size);
    for (i = 0; i < size; i++) {
        append_d_array(v1, (i+1));
        append_d_array(v2, (size-i));
    }
    double d = get_euclidean_distance(v1, v2);
    int ret = almost_equal(d, sqrt(40.0), e);
    ck_assert_msg((ret != 0), "euclidean distance is %lf, expecting "
            "%lf", d, sqrt(40.0));
    free_d_array(v1);
    free_d_array(v2);

    v1 = init_d_array(size);
    v2 = init_d_array(size);
    for (i = 0; i < size; i++) {
        append_d_array(v1, -(i+1));
        append_d_array(v2, -(size-i));
    }
    d = get_euclidean_distance(v1, v2);
    ret = almost_equal(d, sqrt(40.0), e);
    ck_assert_msg((ret != 0), "euclidean distance is %lf, expecting "
            "%lf", d, sqrt(40.0));
    free_d_array(v1);
    free_d_array(v2);
}
END_TEST

Suite * math_utils_suite(void) {
    Suite * s = suite_create("math_utils");

    TCase * tc_get_euclidean_distance = tcase_create(
            "euclidean_distance_test_case");
    tcase_add_test(tc_get_euclidean_distance, test_get_euclidean_distance);
    suite_add_tcase(s, tc_get_euclidean_distance);

    return s;
}

int main(void) {
    int number_failed;
    Suite * s = math_utils_suite();
    SRunner * sr = srunner_create(s);
    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

