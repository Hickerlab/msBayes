#include <stdlib.h>
#include <check.h>
#include <signal.h>
#include "../src/partition_combinatorics.c"
#include "test_utils.h"

START_TEST (test_cumulative_number_of_int_partitions_by_k_n0) {
    int n, ip;
    i_array * counts;
    counts = init_i_array(1);
    n = 0;
    ip = cumulative_number_of_int_partitions_by_k(n, counts); // SIGABRT
    free_i_array(counts);
}
END_TEST

START_TEST (test_cumulative_number_of_int_partitions_by_k_n1) {
    int n, ip;
    i_array * counts;
    counts = init_i_array(1);
    n = 1;
    ip = cumulative_number_of_int_partitions_by_k(n, counts);
    ck_assert_int_eq(ip, 1);
    ck_assert_int_eq(counts->length, n);
    ck_assert_int_eq(get_i_array(counts, 0), 1);
    free_i_array(counts);
}
END_TEST
    
START_TEST (test_cumulative_number_of_int_partitions_by_k_n2) {
    int n, ip;
    i_array * counts;
    counts = init_i_array(1);
    n = 2;
    ip = cumulative_number_of_int_partitions_by_k(n, counts);
    ck_assert_int_eq(ip, 2);
    ck_assert_int_eq(counts->length, n);
    ck_assert_int_eq(get_i_array(counts, 0), 1);
    ck_assert_int_eq(get_i_array(counts, 1), 2);
    free_i_array(counts);
}
END_TEST

START_TEST (test_cumulative_number_of_int_partitions_by_k_n3) {
    int n, ip;
    i_array * counts;
    counts = init_i_array(1);
    n = 3;
    ip = cumulative_number_of_int_partitions_by_k(n, counts);
    ck_assert_int_eq(ip, 3);
    ck_assert_int_eq(counts->length, n);
    ck_assert_int_eq(get_i_array(counts, 0), 1);
    ck_assert_int_eq(get_i_array(counts, 1), 2);
    ck_assert_int_eq(get_i_array(counts, 2), 3);
    free_i_array(counts);
}
END_TEST

START_TEST (test_cumulative_number_of_int_partitions_by_k_n4) {
    int n, ip;
    i_array * counts;
    counts = init_i_array(1);
    n = 4;
    ip = cumulative_number_of_int_partitions_by_k(n, counts);
    ck_assert_int_eq(ip, 5);
    ck_assert_int_eq(counts->length, n);
    ck_assert_int_eq(get_i_array(counts, 0), 1);
    ck_assert_int_eq(get_i_array(counts, 1), 3);
    ck_assert_int_eq(get_i_array(counts, 2), 4);
    ck_assert_int_eq(get_i_array(counts, 3), 5);
    free_i_array(counts);
}
END_TEST

START_TEST (test_cumulative_number_of_int_partitions_by_k_n7) {
    int n, ip;
    i_array * counts;
    counts = init_i_array(1);
    n = 7;
    ip = cumulative_number_of_int_partitions_by_k(n, counts);
    ck_assert_int_eq(ip, 15);
    ck_assert_int_eq(counts->length, n);
    ck_assert_int_eq(get_i_array(counts, 0), 1);
    ck_assert_int_eq(get_i_array(counts, 1), 4);
    ck_assert_int_eq(get_i_array(counts, 2), 8);
    ck_assert_int_eq(get_i_array(counts, 3), 11);
    ck_assert_int_eq(get_i_array(counts, 4), 13);
    ck_assert_int_eq(get_i_array(counts, 5), 14);
    ck_assert_int_eq(get_i_array(counts, 6), 15);
    free_i_array(counts);
}
END_TEST

START_TEST (test_number_of_int_partitions_by_k_n0) {
    int n, ip;
    i_array * counts;
    counts = init_i_array(1);
    n = 0;
    ip = number_of_int_partitions_by_k(n, counts); // SIGABRT
    free_i_array(counts);
}
END_TEST

START_TEST (test_number_of_int_partitions_by_k_n1) {
    int n, ip;
    i_array * counts;
    counts = init_i_array(1);
    n = 1;
    ip = number_of_int_partitions_by_k(n, counts);
    ck_assert_int_eq(ip, 1);
    ck_assert_int_eq(counts->length, n);
    ck_assert_int_eq(get_i_array(counts, 0), 1);
    free_i_array(counts);
}
END_TEST
    
START_TEST (test_number_of_int_partitions_by_k_n2) {
    int n, ip;
    i_array * counts;
    counts = init_i_array(1);
    n = 2;
    ip = number_of_int_partitions_by_k(n, counts);
    ck_assert_int_eq(ip, 2);
    ck_assert_int_eq(counts->length, n);
    ck_assert_int_eq(get_i_array(counts, 0), 1);
    ck_assert_int_eq(get_i_array(counts, 1), 1);
    free_i_array(counts);
}
END_TEST

START_TEST (test_number_of_int_partitions_by_k_n3) {
    int n, ip;
    i_array * counts;
    counts = init_i_array(1);
    n = 3;
    ip = number_of_int_partitions_by_k(n, counts);
    ck_assert_int_eq(ip, 3);
    ck_assert_int_eq(counts->length, n);
    ck_assert_int_eq(get_i_array(counts, 0), 1);
    ck_assert_int_eq(get_i_array(counts, 1), 1);
    ck_assert_int_eq(get_i_array(counts, 2), 1);
    free_i_array(counts);
}
END_TEST

START_TEST (test_number_of_int_partitions_by_k_n4) {
    int n, ip;
    i_array * counts;
    counts = init_i_array(1);
    n = 4;
    ip = number_of_int_partitions_by_k(n, counts);
    ck_assert_int_eq(ip, 5);
    ck_assert_int_eq(counts->length, n);
    ck_assert_int_eq(get_i_array(counts, 0), 1);
    ck_assert_int_eq(get_i_array(counts, 1), 2);
    ck_assert_int_eq(get_i_array(counts, 2), 1);
    ck_assert_int_eq(get_i_array(counts, 3), 1);
    free_i_array(counts);
}
END_TEST

START_TEST (test_number_of_int_partitions_by_k_n7) {
    int n, ip;
    i_array * counts;
    counts = init_i_array(1);
    n = 7;
    ip = number_of_int_partitions_by_k(n, counts);
    ck_assert_int_eq(ip, 15);
    ck_assert_int_eq(counts->length, n);
    ck_assert_int_eq(get_i_array(counts, 0), 1);
    ck_assert_int_eq(get_i_array(counts, 1), 3);
    ck_assert_int_eq(get_i_array(counts, 2), 4);
    ck_assert_int_eq(get_i_array(counts, 3), 3);
    ck_assert_int_eq(get_i_array(counts, 4), 2);
    ck_assert_int_eq(get_i_array(counts, 5), 1);
    ck_assert_int_eq(get_i_array(counts, 6), 1);
    free_i_array(counts);
}
END_TEST

START_TEST (test_cumulative_frequency_of_int_partitions_by_k_n0) {
    int n, tp;
    d_array * probs;
    probs = init_d_array(1);
    n = 0;
    tp = cumulative_frequency_of_int_partitions_by_k(n, probs); // SIGABRT
    free_d_array(probs);
}
END_TEST

START_TEST (test_cumulative_frequency_of_int_partitions_by_k_n1) {
    double e, ep;
    e = 0.000001;
    int n, tp, i;
    d_array * probs;
    probs = init_d_array(1);
    n = 1;
    tp = cumulative_frequency_of_int_partitions_by_k(n, probs);
    ep = 1.0;
    ck_assert_msg((almost_equal(tp, ep, e) != 0),
            "total probability is %lf, not %lf", tp, ep);
    ck_assert_int_eq(probs->length, n);
    i = 0; ep = 1.0;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    free_d_array(probs);
}
END_TEST
    
START_TEST (test_cumulative_frequency_of_int_partitions_by_k_n2) {
    double e, ep;
    e = 0.000001;
    int n, tp, i;
    d_array * probs;
    probs = init_d_array(1);
    n = 2;
    tp = cumulative_frequency_of_int_partitions_by_k(n, probs);
    ep = 1.0;
    ck_assert_msg((almost_equal(tp, ep, e) != 0),
            "total probability is %lf, not %lf", tp, ep);
    ck_assert_int_eq(probs->length, n);
    i = 0; ep = 0.5;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 1; ep = 1.0;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    free_d_array(probs);
}
END_TEST

START_TEST (test_cumulative_frequency_of_int_partitions_by_k_n3) {
    double e, ep;
    e = 0.000001;
    int n, tp, i;
    d_array * probs;
    probs = init_d_array(1);
    n = 3;
    tp = cumulative_frequency_of_int_partitions_by_k(n, probs);
    ep = 1.0;
    ck_assert_msg((almost_equal(tp, ep, e) != 0),
            "total probability is %lf, not %lf", tp, ep);
    ck_assert_int_eq(probs->length, n);
    i = 0; ep = 0.333333333;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 1; ep = 0.666666666;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 2; ep = 1.0;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    free_d_array(probs);
}
END_TEST

START_TEST (test_cumulative_frequency_of_int_partitions_by_k_n4) {
    double e, ep;
    e = 0.000001;
    int n, tp, i, ip;
    d_array * probs;
    probs = init_d_array(1);
    n = 4;
    tp = cumulative_frequency_of_int_partitions_by_k(n, probs);
    ep = 1.0;
    ip = 5;
    ck_assert_msg((almost_equal(tp, ep, e) != 0),
            "total probability is %lf, not %lf", tp, ep);
    ck_assert_int_eq(probs->length, n);
    i = 0; ep = 1/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 1; ep = 3/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 2; ep = 4/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 3; ep = 5/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    free_d_array(probs);
}
END_TEST

START_TEST (test_cumulative_frequency_of_int_partitions_by_k_n7) {
    double e, ep;
    e = 0.000001;
    int n, tp, i, ip;
    d_array * probs;
    probs = init_d_array(1);
    n = 7;
    tp = cumulative_frequency_of_int_partitions_by_k(n, probs);
    ep = 1.0;
    ip = 15;
    ck_assert_msg((almost_equal(tp, ep, e) != 0),
            "total probability is %lf, not %lf", tp, ep);
    ck_assert_int_eq(probs->length, n);
    i = 0; ep = 1/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 1; ep = 4/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 2; ep = 8/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 3; ep = 11/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 4; ep = 13/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 5; ep = 14/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 6; ep = 15/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    free_d_array(probs);
}
END_TEST

START_TEST (test_frequency_of_int_partitions_by_k_n0) {
    int n, tp;
    d_array * probs;
    probs = init_d_array(1);
    n = 0;
    tp = frequency_of_int_partitions_by_k(n, probs); // SIGABRT
    free_d_array(probs);
}
END_TEST

START_TEST (test_frequency_of_int_partitions_by_k_n1) {
    double e, ep;
    e = 0.000001;
    int n, tp, i;
    d_array * probs;
    probs = init_d_array(1);
    n = 1;
    tp = frequency_of_int_partitions_by_k(n, probs);
    ep = 1.0;
    ck_assert_msg((almost_equal(tp, ep, e) != 0),
            "total probability is %lf, not %lf", tp, ep);
    ck_assert_int_eq(probs->length, n);
    i = 0; ep = 1.0;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    free_d_array(probs);
}
END_TEST
    
START_TEST (test_frequency_of_int_partitions_by_k_n2) {
    double e, ep;
    e = 0.000001;
    int n, tp, i;
    d_array * probs;
    probs = init_d_array(1);
    n = 2;
    tp = frequency_of_int_partitions_by_k(n, probs);
    ep = 1.0;
    ck_assert_msg((almost_equal(tp, ep, e) != 0),
            "total probability is %lf, not %lf", tp, ep);
    ck_assert_int_eq(probs->length, n);
    i = 0; ep = 0.5;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 1; ep = 0.5;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    free_d_array(probs);
}
END_TEST

START_TEST (test_frequency_of_int_partitions_by_k_n3) {
    double e, ep;
    e = 0.000001;
    int n, tp, i;
    d_array * probs;
    probs = init_d_array(1);
    n = 3;
    tp = frequency_of_int_partitions_by_k(n, probs);
    ep = 1.0;
    ck_assert_msg((almost_equal(tp, ep, e) != 0),
            "total probability is %lf, not %lf", tp, ep);
    ck_assert_int_eq(probs->length, n);
    i = 0; ep = 0.333333333;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 1; ep = 0.333333333;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 2; ep = 0.333333333;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    free_d_array(probs);
}
END_TEST

START_TEST (test_frequency_of_int_partitions_by_k_n4) {
    double e, ep;
    e = 0.000001;
    int n, tp, i, ip;
    d_array * probs;
    probs = init_d_array(1);
    n = 4;
    tp = frequency_of_int_partitions_by_k(n, probs);
    ep = 1.0;
    ip = 5;
    ck_assert_msg((almost_equal(tp, ep, e) != 0),
            "total probability is %lf, not %lf", tp, ep);
    ck_assert_int_eq(probs->length, n);
    i = 0; ep = 1/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 1; ep = 2/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 2; ep = 1/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 3; ep = 1/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    free_d_array(probs);
}
END_TEST

START_TEST (test_frequency_of_int_partitions_by_k_n7) {
    double e, ep;
    e = 0.000001;
    int n, tp, i, ip;
    d_array * probs;
    probs = init_d_array(1);
    n = 7;
    tp = frequency_of_int_partitions_by_k(n, probs);
    ep = 1.0;
    ip = 15;
    ck_assert_msg((almost_equal(tp, ep, e) != 0),
            "total probability is %lf, not %lf", tp, ep);
    ck_assert_int_eq(probs->length, n);
    i = 0; ep = 1/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 1; ep = 3/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 2; ep = 4/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 3; ep = 3/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 4; ep = 2/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 5; ep = 1/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    i = 6; ep = 1/(double)ip;
    ck_assert_msg((almost_equal(get_d_array(probs, i), ep, e) != 0),
            "cat %d has probability is %lf, not %lf",
            i+1, get_d_array(probs, i), ep);
    free_d_array(probs);
}
END_TEST


START_TEST (test_number_of_int_partitions_neg) {
    int n, ret;
    n = -1;
    ret = number_of_int_partitions(n);
    ck_assert_int_eq(ret, 0);
}
END_TEST

START_TEST (test_number_of_int_partitions_n0) {
    int n, ret;
    n = 0;
    ret = number_of_int_partitions(n);
    ck_assert_int_eq(ret, 1);
}
END_TEST

START_TEST (test_number_of_int_partitions_n1) {
    int n, ret;
    n = 1;
    ret = number_of_int_partitions(n);
    ck_assert_int_eq(ret, 1);
}
END_TEST

START_TEST (test_number_of_int_partitions_n7) {
    int n, ret;
    n = 7;
    ret = number_of_int_partitions(n);
    ck_assert_int_eq(ret, 15);
}
END_TEST

START_TEST (test_number_of_int_partitions_n22) {
    int n, ret;
    n = 22;
    ret = number_of_int_partitions(n);
    ck_assert_int_eq(ret, 1002);
}
END_TEST


START_TEST (test_generate_int_partitions_n0) {
    int n;
    n = 0;
    i_array_2d * partitions = init_i_array_2d(1,1);
    partitions = generate_int_partitions(n);
    free_i_array_2d(partitions);
}
END_TEST

START_TEST (test_generate_int_partitions_n1) {
    int n;
    n = 1;
    i_array * expected = init_i_array(1);
    i_array_2d * partitions = init_i_array_2d(1,1);
    partitions = generate_int_partitions(n);
    ck_assert_int_eq(partitions->length, 1);
    append_i_array(expected, 1);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, 0),
            expected) != 0),
            "unexpected partition");
    free_i_array(expected);
    free_i_array_2d(partitions);
}
END_TEST

START_TEST (test_generate_int_partitions_n4) {
    int n, i;
    n = 4;
    i_array * expected = init_i_array(1);
    i_array_2d * partitions = init_i_array_2d(1,1);
    partitions = generate_int_partitions(n);
    ck_assert_int_eq(partitions->length, 5);
    i = 0;
    expected->length = 0;
    append_i_array(expected, 4);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 1;
    expected->length = 0;
    append_i_array(expected, 3);
    append_i_array(expected, 1);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 2;
    expected->length = 0;
    append_i_array(expected, 2);
    append_i_array(expected, 2);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 3;
    expected->length = 0;
    append_i_array(expected, 2);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 4;
    expected->length = 0;
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    free_i_array(expected);
    free_i_array_2d(partitions);
}
END_TEST

START_TEST (test_generate_int_partitions_n6) {
    int n, i;
    n = 6;
    i_array * expected = init_i_array(1);
    i_array_2d * partitions = init_i_array_2d(1,1);
    partitions = generate_int_partitions(n);
    ck_assert_int_eq(partitions->length, 11);
    i = 0;
    expected->length = 0;
    append_i_array(expected, 6);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 1;
    expected->length = 0;
    append_i_array(expected, 5);
    append_i_array(expected, 1);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 2;
    expected->length = 0;
    append_i_array(expected, 4);
    append_i_array(expected, 2);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 3;
    expected->length = 0;
    append_i_array(expected, 4);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 4;
    expected->length = 0;
    append_i_array(expected, 3);
    append_i_array(expected, 3);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 5;
    expected->length = 0;
    append_i_array(expected, 3);
    append_i_array(expected, 2);
    append_i_array(expected, 1);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 6;
    expected->length = 0;
    append_i_array(expected, 3);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 7;
    expected->length = 0;
    append_i_array(expected, 2);
    append_i_array(expected, 2);
    append_i_array(expected, 2);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 8;
    expected->length = 0;
    append_i_array(expected, 2);
    append_i_array(expected, 2);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 9;
    expected->length = 0;
    append_i_array(expected, 2);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    i = 10;
    expected->length = 0;
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    append_i_array(expected, 1);
    ck_assert_msg((i_arrays_equal(get_i_array_2d(partitions, i),
            expected) != 0),
            "unexpected partition");
    free_i_array(expected);
    free_i_array_2d(partitions);
}
END_TEST


Suite * partition_combinatorics_suite(void) {
    Suite * s = suite_create("partition_combinatorics");

    TCase * tc_partition_counts = tcase_create("partition_count");
    tcase_add_test_raise_signal(tc_partition_counts,
            test_cumulative_number_of_int_partitions_by_k_n0,
            SIGABRT);
    tcase_add_test(tc_partition_counts,
            test_cumulative_number_of_int_partitions_by_k_n1);
    tcase_add_test(tc_partition_counts,
            test_cumulative_number_of_int_partitions_by_k_n2);
    tcase_add_test(tc_partition_counts,
            test_cumulative_number_of_int_partitions_by_k_n3);
    tcase_add_test(tc_partition_counts,
            test_cumulative_number_of_int_partitions_by_k_n4);
    tcase_add_test(tc_partition_counts,
            test_cumulative_number_of_int_partitions_by_k_n7);
    tcase_add_test_raise_signal(tc_partition_counts,
            test_number_of_int_partitions_by_k_n0,
            SIGABRT);
    tcase_add_test(tc_partition_counts,
            test_number_of_int_partitions_by_k_n1);
    tcase_add_test(tc_partition_counts,
            test_number_of_int_partitions_by_k_n2);
    tcase_add_test(tc_partition_counts,
            test_number_of_int_partitions_by_k_n3);
    tcase_add_test(tc_partition_counts,
            test_number_of_int_partitions_by_k_n4);
    tcase_add_test(tc_partition_counts,
            test_number_of_int_partitions_by_k_n7);
    suite_add_tcase(s, tc_partition_counts);

    TCase * tc_partition_freqs = tcase_create("partition_freq");
    tcase_add_test_raise_signal(tc_partition_freqs,
            test_cumulative_frequency_of_int_partitions_by_k_n0,
            SIGABRT);
    tcase_add_test(tc_partition_freqs,
            test_cumulative_frequency_of_int_partitions_by_k_n1);
    tcase_add_test(tc_partition_freqs,
            test_cumulative_frequency_of_int_partitions_by_k_n2);
    tcase_add_test(tc_partition_freqs,
            test_cumulative_frequency_of_int_partitions_by_k_n3);
    tcase_add_test(tc_partition_freqs,
            test_cumulative_frequency_of_int_partitions_by_k_n4);
    tcase_add_test(tc_partition_freqs,
            test_cumulative_frequency_of_int_partitions_by_k_n7);
    tcase_add_test_raise_signal(tc_partition_freqs,
            test_frequency_of_int_partitions_by_k_n0,
            SIGABRT);
    tcase_add_test(tc_partition_freqs,
            test_frequency_of_int_partitions_by_k_n1);
    tcase_add_test(tc_partition_freqs,
            test_frequency_of_int_partitions_by_k_n2);
    tcase_add_test(tc_partition_freqs,
            test_frequency_of_int_partitions_by_k_n3);
    tcase_add_test(tc_partition_freqs,
            test_frequency_of_int_partitions_by_k_n4);
    tcase_add_test(tc_partition_freqs,
            test_frequency_of_int_partitions_by_k_n7);
    suite_add_tcase(s, tc_partition_freqs);

    TCase * tc_num_partitions = tcase_create("number_of_partitions");
    tcase_add_test(tc_num_partitions,
            test_number_of_int_partitions_neg);
    tcase_add_test(tc_num_partitions,
            test_number_of_int_partitions_n0);
    tcase_add_test(tc_num_partitions,
            test_number_of_int_partitions_n1);
    tcase_add_test(tc_num_partitions,
            test_number_of_int_partitions_n7);
    tcase_add_test(tc_num_partitions,
            test_number_of_int_partitions_n22);
    suite_add_tcase(s, tc_num_partitions);

    TCase * tc_generate_int_partitions = tcase_create(
            "generate_int_partitions");
    tcase_add_test_raise_signal(tc_generate_int_partitions,
            test_generate_int_partitions_n0,
            SIGABRT);
    tcase_add_test(tc_generate_int_partitions,
            test_generate_int_partitions_n1);
    tcase_add_test(tc_generate_int_partitions,
            test_generate_int_partitions_n4);
    tcase_add_test(tc_generate_int_partitions,
            test_generate_int_partitions_n6);
    suite_add_tcase(s, tc_generate_int_partitions);

    return s;
}

int main(void) {
    int number_failed;
    Suite * s = partition_combinatorics_suite();
    SRunner * sr = srunner_create(s);
    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

