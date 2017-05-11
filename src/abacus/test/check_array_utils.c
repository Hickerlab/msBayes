#include <stdlib.h>
#include <check.h>
#include <signal.h>
#include "test_utils.h"
#include "../src/array_utils.c"

/**
 * `d_array` tests
 */
START_TEST (test_init_free_d_array) {
    int i, size;
    d_array * v;
    size = 10;
    v = init_d_array(size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 0);
    ck_assert_msg((sizeof(v->length) == sizeof(int)), "`d_array.length` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->capacity) == sizeof(int)), "`d_array.capacity` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->a) == sizeof(long)), "Array attribute "
            "of `d_array` is not a pointer");
    ck_assert_msg((sizeof(*v->a) == sizeof(double)), "Pointer to array"
            "of `d_array` is not pointing to double");
    for (i = 0; i < v->capacity; i++) {
        ck_assert_msg((v->a[i] == 0), "Element %d of %d was not initialized "
            "to zero", i, size);
    }
    free_d_array(v);
}
END_TEST

START_TEST (test_init_d_array_fail) {
    d_array * v;
    v = init_d_array(0); // SIGABRT
}
END_TEST

START_TEST (test_expand_d_array) {
    int i, size;
    d_array * v;
    size = 1;
    v = init_d_array(size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 0);
    for (i = 0; i < 5; i++) {
        expand_d_array(v);
        size *= 2;
        ck_assert_int_eq(v->capacity, size);
        ck_assert_int_eq(v->length, 0);
    }
    free_d_array(v);
}
END_TEST

START_TEST (test_append_d_array) {
    int size;
    d_array * v;
    size = 1;
    v = init_d_array(size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 0);
    append_d_array(v, 0.1);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 1);
    ck_assert_msg((v->a[0] == 0.1), "`d_array` element is %lf, expecting %lf",
            v->a[0], 0.1);
    append_d_array(v, 1.9);
    ck_assert_int_eq(v->capacity, size*2);
    ck_assert_int_eq(v->length, 2);
    ck_assert_msg((v->a[1] == 1.9), "`d_array` element is %lf, expecting %lf",
            v->a[0], 1.9);
}
END_TEST

START_TEST (test_set_d_array) {
    int size;
    d_array * v;
    size = 1;
    v = init_d_array(size);
    append_d_array(v, 0.1);
    ck_assert_msg((v->a[0] == 0.1), "`d_array` element is %lf, expecting %lf",
            v->a[0], 0.1);
    set_d_array(v, 0, 0.2);
    ck_assert_msg((v->a[0] == 0.2), "`d_array` element is %lf, expecting %lf",
            v->a[0], 0.2);
}
END_TEST

START_TEST (test_set_d_array_fail) {
    int size;
    d_array * v;
    size = 1;
    v = init_d_array(size);
    append_d_array(v, 0.1);
    set_d_array(v, 1, 0.3); // SIGABRT
}
END_TEST

START_TEST (test_get_d_array) {
    int size;
    d_array * v;
    size = 1;
    v = init_d_array(size);
    append_d_array(v, 0.1);
    ck_assert_msg((get_d_array(v, 0) == 0.1), "`d_array` element is %lf, expecting %lf",
            get_d_array(v, 0), 0.1);
    append_d_array(v, 0.2);
    ck_assert_msg((get_d_array(v, 1) == 0.2), "`d_array` element is %lf, expecting %lf",
            get_d_array(v, 1), 0.2);
}
END_TEST

START_TEST (test_get_d_array_fail) {
    int size;
    d_array * v;
    size = 1;
    v = init_d_array(size);
    append_d_array(v, 0.1);
    get_d_array(v, 1); // SIGABRT
}
END_TEST

/**
 * `i_array` tests
 */
START_TEST (test_init_free_i_array) {
    int i, size;
    i_array * v;
    size = 10;
    v = init_i_array(size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 0);
    ck_assert_msg((sizeof(v->length) == sizeof(int)), "`i_array.length` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->capacity) == sizeof(int)), "`i_array.capacity` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->a) == sizeof(long)), "Array attribute "
            "of `i_array` is not a pointer");
    ck_assert_msg((sizeof(*v->a) == sizeof(int)), "Pointer to array"
            "of `i_array` is not pointing to int");
    for (i = 0; i < v->capacity; i++) {
        ck_assert_msg((v->a[i] == 0), "Element %d of %d was not initialized "
            "to zero", i, size);
    }
    free_i_array(v);
}
END_TEST

START_TEST (test_init_i_array_fail) {
    i_array * v;
    v = init_i_array(0); // SIGABRT
}
END_TEST

START_TEST (test_expand_i_array) {
    int i, size;
    i_array * v;
    size = 1;
    v = init_i_array(size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 0);
    for (i = 0; i < 5; i++) {
        expand_i_array(v);
        size *= 2;
        ck_assert_int_eq(v->capacity, size);
        ck_assert_int_eq(v->length, 0);
    }
    free_i_array(v);
}
END_TEST

START_TEST (test_append_i_array) {
    int size;
    i_array * v;
    size = 1;
    v = init_i_array(size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 0);
    append_i_array(v, 1);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 1);
    ck_assert_msg((v->a[0] == 1), "`i_array` element is %d, expecting %d",
            v->a[0], 1);
    append_i_array(v, -99);
    ck_assert_int_eq(v->capacity, size*2);
    ck_assert_int_eq(v->length, 2);
    ck_assert_msg((v->a[1] == -99), "`i_array` element is %d, expecting %d",
            v->a[0], -99);
}
END_TEST

START_TEST (test_set_i_array) {
    int size;
    i_array * v;
    size = 1;
    v = init_i_array(size);
    append_i_array(v, 1);
    ck_assert_msg((v->a[0] == 1), "`i_array` element is %d, expecting %d",
            v->a[0], 1);
    set_i_array(v, 0, 2);
    ck_assert_msg((v->a[0] == 2), "`i_array` element is %d, expecting %d",
            v->a[0], 2);
}
END_TEST

START_TEST (test_set_i_array_fail) {
    int size;
    i_array * v;
    size = 1;
    v = init_i_array(size);
    append_i_array(v, 1);
    set_i_array(v, 1, 3); // SIGABRT
}
END_TEST

START_TEST (test_get_i_array) {
    int size;
    i_array * v;
    size = 1;
    v = init_i_array(size);
    append_i_array(v, 1);
    ck_assert_msg((get_i_array(v, 0) == 1), "`i_array` element is %d, expecting %d",
            get_i_array(v, 0), 1);
    append_i_array(v, 2);
    ck_assert_msg((get_i_array(v, 1) == 2), "`i_array` element is %d, expecting %d",
            get_i_array(v, 1), 2);
}
END_TEST

START_TEST (test_get_i_array_fail) {
    int size;
    i_array * v;
    size = 1;
    v = init_i_array(size);
    append_i_array(v, 1);
    get_i_array(v, 1); // SIGABRT
}
END_TEST

/**
 * `c_array` tests
 */
START_TEST (test_init_free_c_array) {
    int i, size;
    c_array * v;
    size = 10;
    v = init_c_array(size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_msg((sizeof(v->capacity) == sizeof(int)), "`c_array.capacity` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->a) == sizeof(long)), "Array attribute "
            "of `c_array` is not a pointer");
    ck_assert_msg((sizeof(*v->a) == sizeof(char)), "Pointer to array"
            "of `c_array` is not pointing to char");
    for (i = 0; i < v->capacity; i++) {
        ck_assert_msg((v->a[i] == 0), "Element %d of %d was not initialized "
            "to zero", i, size);
    }
    ck_assert_msg((v->a[v->capacity] == '\0'), "`c_array` does not have "
            "extra NULL terminating character");
    free_c_array(v);
}
END_TEST

START_TEST (test_init_c_array_fail) {
    c_array * v;
    v = init_c_array(0); // SIGABRT
}
END_TEST

START_TEST (test_expand_c_array) {
    int i, size;
    c_array * v;
    size = 7;
    v = init_c_array(size);
    ck_assert_int_eq(v->capacity, size);
    for (i = 0; i < 3; i++) {
        expand_c_array(v);
        size = ((size + 1) * 2) - 1;
        ck_assert_int_eq(v->capacity, size);
        ck_assert_msg((v->a[v->capacity] == '\0'), "`c_array` does not "
            "have extra NULL terminating character");
    }
    free_c_array(v);
}
END_TEST

START_TEST (test_assign_c_array) {
    int size;
    c_array * v;
    char * s = "blah";
    size = 1;
    v = init_c_array(size);
    assign_c_array(v, s);
    ck_assert_int_eq(v->capacity, 7);
    ck_assert_msg((v->a[8] == '\0'), "`c_array` does not have "
            "extra NULL terminating character");
    ck_assert_msg((v->a[5] == '\0'), "`c_array` did not copy "
            "extra NULL terminating character of target");
    ck_assert_msg((v->a[0] == 'b'), "`c_array` element is %c, expecting %c",
            v->a[0], 'b');
    ck_assert_msg((v->a[1] == 'l'), "`c_array` element is %c, expecting %c",
            v->a[1], 'l');
    ck_assert_msg((v->a[2] == 'a'), "`c_array` element is %c, expecting %c",
            v->a[2], 'a');
    ck_assert_msg((v->a[3] == 'h'), "`c_array` element is %c, expecting %c",
            v->a[3], 'h');
    free_c_array(v);
}
END_TEST

START_TEST (test_get_c_array) {
    int size;
    c_array * v;
    char * s = "blah";
    size = 1;
    v = init_c_array(size);
    assign_c_array(v, s);
    ck_assert_msg((*(get_c_array(v)) == *s), "`c_array` is %s, expecting %s",
            get_c_array(v), s);
    free_c_array(v);
}
END_TEST

/**
 * `s_array` tests
 */
START_TEST (test_init_free_s_array) {
    int i, size;
    s_array * v;
    size = 10;
    v = init_s_array(size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 0);
    ck_assert_msg((sizeof(v->length) == sizeof(int)), "`s_array.length` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->capacity) == sizeof(int)), "`s_array.capacity` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->a) == sizeof(long)), "Array attribute "
            "of `s_array` is not a pointer");
    ck_assert_msg((sizeof(*v->a) == sizeof(c_array *)), "Pointer to array"
            "of `s_array` is not pointing to array of `c_array`s");
    for (i = 0; i < v->capacity; i++) {
        ck_assert_msg((v->a[i]->capacity > 0), "`c_array` element %d of %d was "
                "not initialized ", i, size);
    }
    free_s_array(v);
}
END_TEST

START_TEST (test_init_s_array_fail) {
    s_array * v;
    v = init_s_array(0); // SIGABRT
}
END_TEST

START_TEST (test_expand_s_array) {
    int i, size;
    s_array * v;
    size = 1;
    v = init_s_array(size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 0);
    for (i = 0; i < 5; i++) {
        expand_s_array(v);
        size *= 2;
        ck_assert_int_eq(v->capacity, size);
        ck_assert_int_eq(v->length, 0);
    }
    for (i = 0; i < v->capacity; i++) {
        ck_assert_msg((v->a[i]->capacity > 0), "`c_array` element %d of %d was "
                "not initialized ", i, size);
    }
    free_s_array(v);
}
END_TEST

START_TEST (test_append_s_array) {
    int size;
    s_array * v;
    char * s = "foo";
    char * s2 = "bar";
    size = 1;
    v = init_s_array(size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 0);
    append_s_array(v, s);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 1);
    ck_assert_msg((*v->a[0]->a == *s), "`s_array` element is %s, expecting %s",
            v->a[0]->a, s);
    append_s_array(v, s2);
    ck_assert_int_eq(v->capacity, size*2);
    ck_assert_int_eq(v->length, 2);
    ck_assert_msg((*v->a[1]->a == *s2), "`s_array` element is %s, expecting %s",
            v->a[0]->a, s2);
}
END_TEST

START_TEST (test_set_s_array) {
    int size;
    s_array * v;
    char * s = "foo";
    char * s2 = "bar";
    size = 1;
    v = init_s_array(size);
    append_s_array(v, s);
    ck_assert_msg((*v->a[0]->a == *s), "`s_array` element is %s, expecting %s",
            v->a[0], s);
    set_s_array(v, 0, s2);
    ck_assert_msg((*v->a[0]->a == *s2), "`s_array` element is %d, expecting %d",
            v->a[0], s2);
}
END_TEST

START_TEST (test_set_s_array_fail) {
    int size;
    s_array * v;
    char * s = "foo";
    size = 1;
    v = init_s_array(size);
    append_s_array(v, s);
    set_s_array(v, 1, s); // SIGABRT
}
END_TEST

START_TEST (test_get_s_array) {
    int size;
    s_array * v;
    char * s = "foo";
    char * s2 = "bar";
    size = 1;
    v = init_s_array(size);
    append_s_array(v, s);
    ck_assert_msg((*(get_s_array(v, 0)) == *s), "`s_array` element is %s, "
            "expecting %s", get_s_array(v, 0), s);
    append_s_array(v, s2);
    ck_assert_msg((*(get_s_array(v, 1)) == *s2), "`s_array` element is %s, "
            "expecting %s", get_s_array(v, 1), s2);
}
END_TEST

START_TEST (test_get_s_array_fail) {
    int size;
    s_array * v;
    char * s = "foo";
    size = 1;
    v = init_s_array(size);
    append_s_array(v, s);
    get_s_array(v, 1); // SIGABRT
}
END_TEST

START_TEST (test_d_arrays_equal) {
    int size, ret;
    d_array * v1;
    d_array * v2;
    double s = 1.0;
    double e = 0.0000001;
    size = 3;
    v1 = init_d_array(size);
    v2 = init_d_array(size);
    ret = d_arrays_equal(v1, v2, e);
    ck_assert_msg((ret != 0), "`d_array`s not equal, but should be");
    append_d_array(v1, s);
    ret = d_arrays_equal(v1, v2, e);
    ck_assert_msg((ret == 0), "`d_array`s equal, but should not be");
    append_d_array(v2, s);
    ret = d_arrays_equal(v1, v2, e);
    ck_assert_msg((ret != 0), "`d_array`s not equal, but should be");
    s = 2.0;
    append_d_array(v1, s);
    ret = d_arrays_equal(v1, v2, e);
    ck_assert_msg((ret == 0), "`d_array`s equal, but should not be");
    append_d_array(v2, s);
    ret = d_arrays_equal(v1, v2, e);
    ck_assert_msg((ret != 0), "`d_array`s not equal, but should be");
    append_d_array(v1, s);
    ret = d_arrays_equal(v1, v2, e);
    ck_assert_msg((ret == 0), "`d_array`s equal, but should not be");
    s = 3.0;
    append_d_array(v2, s);
    ret = d_arrays_equal(v1, v2, e);
    ck_assert_msg((ret == 0), "`d_array`s equal, but should not be");
    free_d_array(v1);
    free_d_array(v2);
}
END_TEST

START_TEST (test_i_arrays_equal) {
    int size, ret;
    i_array * v1;
    i_array * v2;
    int s = 1;
    size = 3;
    v1 = init_i_array(size);
    v2 = init_i_array(size);
    ret = i_arrays_equal(v1, v2);
    ck_assert_msg((ret != 0), "`i_array`s not equal, but should be");
    append_i_array(v1, s);
    ret = i_arrays_equal(v1, v2);
    ck_assert_msg((ret == 0), "`i_array`s equal, but should not be");
    append_i_array(v2, s);
    ret = i_arrays_equal(v1, v2);
    ck_assert_msg((ret != 0), "`i_array`s not equal, but should be");
    s = 2;
    append_i_array(v1, s);
    ret = i_arrays_equal(v1, v2);
    ck_assert_msg((ret == 0), "`i_array`s equal, but should not be");
    append_i_array(v2, s);
    ret = i_arrays_equal(v1, v2);
    ck_assert_msg((ret != 0), "`i_array`s not equal, but should be");
    append_i_array(v1, s);
    ret = i_arrays_equal(v1, v2);
    ck_assert_msg((ret == 0), "`i_array`s equal, but should not be");
    s = 3;
    append_i_array(v2, s);
    ret = i_arrays_equal(v1, v2);
    ck_assert_msg((ret == 0), "`i_array`s equal, but should not be");
    free_i_array(v1);
    free_i_array(v2);
}
END_TEST

START_TEST (test_s_arrays_equal) {
    int size, ret;
    s_array * v1;
    s_array * v2;
    char * s = "foo";
    size = 3;
    v1 = init_s_array(size);
    v2 = init_s_array(size);
    ret = s_arrays_equal(v1, v2);
    ck_assert_msg((ret != 0), "`s_array`s not equal, but should be");
    append_s_array(v1, s);
    ret = s_arrays_equal(v1, v2);
    ck_assert_msg((ret == 0), "`s_array`s equal, but should not be");
    append_s_array(v2, s);
    ret = s_arrays_equal(v1, v2);
    ck_assert_msg((ret != 0), "`s_array`s not equal, but should be");
    s = "bar";
    append_s_array(v1, s);
    ret = s_arrays_equal(v1, v2);
    ck_assert_msg((ret == 0), "`s_array`s equal, but should not be");
    append_s_array(v2, s);
    ret = s_arrays_equal(v1, v2);
    ck_assert_msg((ret != 0), "`s_array`s not equal, but should be");
    append_s_array(v1, s);
    ret = s_arrays_equal(v1, v2);
    ck_assert_msg((ret == 0), "`s_array`s equal, but should not be");
    s = "blah";
    append_s_array(v2, s);
    ret = s_arrays_equal(v1, v2);
    ck_assert_msg((ret == 0), "`s_array`s equal, but should not be");
    free_s_array(v1);
    free_s_array(v2);
}
END_TEST

START_TEST (test_get_matching_indices) {
    s_array * search_strings;
    s_array * target_strings;
    i_array * indices;
    search_strings = init_s_array(1);
    target_strings = init_s_array(1);
    indices = init_i_array(1);
    append_s_array(search_strings, "foo");
    append_s_array(search_strings, "bar");
    append_s_array(target_strings, "boo");
    append_s_array(target_strings, "foo");
    append_s_array(target_strings, "bar");
    get_matching_indices(search_strings, target_strings, indices);
    ck_assert_int_eq(indices->length, 2);
    ck_assert_int_eq(indices->capacity, 2);
    ck_assert_int_eq(get_i_array(indices, 0), 1);
    ck_assert_int_eq(get_i_array(indices, 1), 2);
}
END_TEST

START_TEST (test_get_matching_indices_fail) {
    s_array * search_strings;
    s_array * target_strings;
    i_array * indices;
    search_strings = init_s_array(1);
    target_strings = init_s_array(1);
    indices = init_i_array(1);
    append_s_array(search_strings, "foo");
    append_s_array(search_strings, "bar");
    append_s_array(target_strings, "boo");
    append_s_array(target_strings, "foo");
    append_s_array(target_strings, "bar");
    append_s_array(search_strings, "fail");
    get_matching_indices(search_strings, target_strings, indices); // exit(1)
}
END_TEST

START_TEST (test_get_doubles) {
    int ret;
    s_array * search_strings;
    i_array * indices;
    d_array * dest;
    search_strings = init_s_array(1);
    indices = init_i_array(1);
    dest = init_d_array(1);
    append_s_array(search_strings, "foo");
    append_s_array(search_strings, "0.345");
    append_s_array(search_strings, "156.345");
    append_s_array(search_strings, "bar");
    append_s_array(search_strings, "1.23e3");
    append_s_array(search_strings, "1.2");
    append_i_array(indices, 1);
    append_i_array(indices, 2);
    append_i_array(indices, 4);
    ret = get_doubles(search_strings, indices, dest);
    ck_assert_int_eq(ret, 0);
    ck_assert_int_eq(dest->length, 3);
    ck_assert_msg((almost_equal(get_d_array(dest, 0), 0.345, 0.000001)),
            "extracted double is %lf, expected %lf",
            get_d_array(dest, 0), 0.345);
    ck_assert_msg((almost_equal(get_d_array(dest, 1), 156.345, 0.000001)),
            "extracted double is %lf, expected %lf",
            get_d_array(dest, 1), 156.345);
    ck_assert_msg((almost_equal(get_d_array(dest, 2), 1230.0, 0.000001)),
            "extracted double is %lf, expected %lf",
            get_d_array(dest, 2), 1230.0);
    append_i_array(indices, 0);
    ret = get_doubles(search_strings, indices, dest);
    ck_assert_int_eq(ret, 1);
}
END_TEST

/**
 * `i_array_2d` tests
 */
START_TEST (test_init_free_i_array_2d) {
    int i, size, el_size;
    i_array_2d * v;
    size = 10;
    el_size = 5;
    v = init_i_array_2d(size, el_size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->initial_element_capacity, el_size);
    ck_assert_int_eq(v->length, 0);
    ck_assert_msg((sizeof(v->length) == sizeof(int)), "`i_array_2d.length` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->capacity) == sizeof(int)), "`i_array_2d.capacity` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->a) == sizeof(long)), "Array attribute "
            "of `i_array_2d` is not a pointer");
    ck_assert_msg((sizeof(*v->a) == sizeof(i_array *)), "Pointer to array"
            "of `i_array_2d` is not pointing to array of `i_array`s");
    for (i = 0; i < v->capacity; i++) {
        ck_assert_msg((v->a[i]->capacity == el_size),
                "`i_array` element %d of %d was initialized with capacity %d, "
                "expecting %d", i, size, v->a[i]->capacity, el_size);
    }
    free_i_array_2d(v);
}
END_TEST

START_TEST (test_init_i_array_2d_fail) {
    i_array_2d * v;
    v = init_i_array_2d(0, 1); // SIGABRT
}
END_TEST

START_TEST (test_expand_i_array_2d) {
    int i, size, el_size;
    i_array_2d * v;
    size = 1;
    el_size = 1;
    v = init_i_array_2d(size, el_size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 0);
    for (i = 0; i < 5; i++) {
        expand_i_array_2d(v);
        size *= 2;
        ck_assert_int_eq(v->capacity, size);
        ck_assert_int_eq(v->length, 0);
    }
    for (i = 0; i < v->capacity; i++) {
        ck_assert_msg((v->a[i]->capacity == el_size), "`i_array` element %d of "
                "%d was initialized with capacity %d, expecting %d",
                i, size, v->a[i]->capacity, el_size);
    }
    free_i_array_2d(v);
}
END_TEST

START_TEST (test_append_i_array_2d) {
    int size, el_size, i, ret;
    i_array_2d * v;
    i_array * ia1;
    i_array * ia2;
    el_size = 3;
    ia1 = init_i_array(el_size);
    ia2 = init_i_array(el_size);
    for (i = 0; i < el_size; i++) {
        append_i_array(ia1, i);
        append_i_array(ia2, (el_size + i));
    }
    size = 1;
    v = init_i_array_2d(size, el_size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 0);
    append_i_array_2d(v, ia1);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 1);
    ret = i_arrays_equal(v->a[0], ia1);
    ck_assert_msg((ret != 0), "`i_array_2d` element is incorrect");
    append_i_array_2d(v, ia2);
    ck_assert_int_eq(v->capacity, size*2);
    ck_assert_int_eq(v->length, 2);
    ret = i_arrays_equal(v->a[1], ia2);
    ck_assert_msg((ret != 0), "`i_array_2d` element is incorrect");

    free_i_array_2d(v);
    free_i_array(ia1);
    free_i_array(ia2);
}
END_TEST

START_TEST (test_append_el_i_array_2d) {
    int size, el_size, i, ret;
    i_array_2d * v;
    i_array * ia1;
    i_array * ia2;
    el_size = 3;
    size = 1;
    v = init_i_array_2d(1, 1);
    ia1 = init_i_array(el_size);
    ia2 = init_i_array(el_size);
    ck_assert_int_eq(v->length, 0);
    for (i = 0; i < el_size; i++) {
        append_i_array(ia1, i);
        append_el_i_array_2d(v, 0, i);
    }
    ck_assert_int_eq(v->length, 1);
    for (i = 0; i < el_size; i++) {
        append_i_array(ia2, (el_size + i));
        append_el_i_array_2d(v, 1, (el_size + i));
    }
    ck_assert_int_eq(v->length, 2);
    ret = i_arrays_equal(v->a[0], ia1);
    ck_assert_msg((ret != 0), "`i_array_2d` element %d is incorrect: %d %d %d",
            0, v->a[0]->a[0], v->a[0]->a[1], v->a[0]->a[2]);
    ret = i_arrays_equal(v->a[1], ia2);
    ck_assert_msg((ret != 0), "`i_array_2d` element %d is incorrect: %d %d %d",
            0, v->a[1]->a[0], v->a[1]->a[1], v->a[1]->a[2]);

    free_i_array_2d(v);
    free_i_array(ia1);
    free_i_array(ia2);
}
END_TEST

START_TEST (test_set_i_array_2d) {
    int size, el_size, i, ret;
    i_array_2d * v;
    i_array * ia1;
    i_array * ia2;
    el_size = 3;
    ia1 = init_i_array(el_size);
    ia2 = init_i_array(el_size);
    for (i = 0; i < el_size; i++) {
        append_i_array(ia1, i);
        append_i_array(ia2, (el_size + i));
    }
    size = 1;
    v = init_i_array_2d(size, el_size);
    append_i_array_2d(v, ia1);
    ret = i_arrays_equal(v->a[0], ia1);
    ck_assert_msg((ret != 0), "`i_array_2d` element %d is incorrect: %d %d %d",
            0, v->a[0]->a[0], v->a[0]->a[1], v->a[0]->a[2]);
    set_i_array_2d(v, 0, ia2);
    ret = i_arrays_equal(v->a[0], ia2);
    ck_assert_msg((ret != 0), "`i_array_2d` element %d is incorrect: %d %d %d",
            0, v->a[0]->a[0], v->a[0]->a[1], v->a[0]->a[2]);

    free_i_array_2d(v);
    free_i_array(ia1);
    free_i_array(ia2);
}
END_TEST

START_TEST (test_set_el_i_array_2d) {
    int size, el_size, i, ret;
    i_array_2d * v;
    i_array * ia1;
    i_array * ia2;
    el_size = 3;
    ia1 = init_i_array(el_size);
    ia2 = init_i_array(el_size);
    for (i = 0; i < el_size; i++) {
        append_i_array(ia1, i);
        append_i_array(ia2, i);
    }
    set_i_array(ia2, 1, 99);
    size = 1;
    v = init_i_array_2d(size, el_size);
    append_i_array_2d(v, ia1);
    set_el_i_array_2d(v, 0, 1, 99);
    ret = i_arrays_equal(v->a[0], ia2);
    ck_assert_msg((ret != 0), "`i_array_2d` element is incorrect");
    ret = i_arrays_equal(v->a[0], ia1);
    ck_assert_msg((ret == 0), "`i_array_2d` element is incorrect");

    free_i_array_2d(v);
    free_i_array(ia1);
    free_i_array(ia2);
}
END_TEST

START_TEST (test_set_i_array_2d_fail) {
    int size, el_size, i, ret;
    i_array_2d * v;
    i_array * ia1;
    el_size = 3;
    ia1 = init_i_array(el_size);
    for (i = 0; i < el_size; i++) {
        append_i_array(ia1, i);
    }
    size = 1;
    v = init_i_array_2d(size, el_size);
    append_i_array_2d(v, ia1);
    set_i_array_2d(v, 1, ia1); // SIGABRT
}
END_TEST

START_TEST (test_get_i_array_2d) {
    int size, el_size, i, ret;
    i_array_2d * v;
    i_array * ia1;
    i_array * ia2;
    el_size = 3;
    ia1 = init_i_array(el_size);
    ia2 = init_i_array(el_size);
    for (i = 0; i < el_size; i++) {
        append_i_array(ia1, i);
        append_i_array(ia2, (el_size + i));
    }
    size = 1;
    v = init_i_array_2d(size, el_size);
    append_i_array_2d(v, ia1);
    ret = i_arrays_equal(get_i_array_2d(v, 0), ia1);
    ck_assert_msg((ret != 0), "`i_array_2d` element is incorrect");
    append_i_array_2d(v, ia2);
    ret = i_arrays_equal(get_i_array_2d(v, 1), ia2);
    ck_assert_msg((ret != 0), "`i_array_2d` element is incorrect");

    free_i_array_2d(v);
    free_i_array(ia1);
    free_i_array(ia2);
}
END_TEST

START_TEST (test_get_el_i_array_2d) {
    int size, el_size, i, ret;
    i_array_2d * v;
    i_array * ia1;
    i_array * ia2;
    el_size = 3;
    ia1 = init_i_array(el_size);
    ia2 = init_i_array(el_size);
    for (i = 0; i < el_size; i++) {
        append_i_array(ia1, i);
        append_i_array(ia2, (el_size + i));
    }
    size = 1;
    v = init_i_array_2d(size, el_size);
    append_i_array_2d(v, ia1);
    ck_assert_msg((get_el_i_array_2d(v, 0, 1) == 1), "`i_array_2d` array %d "
            "element %d is %d, expecting %d", 0, 1,
            get_el_i_array_2d(v, 0, 1), 1);
    append_i_array_2d(v, ia2);
    ck_assert_msg((get_el_i_array_2d(v, 1, 1) == 4), "`i_array_2d` array %d "
            "element %d is %d, expecting %d", 1, 1,
            get_el_i_array_2d(v, 1, 1), 4);

    free_i_array_2d(v);
    free_i_array(ia1);
    free_i_array(ia2);
}
END_TEST

START_TEST (test_get_i_array_2d_fail) {
    int size, el_size, i, ret;
    i_array_2d * v;
    i_array * ia1;
    el_size = 3;
    ia1 = init_i_array(el_size);
    for (i = 0; i < el_size; i++) {
        append_i_array(ia1, i);
    }
    size = 1;
    v = init_i_array_2d(size, el_size);
    append_i_array_2d(v, ia1);
    get_i_array_2d(v, 1); // SIGABRT
}
END_TEST

START_TEST (test_split_str) {
    int ret;
    s_array * words;
    s_array * exp;
    words = init_s_array(1);
    exp = init_s_array(1);
    char * string = "one.1\ttwo.2\tthree.3\0";
    append_s_array(exp, "one.1");
    append_s_array(exp, "two.2");
    append_s_array(exp, "three.3");
    ret = split_str(string, words, 3);
    ck_assert_int_eq(ret, 0);
    ck_assert_int_eq(words->length, 3);
    ck_assert_msg((s_arrays_equal(words, exp) != 0),
            "unexpected result of `split_str`");
    ret = split_str(string, words, 3);
    ck_assert_int_eq(ret, 0);
    ret = split_str(string, words, 4);
    ck_assert_int_eq(ret, 3);
    ck_assert_msg((s_arrays_equal(words, exp) != 0),
            "unexpected result of `split_str`");
    free_s_array(words);
    free_s_array(exp);
}
END_TEST

START_TEST (test_split_str_d) {
    int ret;
    d_array * v;
    d_array * exp;
    v = init_d_array(1);
    exp = init_d_array(1);
    char * string = "1.01\t0.009\t2.0e2";
    append_d_array(exp, 1.01);
    append_d_array(exp, 0.009);
    append_d_array(exp, 200.0);
    ret = split_str_d(string, v, 0);
    ck_assert_int_eq(ret, 0);
    ck_assert_int_eq(v->length, 3);
    ck_assert_msg((d_arrays_equal(v, exp, 0.000001) != 0),
            "unexpected result of `split_str_d`");
    ret = split_str_d(string, v, 3);
    ck_assert_int_eq(ret, 0);
    ret = split_str_d(string, v, 4);
    ck_assert_int_eq(ret, 3);
    ck_assert_msg((d_arrays_equal(v, exp, 0.000001) != 0),
            "unexpected result of `split_str_d`");
    free_d_array(v);
    free_d_array(exp);
}
END_TEST

START_TEST (test_split_str_i) {
    int ret;
    i_array * v;
    i_array * exp;
    v = init_i_array(1);
    exp = init_i_array(1);
    char * string = "1\t22\t-432";
    append_i_array(exp, 1);
    append_i_array(exp, 22);
    append_i_array(exp, -432);
    ret = split_str_i(string, v, 0);
    ck_assert_int_eq(ret, 0);
    ck_assert_int_eq(v->length, 3);
    ck_assert_msg((i_arrays_equal(v, exp) != 0),
            "unexpected result of `split_str_i`");
    ret = split_str_i(string, v, 3);
    ck_assert_int_eq(ret, 0);
    ret = split_str_i(string, v, 4);
    ck_assert_int_eq(ret, 3);
    ck_assert_msg((i_arrays_equal(v, exp) != 0),
            "unexpected result of `split_str_i`");
    free_i_array(v);
    free_i_array(exp);
}
END_TEST


Suite * array_utils_suite(void) {
    Suite * s = suite_create("array_utils");

    TCase * tc_d_array = tcase_create("d_array_test_case");
    tcase_add_test(tc_d_array, test_init_free_d_array);
    tcase_add_test_raise_signal(tc_d_array, test_init_d_array_fail, SIGABRT);
    tcase_add_test(tc_d_array, test_expand_d_array);
    tcase_add_test(tc_d_array, test_append_d_array);
    tcase_add_test(tc_d_array, test_set_d_array);
    tcase_add_test_raise_signal(tc_d_array, test_set_d_array_fail, SIGABRT);
    tcase_add_test(tc_d_array, test_get_d_array);
    tcase_add_test_raise_signal(tc_d_array, test_get_d_array_fail, SIGABRT);
    suite_add_tcase(s, tc_d_array);

    TCase * tc_i_array = tcase_create("i_array_test_case");
    tcase_add_test(tc_i_array, test_init_free_i_array);
    tcase_add_test_raise_signal(tc_i_array, test_init_i_array_fail, SIGABRT);
    tcase_add_test(tc_i_array, test_expand_i_array);
    tcase_add_test(tc_i_array, test_append_i_array);
    tcase_add_test(tc_i_array, test_set_i_array);
    tcase_add_test_raise_signal(tc_i_array, test_set_i_array_fail, SIGABRT);
    tcase_add_test(tc_i_array, test_get_i_array);
    tcase_add_test_raise_signal(tc_i_array, test_get_i_array_fail, SIGABRT);
    suite_add_tcase(s, tc_i_array);
    
    TCase * tc_c_array = tcase_create("c_array_test_case");
    tcase_add_test(tc_c_array, test_init_free_c_array);
    tcase_add_test_raise_signal(tc_c_array, test_init_c_array_fail, SIGABRT);
    tcase_add_test(tc_c_array, test_expand_c_array);
    tcase_add_test(tc_c_array, test_assign_c_array);
    tcase_add_test(tc_c_array, test_get_c_array);
    suite_add_tcase(s, tc_c_array);

    TCase * tc_s_array = tcase_create("s_array_test_case");
    tcase_add_test(tc_s_array, test_init_free_s_array);
    tcase_add_test_raise_signal(tc_s_array, test_init_s_array_fail, SIGABRT);
    tcase_add_test(tc_s_array, test_expand_s_array);
    tcase_add_test(tc_s_array, test_append_s_array);
    tcase_add_test(tc_s_array, test_set_s_array);
    tcase_add_test_raise_signal(tc_s_array, test_set_s_array_fail, SIGABRT);
    tcase_add_test(tc_s_array, test_get_s_array);
    tcase_add_test_raise_signal(tc_s_array, test_get_s_array_fail, SIGABRT);
    suite_add_tcase(s, tc_s_array);

    TCase * tc_i_array_2d = tcase_create("i_array_2d_test_case");
    tcase_add_test(tc_i_array_2d, test_init_free_i_array_2d);
    tcase_add_test_raise_signal(tc_i_array_2d, test_init_i_array_2d_fail,
            SIGABRT);
    tcase_add_test(tc_i_array_2d, test_expand_i_array_2d);
    tcase_add_test(tc_i_array_2d, test_append_i_array_2d);
    tcase_add_test(tc_i_array_2d, test_append_el_i_array_2d);
    tcase_add_test(tc_i_array_2d, test_set_i_array_2d);
    tcase_add_test(tc_i_array_2d, test_set_el_i_array_2d);
    tcase_add_test_raise_signal(tc_i_array_2d, test_set_i_array_2d_fail,
            SIGABRT);
    tcase_add_test(tc_i_array_2d, test_get_i_array_2d);
    tcase_add_test(tc_i_array_2d, test_get_el_i_array_2d);
    tcase_add_test_raise_signal(tc_i_array_2d, test_get_i_array_2d_fail,
            SIGABRT);
    suite_add_tcase(s, tc_i_array_2d);

    TCase * tc_arrays_equal = tcase_create("arrays_equal_test_case");
    tcase_add_test(tc_arrays_equal, test_s_arrays_equal);
    tcase_add_test(tc_arrays_equal, test_i_arrays_equal);
    tcase_add_test(tc_arrays_equal, test_d_arrays_equal);
    suite_add_tcase(s, tc_arrays_equal);

    TCase * tc_get_matching_indices = tcase_create(
            "get_matching_indices_test_case");
    tcase_add_test(tc_get_matching_indices, test_get_matching_indices);
    tcase_add_exit_test(tc_get_matching_indices, test_get_matching_indices_fail,
            1);
    suite_add_tcase(s, tc_get_matching_indices);

    TCase * tc_get_doubles = tcase_create("get_doubles_test_case");
    tcase_add_test(tc_get_doubles, test_get_doubles);
    suite_add_tcase(s, tc_get_doubles);

    TCase * tc_split_str = tcase_create("split_str_test_case");
    tcase_add_test(tc_split_str, test_split_str);
    tcase_add_test(tc_split_str, test_split_str_d);
    tcase_add_test(tc_split_str, test_split_str_i);
    suite_add_tcase(s, tc_split_str);

    return s;
}

int main(void) {
    int number_failed;
    Suite * s = array_utils_suite();
    SRunner * sr = srunner_create(s);
    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

