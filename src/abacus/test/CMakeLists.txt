find_package (Check)

if (NOT CHECK_FOUND)
    message(STATUS "Check library not found; skipping unit tests...")
else (NOT CHECK_FOUND)

    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CHECK_DEFINITIONS}")

    enable_testing()
    include_directories ("${CHECK_INCLUDE_DIRS}")
    include_directories("${PROJECT_SOURCE_DIR}/test" "${PROJECT_SOURCE_DIR}/src")

    find_library (MATH_LIBRARY m)
    find_library (RT_LIBRARY rt)

    list (APPEND C_LIBS "${CHECK_LIBRARIES}")
    if (NOT RT_LIBRARY STREQUAL "RT_LIBRARY-NOTFOUND")
        list (APPEND C_LIBS "${RT_LIBRARY}")
    endif ()
    if (NOT MATH_LIBRARY STREQUAL "MATH_LIBRARY-NOTFOUND")
        list (APPEND C_LIBS "${MATH_LIBRARY}")
    endif ()
    message(STATUS "Check libraries: ${C_LIBS}")
    
    add_executable (check_array_utils EXCLUDE_FROM_ALL
        check_array_utils.c
        test_utils.c
        test_utils.h
        )
    target_link_libraries(check_array_utils
        "${C_LIBS}"
        )
    add_test(check_array_utils "${CMAKE_CURRENT_BINARY_DIR}/check_array_utils")
    add_dependencies (check check_array_utils)

    add_executable (check_stats_utils EXCLUDE_FROM_ALL
        check_stats_utils.c
        ${PROJECT_SOURCE_DIR}/src/array_utils.c
        test_utils.c
        test_utils.h
        )
    target_link_libraries(check_stats_utils
        "${C_LIBS}"
        )
    add_test(check_stats_utils "${CMAKE_CURRENT_BINARY_DIR}/check_stats_utils")
    add_dependencies (check check_stats_utils)

    add_executable (check_math_utils EXCLUDE_FROM_ALL
        check_math_utils.c
        ${PROJECT_SOURCE_DIR}/src/array_utils.c
        test_utils.c
        test_utils.h
        )
    target_link_libraries(check_math_utils
        "${C_LIBS}"
        )
    add_test(check_math_utils "${CMAKE_CURRENT_BINARY_DIR}/check_math_utils")
    add_dependencies (check check_math_utils)

    add_executable (check_parsing EXCLUDE_FROM_ALL
        check_parsing.c
        ${PROJECT_SOURCE_DIR}/src/array_utils.c
        ${PROJECT_SOURCE_DIR}/src/stats_utils.c
        test_utils.c
        test_utils.h
        )
    target_link_libraries(check_parsing
        "${C_LIBS}"
        )
    add_test(check_parsing "${CMAKE_CURRENT_BINARY_DIR}/check_parsing")
    add_dependencies (check check_parsing)

    add_executable (check_partition_combinatorics EXCLUDE_FROM_ALL
        check_partition_combinatorics.c
        ${PROJECT_SOURCE_DIR}/src/array_utils.c
        test_utils.c
        test_utils.h
        )
    target_link_libraries(check_partition_combinatorics
        "${C_LIBS}"
        )
    add_test(check_partition_combinatorics "${CMAKE_CURRENT_BINARY_DIR}/check_partition_combinatorics")
    add_dependencies (check check_partition_combinatorics)

    add_executable (check_eureject EXCLUDE_FROM_ALL
        check_eureject.c
        ${PROJECT_SOURCE_DIR}/src/array_utils.c
        ${PROJECT_SOURCE_DIR}/src/stats_utils.c
        ${PROJECT_SOURCE_DIR}/src/math_utils.c
        ${PROJECT_SOURCE_DIR}/src/parsing.c
        ${PROJECT_SOURCE_DIR}/src/abacus.c
        test_utils.c
        test_utils.h
        )
    target_link_libraries(check_eureject
        "${C_LIBS}"
        )
    add_test(check_eureject "${CMAKE_CURRENT_BINARY_DIR}/check_eureject")
    add_dependencies (check check_eureject)

    find_package(GSL)

    if (GSL_FOUND)
        include_directories ("${GSL_INCLUDE_DIRS}")

        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GSL_DEFINITIONS}")

        add_executable (check_partition_combinatorics_random EXCLUDE_FROM_ALL
            check_partition_combinatorics_random.c
            ${PROJECT_SOURCE_DIR}/src/array_utils.c
            ${PROJECT_SOURCE_DIR}/src/partition_combinatorics.c
	        test_rng.c
	        test_rng.h
            )
        target_link_libraries(check_partition_combinatorics_random
            "${GSL_LIBRARIES}"
            "${C_LIBS}"
            )
        add_test(check_partition_combinatorics_random "${CMAKE_CURRENT_BINARY_DIR}/check_partition_combinatorics_random")
        add_dependencies (check check_partition_combinatorics_random)

        add_executable (check_probability EXCLUDE_FROM_ALL
            check_probability.c
            ${PROJECT_SOURCE_DIR}/src/array_utils.c
	        test_rng.c
	        test_rng.h
            )
        target_link_libraries(check_probability
            "${GSL_LIBRARIES}"
            "${C_LIBS}"
            )
        add_test(check_probability "${CMAKE_CURRENT_BINARY_DIR}/check_probability")
        add_dependencies (check check_probability)

    endif (GSL_FOUND)

endif (NOT CHECK_FOUND)

