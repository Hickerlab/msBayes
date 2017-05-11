find_package (GSL)

include_directories("${PROJECT_SOURCE_DIR}/test" "${PROJECT_SOURCE_DIR}/src")

find_library (MATH_LIBRARY m)
find_library (MATH_STATIC_LIBRARY libm.a)
if ("${MATH_STATIC_LIBRARY}" STREQUAL "MATH_STATIC_LIBRARY-NOTFOUND")
    set (MATH_STATIC_LIBRARY "${MATH_LIBRARY}")
endif ("${MATH_STATIC_LIBRARY}" STREQUAL "MATH_STATIC_LIBRARY-NOTFOUND")

if (STATIC_LINKING)
    set (CMAKE_FIND_LIBRARY_SUFFIXES ".a" ${CMAKE_FIND_LIBRARY_SUFFIXES})
    if (APPLE)
        set (CMAKE_EXE_LINKER_FLAGS "-Bstatic")
    else (APPLE)
        set (CMAKE_EXE_LINKER_FLAGS "-static")
    endif (APPLE)
    set (CMAKE_SKIP_RPATH YES)
    set (M_LIB ${MATH_STATIC_LIBRARY})
    set (G_LIBS ${GSL_STATIC_LIBRARIES})
else (STATIC_LINKING)
    set (M_LIB ${MATH_LIBRARY})
    set (G_LIBS ${GSL_LIBRARIES})
endif (STATIC_LINKING)

add_executable (eureject
	eureject.c
	eureject.h
	eureject_main.c
	array_utils.c
	array_utils.h
	math_utils.c
	math_utils.h
	stats_utils.c
	stats_utils.h
	parsing.c
	parsing.h
	abacus.c
	abacus.h
    )
target_link_libraries(eureject
    "${M_LIB}"
    )

add_executable(intpart
	intpart.c
	intpart.h
	intpart_main.c
	array_utils.c
	array_utils.h
	partition_combinatorics.c
	partition_combinatorics.h
	abacus.c
	abacus.h
    )
target_link_libraries(intpart
    "${M_LIB}"
    )

if (GSL_FOUND)
    include_directories ("${GSL_INCLUDE_DIRS}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${GSL_DEFINITIONS}")

    add_executable(dpdraw
	    dpdraw.c
	    dpdraw.h
	    dpdraw_main.c
	    array_utils.c
	    array_utils.h
	    partition_combinatorics.c
	    partition_combinatorics.h
	    partition_combinatorics_random.c
	    partition_combinatorics_random.h
	    abacus.c
	    abacus.h
        )
    target_link_libraries(dpdraw
        "${G_LIBS}"
        "${M_LIB}"
        )
    install(TARGETS dpdraw DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")
    if (STATIC_LINKING)
        set_target_properties(dpdraw PROPERTIES
            LINK_SEARCH_END_STATIC ON
            )
    endif (STATIC_LINKING)
endif (GSL_FOUND)

if (STATIC_LINKING)
    set_target_properties(eureject intpart PROPERTIES
        LINK_SEARCH_END_STATIC ON
        )
endif (STATIC_LINKING)

install(TARGETS eureject intpart DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

