# Find Check libraries
# Once done this will define
#
# CHECK_FOUND - System has check
# CHECK_INCLUDE_DIRS - The check include directories
# CHECK_LIBRARIES - The libraries needed to use check
# CHECK_DEFINITIONS - Compiler switches required for using check

include (FindPkgConfig)
if (PKG_CONFIG_FOUND)
    if (CHECK_FIND_VERSION)
        set (_PACKAGE_ARGS "check>=${CHECK_FIND_VERSION}")
    else ()
        set (_PACKAGE_ARGS "check")
    endif (CHECK_FIND_VERSION)
    if (CHECK_FIND_REQUIRED)
        set(_PACKAGE_ARGS "${_PACKAGE_ARGS}" REQUIRED)
    endif (CHECK_FIND_REQUIRED)
    pkg_check_modules (PC_CHECK "${_PACKAGE_ARGS}")
    set (CHECK_DEFINITIONS "${PC_CHECK_CFLAGS_OTHER}")
    set (CHECK_VERSION_STRING "${PC_CHECK_VERSION}")
endif (PKG_CONFIG_FOUND)

find_path (CHECK_INCLUDE_DIR
    "check.h"
    HINTS
    "${PC_CHECK_INCUDEDIR}"
    "${PC_CHECK_INCLUDE_DIRS}"
    )
find_library (CHECK_LIBRARY
    NAMES "check"
    HINTS
    "${PC_CHECK_LIBDIR}"
    "${PC_CHECK_LIBRARY_DIRS}"
    )

mark_as_advanced (CHECK_INCLUDE_DIR CHECK_LIBRARY)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (check
    DEFAULT_MSG
    CHECK_LIBRARY
    CHECK_INCLUDE_DIR
    )

if (CHECK_FOUND)
    set (CHECK_LIBRARIES ${CHECK_LIBRARY})
    set (CHECK_INCLUDE_DIRS ${CHECK_INCLUDE_DIR})
endif ()

message(STATUS "CHECK_FOUND: ${CHECK_FOUND}")
message(STATUS "CHECK_LIBRARIES: ${CHECK_LIBRARIES}")
message(STATUS "CHECK_INCLUDE_DIRS: ${CHECK_INCLUDE_DIRS}")
message(STATUS "CHECK_DEFINITIONS: ${CHECK_DEFINITIONS}")
