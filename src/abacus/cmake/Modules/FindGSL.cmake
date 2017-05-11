# Find GSL libraries
# Once done this will define
#
# GSL_FOUND - System has gsl
# GSL_CBLAS_FOUND - System has gsl_cblas library
# GSL_INCLUDE_DIRS - The gsl include directories
# GSL_LIBRARIES - The libraries needed to use gsl
# GSL_DEFINITIONS - Compiler switches required for using gsl

include (FindPkgConfig)
if (PKG_CONFIG_FOUND)
    if (GSL_FIND_VERSION)
        set (_PACKAGE_ARGS "gsl>=${GSL_FIND_VERSION}")
    else ()
        set (_PACKAGE_ARGS "gsl")
    endif (GSL_FIND_VERSION)
    if (GSL_FIND_REQUIRED)
        set(_PACKAGE_ARGS "${_PACKAGE_ARGS}" REQUIRED)
    endif (GSL_FIND_REQUIRED)
    pkg_check_modules (PC_GSL "${_PACKAGE_ARGS}")
    set (GSL_DEFINITIONS "${PC_GSL_CFLAGS_OTHER}")
    set (GSL_VERSION_STRING "${PC_GSL_VERSION}")
endif (PKG_CONFIG_FOUND)

find_path (GSL_INCLUDE_DIR
    "gsl/gsl_randist.h"
    HINTS
    "${PC_GSL_INCUDEDIR}"
    "${PC_GSL_INCLUDE_DIRS}"
    )
find_library (GSL_LIBRARY
    NAMES "gsl"
    HINTS
    "${PC_GSL_LIBDIR}"
    "${PC_GSL_LIBRARY_DIRS}"
    )
find_library (GSL_STATIC_LIBRARY
    NAMES "libgsl.a"
    HINTS
    "${PC_GSL_LIBDIR}"
    "${PC_GSL_LIBRARY_DIRS}"
    )
find_library (GSL_CBLAS_LIBRARY
    NAMES "gslcblas"
    HINTS
    "${PC_GSL_LIBDIR}"
    "${PC_GSL_LIBRARY_DIRS}"
    )
find_library (GSL_CBLAS_STATIC_LIBRARY
    NAMES "libgslcblas.a"
    HINTS
    "${PC_GSL_LIBDIR}"
    "${PC_GSL_LIBRARY_DIRS}"
    )

mark_as_advanced (GSL_INCLUDE_DIR GSL_LIBRARY GSL_CBLAS_LIBRARY GSL_STATIC_LIBRARY GSL_CBLAS_STATIC_LIBRARY)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (gsl
    DEFAULT_MSG
    GSL_LIBRARY
    GSL_INCLUDE_DIR
    )

find_package_handle_standard_args (gsl_cblas
    DEFAULT_MSG
    GSL_CBLAS_LIBRARY
    )

if (GSL_FOUND)
    set (GSL_LIBRARIES ${GSL_LIBRARY})
    set (GSL_STATIC_LIBRARIES ${GSL_STATIC_LIBRARY})
    set (GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR})
endif (GSL_FOUND)

if (GSL_CBLAS_FOUND)
    set (GSL_LIBRARIES ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARY})
    set (GSL_STATIC_LIBRARIES ${GSL_STATIC_LIBRARIES} ${GSL_CBLAS_STATIC_LIBRARY})
endif (GSL_CBLAS_FOUND)

message(STATUS "GSL_FOUND: ${GSL_FOUND}")
message(STATUS "GSL_CBLAS_FOUND: ${GSL_CBLAS_FOUND}")
message(STATUS "GSL_LIBRARIES: ${GSL_LIBRARIES}")
message(STATUS "GSL_STATIC_LIBRARIES: ${GSL_STATIC_LIBRARIES}")
message(STATUS "GSL_INCLUDE_DIRS: ${GSL_INCLUDE_DIRS}")
message(STATUS "GSL_DEFINITIONS: ${GSL_DEFINITIONS}")

