include(cmake/utils.cmake)  # for my_find_include, my_find_lib

if(MSVC)
	set(MPFR_ROOT_DEBUG "${MPFR_ROOT}/Debug")
	set(MPFR_ROOT_RELEASE "${MPFR_ROOT}/Release")
else()
	set(MPFR_ROOT_DEBUG "${MPFR_ROOT}")
	set(MPFR_ROOT_RELEASE "${MPFR_ROOT}")
endif()

my_find_include(MPFR_INCLUDE_DIR mpfr.h ${MPFR_ROOT_DEBUG})

my_find_lib(MPFR_LIBRARY_DEBUG mpfr ${MPFR_ROOT_DEBUG})
my_find_lib(MPFR_LIBRARY_RELEASE mpfr ${MPFR_ROOT_RELEASE})

set(MPFR_LIBRARY debug ${MPFR_LIBRARY_DEBUG} optimized ${MPFR_LIBRARY_RELEASE})

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)
mark_as_advanced(MPFR_LIBRARY_DEBUG MPFR_LIBRARY_RELEASE)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(MPFR
	REQUIRED_VARS
		MPFR_INCLUDE_DIR
		MPFR_LIBRARY
		MPFR_LIBRARY_DEBUG
		MPFR_LIBRARY_RELEASE
)
