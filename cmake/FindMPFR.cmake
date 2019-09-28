include(cmake/utils.cmake)  # for my_find_include, my_find_lib

if(MSVC)
	my_find_include(MPFR_INCLUDE_DIR mpfr.h ${MPFR_ROOT}/Debug)
	my_find_lib(MPFR_LIBRARY_DEBUG mpfr ${MPFR_ROOT}/Debug)
	my_find_lib(MPFR_LIBRARY_RELEASE mpfr ${MPFR_ROOT}/Release)
	set(MPFR_LIBRARY debug ${MPFR_LIBRARY_DEBUG} optimized ${MPFR_LIBRARY_RELEASE})
	mark_as_advanced(MPFR_LIBRARY_DEBUG MPFR_LIBRARY_RELEASE)
else()
	my_find_include(MPFR_INCLUDE_DIR mpfr.h ${MPFR_ROOT})
	my_find_lib(MPFR_LIBRARY mpfr ${MPFR_ROOT})
endif()

mark_as_advanced(MPFR_INCLUDE_DIR MPFR_LIBRARY)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(MPFR
	REQUIRED_VARS
		MPFR_LIBRARY
		MPFR_INCLUDE_DIR
)
