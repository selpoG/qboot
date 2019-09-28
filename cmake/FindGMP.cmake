include(cmake/utils.cmake)  # for my_find_include, my_find_lib

if(MSVC)
	my_find_include(GMP_INCLUDE_DIR gmp.h ${GMP_ROOT}/Debug)
	my_find_lib(GMP_LIBRARY_DEBUG mpir ${GMP_ROOT}/Debug)
	my_find_lib(GMP_LIBRARY_RELEASE mpir ${GMP_ROOT}/Release)
	my_find_lib(GMPXX_LIBRARY_DEBUG mpirxx ${GMP_ROOT}/Debug)
	my_find_lib(GMPXX_LIBRARY_RELEASE mpirxx ${GMP_ROOT}/Release)
	set(GMP_LIBRARY debug ${GMP_LIBRARY_DEBUG} optimized ${GMP_LIBRARY_RELEASE})
	set(GMPXX_LIBRARY debug ${GMPXX_LIBRARY_DEBUG} optimized ${GMPXX_LIBRARY_RELEASE})
	mark_as_advanced(GMP_LIBRARY_DEBUG GMPXX_LIBRARY_DEBUG)
	mark_as_advanced(GMP_LIBRARY_RELEASE GMPXX_LIBRARY_RELEASE)
else()
	my_find_include(GMP_INCLUDE_DIR gmp.h ${GMP_ROOT})
	my_find_lib(GMP_LIBRARY gmp ${GMP_ROOT})
	my_find_lib(GMPXX_LIBRARY gmpxx ${GMP_ROOT})
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY GMPXX_LIBRARY)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(GMP
	REQUIRED_VARS
		GMP_LIBRARY
		GMPXX_LIBRARY
		GMP_INCLUDE_DIR
)
