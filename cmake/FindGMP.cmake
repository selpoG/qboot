include(cmake/utils.cmake)  # for my_find_include, my_find_lib

if(MSVC)
	set(GMP_ROOT_DEBUG "${GMP_ROOT}/Debug")
	set(GMP_ROOT_RELEASE "${GMP_ROOT}/Release")
else()
	set(GMP_ROOT_DEBUG "${GMP_ROOT}")
	set(GMP_ROOT_RELEASE "${GMP_ROOT}")
endif()

my_find_include(GMP_INCLUDE_DIR gmp.h ${GMP_ROOT_DEBUG})

if(MSVC)
	set(_GMP_LIB_NAME mpir)
else()
	set(_GMP_LIB_NAME gmp)
endif()

my_find_lib(GMP_LIBRARY_DEBUG ${_GMP_LIB_NAME} ${GMP_ROOT_DEBUG})
my_find_lib(GMP_LIBRARY_RELEASE ${_GMP_LIB_NAME} ${GMP_ROOT_RELEASE})
my_find_lib(GMPXX_LIBRARY_DEBUG ${_GMP_LIB_NAME}xx ${GMP_ROOT_DEBUG})
my_find_lib(GMPXX_LIBRARY_RELEASE ${_GMP_LIB_NAME}xx ${GMP_ROOT_RELEASE})

set(GMP_LIBRARY debug ${GMP_LIBRARY_DEBUG} optimized ${GMP_LIBRARY_RELEASE})
set(GMPXX_LIBRARY debug ${GMPXX_LIBRARY_DEBUG} optimized ${GMPXX_LIBRARY_RELEASE})

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY GMPXX_LIBRARY)
mark_as_advanced(GMP_LIBRARY_DEBUG GMPXX_LIBRARY_DEBUG)
mark_as_advanced(GMP_LIBRARY_RELEASE GMPXX_LIBRARY_RELEASE)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(GMP
	REQUIRED_VARS
		GMP_INCLUDE_DIR
		GMP_LIBRARY
		GMPXX_LIBRARY
		GMP_LIBRARY_DEBUG
		GMP_LIBRARY_RELEASE
)
