find_path(GMP_INCLUDE_DIR_DEBUG
	NAMES gmp.h
	PATHS
		ENV GMP_ROOT_DEBUG
		ENV GMP_INCLUDE_DIR_DEBUG
		${GMP_ROOT_DEBUG}
		/usr
		/usr/local
	PATH_SUFFIXES include
)

find_path(GMP_INCLUDE_DIR_RELEASE
	NAMES gmp.h
	PATHS
		ENV GMP_ROOT_RELEASE
		ENV GMP_INCLUDE_DIR_RELEASE
		${GMP_ROOT_RELEASE}
		/usr
		/usr/local
	PATH_SUFFIXES include
)

if(MSVC)
	set(_GMP_LIB_NAME mpir)
else()
	set(_GMP_LIB_NAME gmp)
endif()

find_library(GMP_LIBRARY_DEBUG
	NAMES ${_GMP_LIB_NAME}
	PATHS
		ENV GMP_ROOT_DEBUG
		ENV GMP_LIB_DIR_DEBUG
		${GMP_ROOT_DEBUG}
		/usr
		/usr/local
	PATH_SUFFIXES lib
)

find_library(GMP_LIBRARY_RELEASE
	NAMES ${_GMP_LIB_NAME}
	PATHS
		ENV GMP_ROOT_RELEASE
		ENV GMP_LIB_DIR_RELEASE
		${GMP_ROOT_RELEASE}
		/usr
		/usr/local
	PATH_SUFFIXES lib
)

find_library(GMPXX_LIBRARY_DEBUG
	NAMES ${_GMP_LIB_NAME}xx
	PATHS
		ENV GMP_ROOT_DEBUG
		ENV GMP_LIB_DIR_DEBUG
		${GMP_ROOT_DEBUG}
		/usr
		/usr/local
	PATH_SUFFIXES lib
)

find_library(GMPXX_LIBRARY_RELEASE
	NAMES ${_GMP_LIB_NAME}xx
	PATHS
		ENV GMP_ROOT_RELEASE
		ENV GMP_LIB_DIR_RELEASE
		${GMP_ROOT_RELEASE}
		/usr
		/usr/local
	PATH_SUFFIXES lib
)

mark_as_advanced(GMP_INCLUDE_DIR_DEBUG GMP_LIBRARY_DEBUG GMPXX_LIBRARY_DEBUG)
mark_as_advanced(GMP_INCLUDE_DIR_RELEASE GMP_LIBRARY_RELEASE GMPXX_LIBRARY_RELEASE)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(GMP
	REQUIRED_VARS
		GMP_INCLUDE_DIR_DEBUG
		GMP_INCLUDE_DIR_RELEASE
		GMP_LIBRARY_DEBUG
		GMP_LIBRARY_RELEASE
)
