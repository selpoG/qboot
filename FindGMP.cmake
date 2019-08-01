find_path(GMP_INCLUDE_DIR
	NAMES gmp.h
	PATHS
		ENV GMP_ROOT
		ENV GMP_INCLUDE_DIR
		${GMP_ROOT}
		/usr
		/usr/local
	PATH_SUFFIXES include
)

if(MSVC)
	set(_GMP_LIB_NAME mpir)
else()
	set(_GMP_LIB_NAME gmp)
endif()

find_library(GMP_LIBRARY
	NAMES ${_GMP_LIB_NAME}
	PATHS
		ENV GMP_ROOT
		ENV GMP_LIB_DIR
		${GMP_ROOT}
		/usr
		/usr/local
	PATH_SUFFIXES lib
)

find_library(GMPXX_LIBRARY
	NAMES ${_GMP_LIB_NAME}xx
	PATHS
		ENV GMP_ROOT
		ENV GMP_LIB_DIR
		${GMP_ROOT}
		/usr
		/usr/local
	PATH_SUFFIXES lib
)

mark_as_advanced(GMP_INCLUDE_DIR GMP_LIBRARY)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(GMP
	REQUIRED_VARS
		GMP_INCLUDE_DIR
		GMP_LIBRARY
)

if(GMP_FOUND AND NOT TARGET GMP::GMP)
	add_library(GMP::GMP UNKNOWN IMPORTED)
	set_target_properties(GMP::GMP PROPERTIES
		IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
		IMPORTED_LOCATION "${GMP_LIBRARY}"
		INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
	)
endif()
