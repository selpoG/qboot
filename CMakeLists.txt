cmake_minimum_required(VERSION 3.13.0)
project(qboot CXX)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED on)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}" ${CMAKE_MODULE_PATH})

find_package(MPFR REQUIRED)

if(MSVC)
	string(REPLACE "/W3" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	string(REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
	string(REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
	if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.14)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /experimental:external /external:W0")
		if (IS_DIRECTORY "${MPFR_INCLUDE_DIR}")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /external:I ${MPFR_INCLUDE_DIR}")
		endif()
		if (IS_DIRECTORY "${GMP_INCLUDE_DIR}")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /external:I ${GMP_INCLUDE_DIR}")
		endif()
	endif()
else()
	if (IS_DIRECTORY "${MPFR_INCLUDE_DIR}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${MPFR_INCLUDE_DIR}")
	endif()
	if (IS_DIRECTORY "${GMP_INCLUDE_DIR}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${GMP_INCLUDE_DIR}")
	endif()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -Wpedantic -g -Wno-padded")
	if(CMAKE_COMPILER_IS_GNUCXX)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++ -Wno-noexcept -Walloca -Wcast-align -Wcast-qual -Wconditionally-supported -Wconversion -Wctor-dtor-privacy -Wdate-time -Wdisabled-optimization -Wdouble-promotion -Wduplicated-cond -Wfloat-conversion -Wfloat-equal -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wmultichar -Wmultiple-inheritance -Wnull-dereference -Wold-style-cast -Woverlength-strings -Woverloaded-virtual -Wpacked -Wredundant-decls -Wrestrict -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wsuggest-override -Wswitch-default -Wswitch-enum -Wuninitialized -Wunused-macros -Wunused-parameter -Wuseless-cast -Wvariadic-macros -Wvector-operation-performance -Wvirtual-inheritance")
	else()
		if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 9.0.0)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ctad-maybe-unsupported")
		endif()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything -Wno-c++98-compat -Wno-shadow-field-in-constructor -Wno-c++98-compat-pedantic -Wno-global-constructors -Wno-exit-time-destructors -Wno-covered-switch-default")
	endif()
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ../bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ../lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ../lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG ../lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ../lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ../lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ../lib)

add_subdirectory(util)
add_subdirectory(algebra)
add_subdirectory(qboot)
add_subdirectory(main)

set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT main)

add_custom_target(run
	COMMAND bin/main
	DEPENDS main
	WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
