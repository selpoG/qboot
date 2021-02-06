cmake_minimum_required(VERSION 3.12.4)
set(QBOOT_MAJOR_VERSION 0)
set(QBOOT_MINOR_VERSION 8)
set(QBOOT_PATCH_VERSION 1)
set(QBOOT_VERSION ${QBOOT_MAJOR_VERSION}.${QBOOT_MINOR_VERSION}.${QBOOT_PATCH_VERSION})
project(qboot LANGUAGES CXX VERSION ${QBOOT_VERSION})

set(CMAKE_CONFIGURATION_TYPES "Debug;Release")

option(qboot_build_tests "Build test targets (for developers)" OFF)
mark_as_advanced(qboot_build_tests)

if(MSVC)
	if(CMAKE_BUILD_TYPE)
		message("CMAKE_BUILD_TYPE is set to `${CMAKE_BUILD_TYPE}`, but would be ignored.")
	endif()
else()
	if(NOT CMAKE_BUILD_TYPE)
		set(CMAKE_BUILD_TYPE "Release")
		message("CMAKE_BUILD_TYPE is not specified: Use Release by default.")
	endif()
	if(NOT((CMAKE_BUILD_TYPE STREQUAL "Debug") OR (CMAKE_BUILD_TYPE STREQUAL "Release")))
		message(FATAL_ERROR "CMAKE_BUILD_TYPE must be 'Debug' or 'Release' (case-sensitive)")
	endif()
endif()

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
include(cmake/utils.cmake)  # for fix_default_cxx_flags, mark_external, set_flags, set_default_property

set(GMP_ROOT "" CACHE PATH "path of GMP. If you have gmp.h in /some/where/include/ and libgmp.a in /some/where/lib/, please set `/some/where`")
set(MPFR_ROOT "" CACHE PATH "path of MPFR. If you have mpfr.h in /some/where/include/ and libmpfr.a in /some/where/lib/, please set `/some/where`")

find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

fix_default_cxx_flags()
mark_external(CMAKE_CXX_FLAGS ${MPFR_INCLUDE_DIR})
mark_external(CMAKE_CXX_FLAGS ${GMP_INCLUDE_DIR})
set_flags()

set(MY_TARGET qboot)

add_library(${MY_TARGET} STATIC
	src/mp/integer.cpp
	src/mp/rational.cpp
	src/mp/real.cpp
	src/algebra/complex_function.cpp
	src/algebra/matrix.cpp
	src/algebra/polynomial.cpp
	src/algebra/real_function.cpp
	src/block.cpp
	src/bootstrap_equation.cpp
	src/conformal_scale.cpp
	src/context.cpp
	src/hor_formula.cpp
	src/hor_recursion.cpp
	src/polynomial_program.cpp
	src/primary_op.cpp
	src/scale_factor.cpp
	src/sdpb_input.cpp
	src/task_queue.cpp
	src/xml_input.cpp
	include/qboot/mp/integer.hpp
	include/qboot/mp/rational.hpp
	include/qboot/mp/real.hpp
	include/qboot/algebra/complex_function.hpp
	include/qboot/algebra/matrix.hpp
	include/qboot/algebra/polynomial.hpp
	include/qboot/algebra/real_function.hpp
	include/qboot/block.hpp
	include/qboot/bootstrap_equation.hpp
	include/qboot/conformal_scale.hpp
	include/qboot/context.hpp
	include/qboot/hor_formula.hpp
	include/qboot/hor_recursion.hpp
	include/qboot/my_concepts.hpp
	include/qboot/my_filesystem.hpp
	include/qboot/polynomial_program.hpp
	include/qboot/primary_op.hpp
	include/qboot/qboot.hpp
	include/qboot/scale_factor.hpp
	include/qboot/sdpb_input.hpp
	include/qboot/task_queue.hpp
	include/qboot/xml_input.hpp
)

find_package(Threads)
target_link_libraries(${MY_TARGET} INTERFACE Threads::Threads)
target_link_libraries(
	${MY_TARGET} INTERFACE
	${MPFR_LIBRARY}
	${GMP_LIBRARY}
	${GMPXX_LIBRARY}
)
if(NOT MSVC)
	target_link_libraries(${MY_TARGET} INTERFACE stdc++fs)
	if(NOT CMAKE_COMPILER_IS_GNUCXX)
		target_link_libraries(${MY_TARGET} INTERFACE c++abi)
	endif()
endif()

target_include_directories(${MY_TARGET} PUBLIC
	$<BUILD_INTERFACE:${MPFR_INCLUDE_DIR}>
	$<BUILD_INTERFACE:${GMP_INCLUDE_DIR}>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	INTERFACE
	$<INSTALL_INTERFACE:${MPFR_INCLUDE_DIR}>
	$<INSTALL_INTERFACE:${GMP_INCLUDE_DIR}>
	$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

set_default_property(${MY_TARGET} .)

include(CMakePackageConfigHelpers)
install(DIRECTORY "include/" DESTINATION include)
set(cmake_package_name QBoot)
set(targets_export_name ${cmake_package_name}Targets CACHE INTERNAL "")
set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated" CACHE INTERNAL "")
set(cmake_files_install_dir "lib/cmake/${cmake_package_name}")
set(version_file "${generated_dir}/${cmake_package_name}ConfigVersion.cmake")
write_basic_package_version_file(${version_file} VERSION ${QBOOT_VERSION} COMPATIBILITY AnyNewerVersion)
install(EXPORT ${targets_export_name}
	NAMESPACE ${cmake_package_name}::
	DESTINATION ${cmake_files_install_dir}
)
set(config_file "${generated_dir}/${cmake_package_name}Config.cmake")
configure_package_config_file("cmake/Config.cmake.in"
	"${config_file}" INSTALL_DESTINATION ${cmake_files_install_dir})
install(FILES ${version_file} ${config_file}
	DESTINATION ${cmake_files_install_dir})
install(TARGETS ${MY_TARGET}
	EXPORT ${targets_export_name}
	CONFIGURATIONS Debug Release
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib)

if(qboot_build_tests)
	add_subdirectory(test)
	set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT test)
else()
	set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT qboot)
endif()
