cmake_minimum_required(VERSION 3.12.4)
project(sample CXX)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
find_package(QBoot 0.5.0 REQUIRED)

add_executable(sample main.cpp)

target_link_libraries(sample QBoot::qboot)

if(MSVC)
	string(REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
	string(REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
else()
	set(cxx_flags -pthread -fPIC)
endif()

set_target_properties(sample PROPERTIES
	DEBUG_POSTFIX "-debug"
	COMPILE_OPTIONS "${cxx_flags}"
	CXX_EXTENSIONS OFF
	CXX_STANDARD 17
	CXX_STANDARD_REQUIRED ON
	RUNTIME_OUTPUT_DIRECTORY bin
	RUNTIME_OUTPUT_DIRECTORY_DEBUG bin
	RUNTIME_OUTPUT_DIRECTORY_RELEASE bin
)

set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT sample)
