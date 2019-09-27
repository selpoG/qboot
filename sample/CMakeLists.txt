cmake_minimum_required(VERSION 3.13.0)
project(qboot CXX)

set(CMAKE_PREFIX_PATH ${QBoot_DIR} ${CMAKE_PREFIX_PATH})
find_package(QBoot 0.5.0 REQUIRED)

set(MY_TARGET sample)
add_executable(${MY_TARGET} main.cpp)

target_link_libraries(${MY_TARGET} QBoot::qboot)

if(MSVC)
	string(REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
	string(REPLACE "/MD" "/MT" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
else()
	set(cxx_flags -pthread -fPIC)
endif()

set_target_properties(${MY_TARGET} PROPERTIES
	DEBUG_POSTFIX "-debug"
	COMPILE_OPTIONS "${cxx_flags}"
	CXX_EXTENSIONS OFF
	CXX_STANDARD 17
	CXX_STANDARD_REQUIRED ON
	RUNTIME_OUTPUT_DIRECTORY_DEBUG bin
	RUNTIME_OUTPUT_DIRECTORY_RELEASE bin
)
