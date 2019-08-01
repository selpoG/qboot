set(MY_TARGET util)

file(GLOB SOURCES "src/*.cpp" "src/*.c")
file(GLOB HEADERS "include/*.hpp" "include/*.h")

add_library(${MY_TARGET} STATIC ${SOURCES} ${HEADERS})

target_include_directories(${MY_TARGET} PUBLIC include)

set_target_properties(
	${MY_TARGET} PROPERTIES
	OUTPUT_NAME_DEBUG ${MY_TARGET}-debug
	OUTPUT_NAME_RELEASE ${MY_TARGET}
)
