set(MY_TARGET main)

file(GLOB SOURCES "src/*.cpp" "src/*.c")
file(GLOB HEADERS "include/*.hpp" "include/*.h")

add_executable(${MY_TARGET} ${SOURCES} ${HEADERS})

target_link_libraries(${MY_TARGET} qboot)
target_include_directories(${MY_TARGET} PUBLIC
    ../algebra/include
    ../qboot/include
)

set_target_properties(
    ${MY_TARGET} PROPERTIES
    OUTPUT_NAME_DEBUG ${MY_TARGET}-debug
    OUTPUT_NAME_RELEASE ${MY_TARGET}
)
