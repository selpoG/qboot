set(MY_TARGET test)

add_executable(${MY_TARGET} src/main.cpp)

target_link_libraries(${MY_TARGET} PRIVATE qboot)
target_include_directories(${MY_TARGET} PUBLIC
    ../include
)

set_default_property(${MY_TARGET} ..)
