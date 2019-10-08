set(MY_TARGET test)
add_executable(${MY_TARGET} src/main.cpp)
target_link_libraries(${MY_TARGET} PRIVATE qboot)
set_default_property(${MY_TARGET} ..)

set(MY_TARGET comp)
add_executable(${MY_TARGET} src/comp.cpp)
target_link_libraries(${MY_TARGET} PRIVATE qboot)
set_default_property(${MY_TARGET} ..)
