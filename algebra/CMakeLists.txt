set(MY_TARGET algebra)

file(GLOB SOURCES "src/*.cpp" "src/*.c")
file(GLOB HEADERS "include/*.hpp" "include/*.h")

add_library(${MY_TARGET} STATIC ${SOURCES} ${HEADERS})

target_link_libraries(
	${MY_TARGET} INTERFACE
	$<$<CONFIG:Debug>:${MPFR_LIBRARY_DEBUG}>
	$<$<CONFIG:Release>:${MPFR_LIBRARY_RELEASE}>
	$<$<CONFIG:Debug>:${GMP_LIBRARY_DEBUG}>
	$<$<CONFIG:Release>:${GMP_LIBRARY_RELEASE}>
	$<$<CONFIG:Debug>:${GMPXX_LIBRARY_DEBUG}>
	$<$<CONFIG:Release>:${GMPXX_LIBRARY_RELEASE}>
)
target_include_directories(${MY_TARGET} PUBLIC
	$<$<CONFIG:Debug>:${MPFR_INCLUDE_DIR_DEBUG}>
	$<$<CONFIG:Release>:${MPFR_INCLUDE_DIR_RELEASE}>
	$<$<CONFIG:Debug>:${GMP_INCLUDE_DIR_DEBUG}>
	$<$<CONFIG:Release>:${GMP_INCLUDE_DIR_RELEASE}>
	include
)

set_target_properties(
	${MY_TARGET} PROPERTIES
	OUTPUT_NAME_DEBUG ${MY_TARGET}-debug
	OUTPUT_NAME_RELEASE ${MY_TARGET}
)
