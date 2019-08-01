set(MY_TARGET algebra)

file(GLOB SOURCES "src/*.cpp" "src/*.c")
file(GLOB HEADERS "include/*.hpp" "include/*.h")

add_library(${MY_TARGET} STATIC ${SOURCES} ${HEADERS})

target_link_libraries(${MY_TARGET} INTERFACE ${MPFR_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARY})
target_include_directories(${MY_TARGET} PUBLIC
	"${MPFR_INCLUDE_DIR}"
	"${GMP_INCLUDE_DIR}"
	include
)

set_target_properties(
	${MY_TARGET} PROPERTIES
	OUTPUT_NAME_DEBUG ${MY_TARGET}-debug
	OUTPUT_NAME_RELEASE ${MY_TARGET}
)
