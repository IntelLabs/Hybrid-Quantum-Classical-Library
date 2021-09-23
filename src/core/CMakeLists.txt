# Set sources
include(GoogleTest)
include(${CMAKE_SOURCE_DIR}/macros.cmake)
set(coredir core)

set(Root_Srcs
    #${Srcdir}/example.hpp
    #${Srcdir}/example.cpp
    ${coredir}/SymbolicOperator.hpp
    ${coredir}/SymbolicOperator.cpp
    ${coredir}/SymbolicOperatorUtils.hpp
    ${coredir}/SymbolicOperatorUtils.cpp
)

set(Complete_Srcs
    ${Root_Srcs}
    #${Srcdir}/main.cpp

    PARENT_SCOPE
)

# Test sources
set(Tests_Srcs
    ${Root_Srcs}
    #${Srcdir}/example_test.cpp
    ${coredir}/SymbolicOperator_Test.cpp
    ${coredir}/SymbolicOperatorUtils_Test.cpp

    PARENT_SCOPE
)

file(GLOB_RECURSE SOURCES LIST_DIRECTORIES true *.hpp *.cpp)

set(BINARY SymbolicOperator)
include_directories(${CMAKE_BINARY_DIR}/third_party/googletest-release-1.11.0/googletest/include ${CMAKE_BINARY_DIR}/third_party/${ARMADILLO_PREFIX}/include)
add_library(${BINARY}_lib STATIC ${SOURCES})

enable_testing()

#package_add_test_with_libraries(example_test example_test.cpp ${BINARY}_lib "core")
package_add_test_with_libraries(SymbolicOperator_Test SymbolicOperator_Test.cpp ${BINARY}_lib "core")
package_add_test_with_libraries(SymbolicOperatorUtils_Test SymbolicOperatorUtils_Test.cpp ${BINARY}_lib "core")

