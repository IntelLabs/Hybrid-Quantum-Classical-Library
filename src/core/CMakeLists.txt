# Set sources
include(GoogleTest)
set(COREDIR core)
set(LIB_NAME "hqcl")

set(Root_Srcs SymbolicOperator.hpp SymbolicOperator.cpp
              SymbolicOperatorUtils.hpp SymbolicOperatorUtils.cpp)

set(Complete_Srcs ${Root_Srcs})

# Test sources
set(Tests_Srcs ${Root_Srcs} SymbolicOperator_Test.cpp
               SymbolicOperatorUtils_Test.cpp)

# inc headers
list(APPEND INC_HEADERS SymbolicOperator.hpp SymbolicOperatorUtils.hpp)

add_library(${LIB_NAME} SHARED ${Root_Srcs})
set_target_properties(
  ${LIB_NAME} PROPERTIES PUBLIC_HEADER "${INC_HEADERS}"
                         LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

target_include_directories(
  ${LIB_NAME}
  PUBLIC ${CMAKE_BINARY_DIR}/third_party/${ARMADILLO_PREFIX}/include)

file(COPY ${INC_HEADERS} DESTINATION ${CMAKE_BINARY_DIR}/include)

if(PACKAGE_TESTS)
  enable_testing()
  include_directories(
    ${CMAKE_BINARY_DIR}/third_party/googletest-release-1.11.0/googletest/include
  )
  package_add_test_with_libraries(CoreTests "${Tests_Srcs}" ${LIB_NAME}
                                  ${COREDIR})
endif()
