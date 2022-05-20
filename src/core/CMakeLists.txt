# Set sources
include(GoogleTest)
set(coredir core)
set(LIB_NAME "hqcl")

set(Root_Srcs
    SymbolicOperator.hpp
    SymbolicOperator.cpp
    SymbolicOperatorUtils.hpp
    SymbolicOperatorUtils.cpp
)

set(Complete_Srcs
    ${Root_Srcs}
)

# Test sources
set(Tests_Srcs
    ${Root_Srcs}
    SymbolicOperator_Test.cpp
    SymbolicOperatorUtils_Test.cpp
)

# inc headers
list(APPEND INC_HEADERS
    SymbolicOperator.hpp
    SymbolicOperatorUtils.hpp
)

include_directories(${CMAKE_BINARY_DIR}/third_party/${ARMADILLO_PREFIX}/include)
add_library(${LIB_NAME} SHARED ${Root_Srcs})
set_target_properties(${LIB_NAME} PROPERTIES PUBLIC_HEADER "${INC_HEADERS}")

install(TARGETS ${LIB_NAME}
    LIBRARY DESTINATION ${CMAKE_BINARY_DIR}/lib
    PUBLIC_HEADER DESTINATION ${CMAKE_BINARY_DIR}/include)

if(PACKAGE_TESTS)
  enable_testing()
  include_directories(${CMAKE_BINARY_DIR}/third_party/googletest-release-1.11.0/googletest/include)
  package_add_test_with_libraries(SymbolicOperator_Test "${Tests_Srcs}" ${LIB_NAME} "core")
endif()
