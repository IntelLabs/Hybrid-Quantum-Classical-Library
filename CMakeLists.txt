cmake_minimum_required(VERSION 3.1)

project(hybrid-quantum)

# C++ version required
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_C_COMPILER ${CMAKE_C_COMPILER})
set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wformat -Wformat-security -Wparentheses")

# Armadillo
set(ARMADILLO_PREFIX armadillo-10.6.2)

function(package_add_test_with_libraries TESTNAME FILES LIBRARIES TEST_WORKING_DIRECTORY)
    add_executable(${TESTNAME} ${FILES})
    target_link_libraries(${TESTNAME} gtest gmock gtest_main ${LIBRARIES})
endfunction()

# build source directory
add_subdirectory(src)

# build third_party directory
add_subdirectory(third_party)

# build tests
option(PACKAGE_TESTS "Build the tests" OFF)
if(PACKAGE_TESTS)
    enable_testing()
    include(GoogleTest)

    if (EnableCodeCoverage)
        add_compile_options(-fprofile-arcs -ftest-coverage)
    endif()
endif()

