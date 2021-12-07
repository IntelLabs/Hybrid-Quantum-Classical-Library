cmake_minimum_required(VERSION 3.1)

project(hybrid-quantum)

# C++ version required
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_C_COMPILER ${CMAKE_C_COMPILER})
set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wformat -Wformat-security -Wparentheses")

# set installation paths
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/build/_inst)

# Armadillo
set(ARMADILLO_PREFIX armadillo-10.6.2)

# Add LLVM include directories
set(CMAKE_CXX_FLAGS="-I$HOME/code/intel-aqcc/_inst/include")

# build source directory
add_subdirectory(src)

# build third_party directory
add_subdirectory(third_party)

# build tests
option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
    enable_testing()
    include(GoogleTest)
    add_subdirectory(utests)
endif()

