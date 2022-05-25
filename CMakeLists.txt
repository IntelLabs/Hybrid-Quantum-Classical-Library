cmake_minimum_required(VERSION 3.1)

project(my_project)

# C++ version required
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_C_COMPILER ${CMAKE_C_COMPILER})
set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wformat -Wformat-security -Wparentheses")

# set installation paths
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/build/_inst)

## Set names
#set(EXECUTABLE_NAME hqc)

# Add LLVM include directories
set(CMAKE_CXX_FLAGS="-I$HOME/code/intel-aqcc/_inst/include")

# Add the source directory
add_subdirectory(src)

# Add the third_party
add_subdirectory(third_party)

# Build tests
#option(PACKAGE_TESTS "Build the tests" ON)
#if(PACKAGE_TESTS)
    enable_testing()
    include(GoogleTest)
    add_subdirectory(utests)
#endif()

## Create the main executable
#add_executable(
#    ${EXECUTABLE_NAME}
#    ${Complete_Srcs}
#)
#
## Create the test executable
#add_executable(
#    ${EXECUTABLE_NAME}_test
#    ${Tests_Srcs}
#)
#target_link_libraries(${EXECUTABLE_NAME}_test PUBLIC gtest)
#target_include_directories(${EXECUTABLE_NAME} PUBLIC $HOME/code/intel-aqcc/_inst/include)

## Resources
#file(COPY "${CMAKE_SOURCE_DIR}/data" DESTINATION .)
