# ------------------------------------------------------------------------------
# Copyright 2021-2024 Intel Corporation
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy of
# the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations under
# the License.
# ------------------------------------------------------------------------------

# Check the version of CMake.
cmake_minimum_required(VERSION 3.24)

# ##############################################################################
# Define our package information.
# ##############################################################################

set(PACKAGE_VERSION_MAJOR "0")
set(PACKAGE_VERSION_MINOR "4")
set(PACKAGE_VERSION_PATCH "4")
set(PACKAGE_NAME "HybridQuantumClassicalLibrary")
set(PACKAGE_DESCRIPTION
    "A high-level library to build Hybrid Quantum-Classical algorithms.")
set(PACKAGE_AUTHOR "Intel Corporation")
set(PACKAGE_AUTHOR_EMAIL "AQUA Group at Intel Corporation")

project(
  hybrid-quantum-classical-library
  VERSION 0.4.4
  LANGUAGES CXX)

# ##############################################################################
# Set some compiler options related to C++.
# ##############################################################################
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ##############################################################################
# Define build type
# ##############################################################################
set(DEFAULT_BUILD_TYPE "Release")

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  message(
    STATUS "CMAKE_BUILD_TYPE: None, set to ${DEFAULT_BUILD_TYPE} by default")
  set(CMAKE_BUILD_TYPE
      ${DEFAULT_BUILD_TYPE}
      CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
else()
  message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif()

set(CMAKE_C_COMPILER ${CMAKE_C_COMPILER})
set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER})

set(COMMON_GNU_DEBUG_FLAGS " -O0 -Wl,-z,relro -fstack-protector ")
set(COMMON_GNU_RELEASE_FLAGS
    " -Wall -Wformat -Wformat-security -Wparentheses -fPIC -fPIE -pie -D_FORTIFY_SOURCE=2 "
)

set(COMMON_CLANG_DEBUG_FLAGS
    " -O0 -fstack-protector -Wall -Wextra -Wformat -Wformat-security -fPIC -fPIE -D_FORTIFY_SOURCE=2 "
)
set(COMMON_CLANG_RELEASE_FLAGS
    " -O3 -DNDEBUG -Wall -Wextra -Wformat -Wformat-security -Werror=format-security -fPIC -fPIE -D_FORTIFY_SOURCE=2 -fstack-protector-strong -fwrapv -fno-delete-null-pointer-checks "
)

set(COMMON_INTEL_DEBUG_FLAGS
    " -O0 -fstack-protector -Wall -Wextra -Wformat -Wformat-security -fPIC -fPIE "
)
set(COMMON_INTEL_RELEASE_FLAGS
    " -O3 -DNDEBUG -Wall -Wremarks -Wformat -Wformat-security -Werror=format-security -fPIC -fPIE "
)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_C_FLAGS_DEBUG "${COMMON_GNU_DEBUG_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG "${COMMON_GNU_DEBUG_FLAGS}")
  set(CMAKE_C_FLAGS_RELEASE "${COMMON_GNU_RELEASE_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${COMMON_GNU_RELEASE_FLAGS}")
  add_compile_options(-fno-strict-overflow -fno-delete-null-pointer-checks
                      -fwrapv -fstack-protector-strong)
elseif(("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
       OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"))
  set(CMAKE_C_FLAGS_DEBUG "${COMMON_CLANG_DEBUG_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG "${COMMON_CLANG_DEBUG_FLAGS}")
  set(CMAKE_C_FLAGS_RELEASE "${COMMON_CLANG_RELEASE_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${COMMON_CLANG_RELEASE_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM")
  set(CMAKE_C_FLAGS_DEBUG "${COMMON_INTEL_DEBUG_FLAGS}")
  set(CMAKE_CXX_FLAGS_DEBUG "${COMMON_INTEL_DEBUG_FLAGS}")
  set(CMAKE_C_FLAGS_RELEASE "${COMMON_INTEL_RELEASE_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${COMMON_INTEL_RELEASE_FLAGS}")
  add_compile_options(-fno-delete-null-pointer-checks -fwrapv
                      -fstack-protector-strong)
endif()

set(CMAKE_CONFIGURATION_TYPES
    "Release"
    CACHE STRING "" FORCE)
set(CMAKE_CONFIGURATION_TYPES
    "Debug"
    CACHE STRING "" FORCE)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

function(package_add_test_with_libraries TESTNAME FILES LIBRARIES
         TEST_WORKING_DIRECTORY)
  add_executable(${TESTNAME} ${FILES})
  target_link_libraries(${TESTNAME} gtest gmock gtest_main ${LIBRARIES})
endfunction()

# Armadillo
set(ARMADILLO_PREFIX armadillo-12.8.0)

# Add ARMA_USE_WRAPPER
add_definitions(-DARMA_USE_WRAPPER)

# build third_party directory
add_subdirectory(third_party)

# build source directory
add_subdirectory(src)

# build tests
option(PACKAGE_TESTS "Build the tests" OFF)
if(PACKAGE_TESTS)
  enable_testing()
  include(GoogleTest)

  if(EnableCodeCoverage)
    add_compile_options(-fprofile-arcs -ftest-coverage)
  endif()
endif()

# include licensing
set(LIC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE
              ${CMAKE_CURRENT_SOURCE_DIR}/third-party-programs.txt)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE
     DESTINATION ${CMAKE_BINARY_DIR}/licensing)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/third-party-programs.txt
     DESTINATION ${CMAKE_BINARY_DIR}/licensing)
