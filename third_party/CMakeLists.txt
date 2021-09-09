include(${CMAKE_CURRENT_LIST_DIR}/googletest.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/armadillo.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ensmallen.cmake)

## Where gtest's .h files can be found.
#include_directories(
#  googletest/include
#  googletest
#)
#
#add_definitions(-DGTEST_LANG_CXX11=1)
#add_definitions(-DGTEST_HAS_TR1_TUPLE=0)
#
#list(APPEND LIBS pthread)
#
## Regrouping all the third_party sources
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/googletest/googletest/include)
#
#set (SOURCES
# googletest/googletest/src/gtest-all.cc
# googletest/googlemock/src/gmock-all.cc
#)
#
#add_library(gtest ${SOURCES})
#
##add_subdirectory(googletest)

