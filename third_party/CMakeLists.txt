#build third-party ensmallen lib
option(BUILD_ENSMALLEN_LIB "Build the third-party ensmallen library" ON)

#build third-party dlib lib
option(BUILD_DLIB_LIB "Build the third-party dlib library" ON)

# Armadillo
include(${CMAKE_CURRENT_LIST_DIR}/armadillo.cmake)

# Ensmallen
if(BUILD_ENSMALLEN_LIB)
  include(${CMAKE_CURRENT_LIST_DIR}/ensmallen.cmake)
endif()

# DLib
if(BUILD_DLIB_LIB)
  include(${CMAKE_CURRENT_LIST_DIR}/dlib.cmake)
endif()

# Boost
include(${CMAKE_CURRENT_LIST_DIR}/boost_headers.cmake)

# GoogleTest
if(PACKAGE_TESTS)
  include(${CMAKE_CURRENT_LIST_DIR}/googletest.cmake)
endif()
