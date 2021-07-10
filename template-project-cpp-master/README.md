# Description
This is a template for a C++ project with CMake, gtest/gmock and a .gitlab.yml.  
There is a data directory that is taken care of by CMake.

# Installation
## Prerequisites
C++ version must be at least C++11.  
CMake must be at least version 3.1.

## Compilation
```
mkdir build
cd build
cmake ..
cmake --build .
```
There are 2 output executable, one for the test and one for the project.