# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/npsawaya/sw/hybrid-quantum

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/npsawaya/sw/hybrid-quantum/build

# Include any dependencies generated for this target.
include CMakeFiles/my_executable.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/my_executable.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/my_executable.dir/flags.make

CMakeFiles/my_executable.dir/src/example.cpp.o: CMakeFiles/my_executable.dir/flags.make
CMakeFiles/my_executable.dir/src/example.cpp.o: ../src/example.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/npsawaya/sw/hybrid-quantum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/my_executable.dir/src/example.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_executable.dir/src/example.cpp.o -c /Users/npsawaya/sw/hybrid-quantum/src/example.cpp

CMakeFiles/my_executable.dir/src/example.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_executable.dir/src/example.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/npsawaya/sw/hybrid-quantum/src/example.cpp > CMakeFiles/my_executable.dir/src/example.cpp.i

CMakeFiles/my_executable.dir/src/example.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_executable.dir/src/example.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/npsawaya/sw/hybrid-quantum/src/example.cpp -o CMakeFiles/my_executable.dir/src/example.cpp.s

CMakeFiles/my_executable.dir/src/main.cpp.o: CMakeFiles/my_executable.dir/flags.make
CMakeFiles/my_executable.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/npsawaya/sw/hybrid-quantum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/my_executable.dir/src/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/my_executable.dir/src/main.cpp.o -c /Users/npsawaya/sw/hybrid-quantum/src/main.cpp

CMakeFiles/my_executable.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/my_executable.dir/src/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/npsawaya/sw/hybrid-quantum/src/main.cpp > CMakeFiles/my_executable.dir/src/main.cpp.i

CMakeFiles/my_executable.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/my_executable.dir/src/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/npsawaya/sw/hybrid-quantum/src/main.cpp -o CMakeFiles/my_executable.dir/src/main.cpp.s

# Object files for target my_executable
my_executable_OBJECTS = \
"CMakeFiles/my_executable.dir/src/example.cpp.o" \
"CMakeFiles/my_executable.dir/src/main.cpp.o"

# External object files for target my_executable
my_executable_EXTERNAL_OBJECTS =

my_executable: CMakeFiles/my_executable.dir/src/example.cpp.o
my_executable: CMakeFiles/my_executable.dir/src/main.cpp.o
my_executable: CMakeFiles/my_executable.dir/build.make
my_executable: CMakeFiles/my_executable.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/npsawaya/sw/hybrid-quantum/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable my_executable"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/my_executable.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/my_executable.dir/build: my_executable

.PHONY : CMakeFiles/my_executable.dir/build

CMakeFiles/my_executable.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/my_executable.dir/cmake_clean.cmake
.PHONY : CMakeFiles/my_executable.dir/clean

CMakeFiles/my_executable.dir/depend:
	cd /Users/npsawaya/sw/hybrid-quantum/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/npsawaya/sw/hybrid-quantum /Users/npsawaya/sw/hybrid-quantum /Users/npsawaya/sw/hybrid-quantum/build /Users/npsawaya/sw/hybrid-quantum/build /Users/npsawaya/sw/hybrid-quantum/build/CMakeFiles/my_executable.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/my_executable.dir/depend

