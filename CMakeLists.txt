cmake_minimum_required(VERSION 3.1)

project(my_project)

# C++ version required
set(CMAKE_CXX_STANDARD 14)

# Set names
set(EXECUTABLE_NAME my_executable)

# Add the source directory
add_subdirectory(src)

# Add the third_party
add_subdirectory(third_party)

# Create the main executable
add_executable(
    ${EXECUTABLE_NAME}

    ${Complete_Srcs}
)

# Create the test executable
add_executable(
    ${EXECUTABLE_NAME}_test

    ${Tests_Srcs}
)
target_link_libraries(${EXECUTABLE_NAME}_test PUBLIC gtest)


# Resources
file(COPY "${CMAKE_SOURCE_DIR}/data" DESTINATION .)