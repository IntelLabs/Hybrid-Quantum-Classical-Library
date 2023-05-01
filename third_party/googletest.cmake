set(GOOGLETEST_PREFIX googletest-1.13.0)
set(GOOGLETEST_LIBNAME googletest)
set(GOOGLETEST_INCNAME googletest)

set(GOOGLETEST_URL ${CMAKE_CURRENT_SOURCE_DIR}/${GOOGLETEST_PREFIX}.tar.gz)
set(GOOGLETEST_EXTRACT_DIR ${CMAKE_BINARY_DIR}/third_party)
set(GOOGLETEST_SRC_DIR ${GOOGLETEST_EXTRACT_DIR}/${GOOGLETEST_PREFIX})

message(STATUS "GOOGLETEST_EXTRACT_DIR: ${GOOGLETEST_EXTRACT_DIR}")

execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${GOOGLETEST_URL}
                 WORKING_DIRECTORY ${GOOGLETEST_EXTRACT_DIR})

# Keeps the cache clean
mark_as_advanced(
                    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
                    gmock_build_tests gtest_build_samples gtest_build_tests
                    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
                )

#
# set_target_properties(gtest PROPERTIES FOLDER third_party)
# set_target_properties(gtest_main PROPERTIES FOLDER third_party)
# set_target_properties(gmock PROPERTIES FOLDER third_party)
# set_target_properties(gmock_main PROPERTIES FOLDER third_party)

option(INSTALL_GTEST "Enable installation of googletest. (Projects embedding googletest may want to turn this OFF.)" OFF)
add_subdirectory(${GOOGLETEST_EXTRACT_DIR}/${GOOGLETEST_PREFIX} ${GOOGLETEST_EXTRACT_DIR}/${GOOGLETEST_PREFIX})


