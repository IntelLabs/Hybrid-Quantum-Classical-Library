set(ARMADILLO_LIBNAME armadillo)
set(ARMADILLO_INCNAME armadillo)

set(ARMADILLO_URL ${CMAKE_CURRENT_SOURCE_DIR}/${ARMADILLO_PREFIX}.tar.xz)
set(ARMADILLO_EXTRACT_DIR ${CMAKE_BINARY_DIR}/third_party)
set(ARMADILLO_SRC_DIR ${ARMADILLO_EXTRACT_DIR}/${ARMADILLO_PREFIX})
set(ARMADILLO_INCLUDE_DIR ${ARMADILLO_EXTRACT_DIR}/${ARMADILLO_PREFIX}/include)

execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${ARMADILLO_URL}
                 WORKING_DIRECTORY ${ARMADILLO_EXTRACT_DIR})

include(ExternalProject)
ExternalProject_Add(armadillo
  URL               ${ARMADILLO_URL}
  SOURCE_DIR        ${ARMADILLO_SRC_DIR}
  BINARY_DIR        "${ARMADILLO_SRC_DIR}"
  CMAKE_CACHE_ARGS  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}
  CMAKE_COMMAND     ${CMAKE_COMMAND}
)

