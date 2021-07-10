# Set sources

set(Srcdir src)

set(Root_Srcs
    ${Srcdir}/example.cpp
    ${Srcdir}/example.hpp
)

set(Complete_Srcs
    ${Root_Srcs}
    ${Srcdir}/main.cpp

    PARENT_SCOPE
)

# Test sources
set(Tests_Srcs
    ${Root_Srcs}
    ${Srcdir}/example_test.cpp
    ${Srcdir}/main_test.cpp

    PARENT_SCOPE
)