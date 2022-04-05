#!/bin/bash
#Run code coverage on GNU build
rm -rf build_codecoverage
mkdir build_codecoverage
cd build_codecoverage
CXX=g++ CXXFLAGS="-fprofile-arcs -ftest-coverage" cmake -DCMAKE_INSTALL_PREFIX=/data/apkulkar/code/applications.quantum.hybrid-quantum/build/_inst -DPACKAGE_TESTS=ON -DEnableCodeCoverage=ON ..
make -j 20
./src/core/SymbolicOperator_Test
cd ..
gcovr --exclude-throw-branches --exclude-unreachable-branches --config tools/gcovr.cfg > gnu_coverage.json

gcovr --add-tracefile gnu_coverage.json --html-details coverage.html > /dev/null
rm -rf coverage
mkdir coverage
mv *.html coverage
mv coverage.css coverage
rm gnu_coverage.json
