#!/bin/bash

display_usage() {
    echo "Usage: `basename $0` -b build_directory -v lib_version"
    echo "Usage: e.g.: ./package.sh -b ./build -v v0.1.0"
}

# Capture command line parameters

cmdline="$0"
for a in "$@"; do
    cmdline="$cmdline '$a'"
done

#####################################################################
## Initialization
#####################################################################

while getopts ":b:v:" opt; do
    case "$opt" in
        b)
            BUILD_DIR_PREFIX=$OPTARG ;;
        v)
            LIB_RELEASE_VERSION=$OPTARG ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            display_usage
            exit 1
    esac
done

# if less than two arguments supplied, display usage
if [  $# -le 2 ]
then
  display_usage
  exit 1
fi

ABS_BUILD_DIR_PREFIX=$(realpath ${BUILD_DIR_PREFIX})
cd ${ABS_BUILD_DIR_PREFIX}
cmake ..
make -j32

echo "Building Hybrid Quantum-Classical Library package ${LIB_RELEASE_VERSION}.tar.gz"
tar czf ${LIB_RELEASE_VERSION}.tar.gz lib include share licensing
echo "Package build complete"

