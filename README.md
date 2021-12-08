Hybrid Quantum API 

Contact: [Atul Kulkarni](atul.kulkarni@intel.com), [Nicolas Sawaya](nicolas.sawaya@intel.com), [Shavi Premaratne](shavindra.premaratne@intel.com)

## **Quantum Stack Build Instructions (Author: Pradnya):**

In order to compile and run any of the hybrid algorithms, please follow the steps described below.

Note that the until the Tall Pine Lake repo’s release build is ready (See
TPL2-468 and 
TPL2-459) we’ll be running in ‘developer mode’.


**Step#1: Compiler build and set-up**

- Clone the intel-aqcc repo in your development environment

```
$ git clone ssh://git@gitlab.devtools.intel.com:29418/aqua/intel-aqcc.git
$ (Optional) git checkout <branch-name>
```

- Build

```
$ mkdir build
$ cd build
$ cmake -DLLVM_ENABLE_PROJECTS="clang;lld" -DLLVM_EXPERIMENTAL_TARGETS_TO_BUILD=Quantum "Unix Makefiles" ../llvm
$ make -j30
```

**Step#2: QRT library build**

- Clone the TPL repo

`$ git clone ssh://git@gitlab.devtools.intel.com:29418/aqua/tall-pine-lake.git`

- Build (Clean-up build directory if already existing so that all the build flags take effect)

```
$ mkdir build
$ cd build
$ CXX=g++ cmake -DEnableIntelQS=ON -DEnableTPLTrace=ON -DEnableUtest=OFF -DEnableMKL=OFF -DEnableBlas=OFF ..
$ make -j30
```

Note the
absolute path of your TPL directory. I will refer to it as 
<base-dir>/tall-pine-lake in all the commands below.

For instance, mine is - /srv/raid0/pkhalate/workspace/tall-pine-lake

**Step#3: Copy configuration files and set up library path (Running from compiler’s build directory)**

- Copy ‘Configurations’ directory to compiler’s build directory

`$ cp -r <base-dir>/tall-pine-lake/Configurations .`

- Copy ‘Calibrations’ directory to compiler’s build directory

`$ cp -r <base-dir>/tall-pine-lake/Calibrations .`

- Environment variables for IQS

`$ export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:<base-dir>/tall-pine-lake/build/iqs-build/lib/"`

**Step#4: Compile a hybrid algorithm (Running from compiler’s build directory)**

```
$ ../intel-aqcc-ui.sh -q -I<base-dir>/tall-pine-lake/QRT_Model/include -L<base-dir>/tall-pine-lake/build/lib/ -L<base-dir>/tall-pine-lake/build/iqs-build/lib/ ../quantum-examples/hybrid_examples/tfd_q4.cpp
```

**Step#5: Run the hybrid algorithm from above step (Running from compiler’s build directory)**

`$ ./tfd_q4`

Repeat steps#4 and #5 after making any changes to the algorithm source.

## **Hybrid-Quantum API Build Instructions (Author: Atul):**

- Compile the API using the following instructions

```
$ git clone https://gitlab.devtools.intel.com/aqua/hybrid-quantum
$ cd hybrid-quantum
$ mkdir build
$ cd build
$ cmake -DCMAKE_INSTALL_PREFIX=<path to build dir>/_inst ..
$ make -j32
$ make install
```

- Once the library is built, compile a hybrid algorithm from the hybrid-quantum examples directory. The examples needs to be compiled from compiler’s build directory.

```
export TPL=/data4/apkulkar/code/tallpine
export OMP_NUM_THREADS=1 # For running w/ IQS
LD_LIBRARY_PATH=/data4/apkulkar/code/hybrid-quantum/build/_inst/lib:$TPL/build/iqs-build/lib:$LD_LIBRARY_PATH
$ ../intel-aqcc-ui.sh -q -I$TPL/QRT_Model/include -L$TPL/build/lib/ -L$TPL/build/iqs-build/lib/ -I /data4/apkulkar/code/hybrid-quantum/build/_inst/include -L/data4/apkulkar/code/hybrid-quantum/build/_inst/lib -larmadillo -lhybrid_quantum -lpthread -lgtest ../../hybrid-quantum/examples/vqe_q2_symbolic_operator.cpp
```

