Hybrid Quantum-Classical Library (HQCL)

A high-level library to automate the execution of Hybrid Quantum-Classical algorithms with the Intel(R) Quantum SDK. The library allows the user a high degree of programmability in defining and manipulating hybrid algorithms during runtime, while maintaining the faster program execution of the fully compiled (as opposed to transpiled) Intel quantum stack.

Contact: [Atul Kulkarni](atul.kulkarni@intel.com), [Nicolas Sawaya](nicolas.sawaya@intel.com), [Shavi Premaratne](shavindra.premaratne@intel.com)

## **Quantum Stack Build Instructions**

In order to compile and run any of the hybrid algorithms, please follow the steps described below.

Note that the until the Tall Pine Lake repo’s release build is ready (See TPL2-468 and TPL2-459) we’ll be running in ‘developer mode’.


**Step#1: Compiler build and set-up**

- Clone the intel-aqcc repo in the development environment

  ```
  $ git clone https://github.com/intel-sandbox/applications.quantum.compiler-llvm10.git
  $ git checkout release
  ```

- Clone the config-json repo in the development environment
  ```
  $ cd applications.quantum.compiler-llvm10
  $ git clone https://github.com/intel-sandbox/applications.quantum.config-json.git extern/config-json
  ```

- Build

  ```
  $ mkdir build
  $ cd build
  $ cmake -DCMAKE_INSTALL_PREFIX=<path to install dir> -DLLVM_ENABLE_PROJECTS="clang;lld"   -DLLVM_EXPERIMENTAL_TARGETS_TO_BUILD=Quantum "Unix Makefiles" ../llvm
  $ make -j30
  ```

**Step#2: QRT library build**

- Clone the TPL repo

  `$ git clone https://github.com/intel-sandbox/applications.quantum.tall-pine-lake.git`

- Build (Clean-up build directory if already existing so that all the build flags take effect)

  ```
  $ mkdir build
  $ cd build
  $ CXX=g++ cmake -DEnableIntelQS=ON -DEnableTPLTrace=ON -DEnableUtest=OFF -DEnableMKL=OFF  -DEnableBlas=OFF ..
  $ make -j30
  ```

  Note the absolute path of your TPL directory. Refer to it as  <base-dir>/tall-pine-lake in all the commands below.

**Step#3: Copy configuration files and set up library path (Running from compiler’s build directory)**

- Copy ‘Configurations’ directory to compiler’s build directory

  `$ cp -r <base-dir>/tall-pine-lake/Configurations .`

- Copy ‘Calibrations’ directory to compiler’s build directory

  `$ cp -r <base-dir>/tall-pine-lake/Calibrations .`

- Environment variables for IQS

  `$ export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:<base-dir>/tall-pine-lake/build/iqs-build/lib/"`

**Step#4: Compile a hybrid algorithm (Running from compiler’s build directory)**

  ```
  $ ../intel-aqcc-ui.sh -q -I<base-dir>/tall-pine-lake/QRT_Model/include -L<base-dir>/  tall-pine-lake/build/lib/ -L<base-dir>/tall-pine-lake/build/iqs-build/lib/ ../quantum-examples/ hybrid_examples/tfd_q4.cpp
  ```

**Step#5: Run the hybrid algorithm from above step (Running from compiler’s build directory)**

  ```
  $ ./tfd_q4`
  ```

Repeat steps#4 and #5 after making any changes to the algorithm source.

## **Hybrid Quantum-Classical Library Build Instructions**

- Compile the library using the following instructions

  ```
  $ git clone https://github.com/intel-sandbox/ applications.quantum.hybrid-quantum
  $ cd hybrid-quantum
  $ mkdir build
  $ cd build
  $ cmake ..
  $ make -j32
  ```

- Once the library is built, compile a hybrid algorithm from the hybrid-quantum examples directory. The examples needs to be compiled from compiler’s build directory.

  ```
  $ export TPL=<path to tallpine project>/tall-pine-lake
  $ export OMP_NUM_THREADS=1 # For running w/ IQS
  $ LD_LIBRARY_PATH=<path to HQCL project>/build/ lib:$TPL/build/iqs-build/lib:$LD_LIBRARY_PATH
  $ ../intel-aqcc-ui.sh -q -I$TPL/QRT_Model/include   -L$TPL/build/lib/ -L$TPL/build/iqs-build/lib/ -I <path  to HQCL project>/build/include -L<path to HQCL   project>/build/lib -larmadillo -lhqcl -lpthread   -lgtest ../../hybrid-quantum/examples/  vqe_q2_symbolic_operator.cpp
  ```

