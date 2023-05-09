Hybrid Quantum-Classical Library (HQCL)

A high-level library to automate the execution of Hybrid Quantum-Classical algorithms with the Intel(R) Quantum SDK. The library allows the user a high degree of programmability in defining and manipulating hybrid algorithms during runtime, while maintaining the faster program execution of the fully compiled (as opposed to transpiled) Intel quantum stack.

Contact: [Atul Kulkarni](mailto:atul.kulkarni@intel.com), [Nicolas Sawaya](mailto:nicolas.sawaya@intel.com), [Shavi Premaratne](mailto:shavindra.premaratne@intel.com)

## **Hybrid Quantum-Classical Library Build Instructions**

- Prequisites to using this library is to first install OpenBLAS and LAPACK, along with the corresponding development/header files or Intel(R) Math Kernel Library (Intel(R) MKL) which is available as part of the Intel(R) oneAPI Suite. This is required to use the third-party library Armadillo.

- Compile the library using the following instructions

  ```
  $ git clone https://github.com/intel-sandbox/ applications.quantum.hybrid-quantum
  $ cd applications.quantum.hybrid-quantum
  $ mkdir build
  $ cd build
  $ cmake ..
  $ make -j32
  ```

- Once the library is built, compile a hybrid algorithm from the hybrid-quantum examples directory. The examples needs to be compiled using the Intel(R) Quantum compiler available as part of Intel(R) Quantum SDK package.

  ```
  $ export OMP_NUM_THREADS=1 # For running w/ Intel(R) Quantum Simulator
  $ ./intel-quantum-compiler -I<path to HQCL project>/build/include -L<path to HQCL project>/build/lib -larmadillo -lhqcl <path to examples in HQCL>/examples/excited_states_q2.cpp
  ```

- Users will have the option to use the third-party ensmallen and dlib c++ optimization libraries. The libraries are built by default but will have the option to build individually using BUILD_ENSMALLEN_LIB=ON/OFF and BUILD_DLIB_LIB=ON/OFF cmake options.
