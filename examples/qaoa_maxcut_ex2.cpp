//===----------------------------------------------------------------------===//
// INTEL CONFIDENTIAL
//
// Copyright  2023 Intel Corporation.
//
// This software and the related documents are Intel copyrighted materials, and
// your use of them is governed by the express license under which they were
// provided to you ("License"). Unless the License provides otherwise, you may
// not use, modify, copy, publish, distribute, disclose or transmit this
// software or the related documents without Intel's prior written permission.
//
// This software and the related documents are provided as is, with no express
// or implied warranties, other than those that are expressly stated in the
// License.
//===----------------------------------------------------------------------===//
//
// QAOA
// Quantum Approximate Optimization Algorithm

//
// 5-qubit Maxcut example
// Reference: https://pennylane.ai/qml/demos/tutorial_qaoa_maxcut.html
//
// For a graph, a maximum cut is a cut whose size is at least the size of any
// other cut. That is, it is a partition of the graph's vertices into two
// complementary sets S and T, such that the number of edges between S and T is
// as large as possible. Finding such a cut is known as the max-cut problem -
// Wikipedia
//
//   _______
//  /       \    A
// /      1__\________2
// \     ╱|   \       |
//  \   / |    \      |
//   \ /  |     \     |
//    /\  |      \    |
// 0 /  \ |       \   |
//   \   \|        \  |
//    \   |\        \ |
//     \  | \    B   \|
//      \ |  \        |\
//       \|___\_______| \
//       3     \      4
//
// Graph = [(0, 1), (0, 3), (1, 2), (1, 3), (2, 4), (3, 4)]
//
// 0:
// ──RX(0.10)─╭●──RZ(0.10)─╭●──RZ(0.10)─────────────────────────────────────────────────────────┤
// <Z> 1:
// ──RX(0.10)─╰────────────┤────────────╭●──RZ(0.10)──╭●──RZ(0.10)──────────────────────────────┤
// <Z> 2:
// ──RX(0.10)──────────────┤────────────╰─────────────┤─────────────╭●──RZ(0.10)────────────────┤
// <Z> 3:
// ──RX(0.10)──────────────╰──────────────────────────╰─────────────┤─────────────╭●──RZ(0.10)──┤
// <Z> 4:
// ──RX(0.10)───────────────────────────────────────────────────────╰─────────────╰─────────────┤
// <Z>
//
//
// Theoretical value C = 5; The max cut is shown in the above graph which
// divides the graph into A & B group resulting into 5 edges being cut through.
//
// Expected Output
// $ export HQ=<path to project directory>
// $ ./intel-quantum-compiler -I $HQ/build/include -L $HQ/build/lib -larmadillo
// -lhqcl $HQ/examples/qaoa_maxcut_ex2.cpp
// $ ./qaoa_maxcut_ex2
//
// Original Hamiltonian:
// -3.000000 [ ]
// 0.500000 [ Z0 Z1 ]
// 0.500000 [ Z0 Z3 ]
// 0.500000 [ Z1 Z2 ]
// 0.500000 [ Z1 Z3 ]
// 0.500000 [ Z2 Z4 ]
// 0.500000 [ Z3 Z4 ]

// Quantum approximate optimization algorithm (QAOA) for the MaxCut Problem
// Optimized Parameters using Simulated Annealing (SA):
//     8.6377
//    19.9916
//    18.8962
//    29.7367
//     0.5858
//     1.0941
//    18.1669
//   -31.3474
//     6.1454
//    -2.4512
//   -18.7885
//    10.6296
//    14.7831
//    17.6196
//     7.3579
//    11.7967
//   -30.7771
//    14.4257
//    55.8855
//     0.7003
// Simulated Annealing (SA) execution count: 166369
// Maximum Cut for the given graph is C = 3.99992
//
// Note: The actual output is not close to the expected output which is okay. Perhaps the 
// the optimization algorithm used (Simulated Annealing in this case) does not get to the
// expected output. Trying out different optimization algorithms may help get to the 
// expected output. Also, adding more number of layers would help but might take longer time
// to execute.

/// Production mode
#include <clang/Quantum/quintrinsics.h>

/// Development mode
// #include "../../clang/include/clang/Quantum/quintrinsics.h"

/// Quantum Runtime Library APIs
#include <quantum.hpp>

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <armadillo>
#include <ensmallen.hpp>

#include "SymbolicOperatorUtils.hpp"

using namespace arma;
using namespace hybrid::quantum::core;
using namespace iqsdk;

// Initial setup
const int N = 5;
const int num_layers = 10;   // Number of layering circuits
qbit QubitReg[N];
cbit CReg[N];

// Special global array to hold dynamic parameters for quantum algorithm
double QuantumVariableParams[2 * num_layers];

// Optimization steps
int steps_count = 0;

///
///@brief  defines the quantum kernel
///
///@return quantum_kernel - quantum kernel object
///
quantum_kernel void qaoaQ5() {
  // Index to loop over later
  int Index = 0;

  // Initialization of the qubits
  for (Index = 0; Index < N; Index++) {
    PrepZ(QubitReg[Index]);
  }

  // Apply Hadamard to get into superposition state
  for (Index = 0; Index < N; Index++) {
    H(QubitReg[Index]);
  }

  int step_diff = 2;
  for (int layer = 0; layer < num_layers; layer++) {

    // Part of original Hamiltonian
    // Graph = [(0, 1), (0, 3), (1, 2), (1, 3), (2, 4), (3, 4)]
    // Z0Z1
    CNOT(QubitReg[0], QubitReg[1]);
    RZ(QubitReg[1], QuantumVariableParams[step_diff * layer]);
    CNOT(QubitReg[0], QubitReg[1]);

    // Z0Z3
    CNOT(QubitReg[0], QubitReg[3]);
    RZ(QubitReg[3], QuantumVariableParams[step_diff * layer]);
    CNOT(QubitReg[0], QubitReg[3]);

    // Z1Z2
    CNOT(QubitReg[1], QubitReg[2]);
    RZ(QubitReg[2], QuantumVariableParams[step_diff * layer]);
    CNOT(QubitReg[1], QubitReg[2]);

    // Z1Z3
    CNOT(QubitReg[1], QubitReg[3]);
    RZ(QubitReg[3], QuantumVariableParams[step_diff * layer]);
    CNOT(QubitReg[1], QubitReg[3]);

    // Z2Z4
    CNOT(QubitReg[2], QubitReg[4]);
    RZ(QubitReg[4], QuantumVariableParams[step_diff * layer]);
    CNOT(QubitReg[2], QubitReg[4]);

    // Z3Z4
    CNOT(QubitReg[3], QubitReg[4]);
    RZ(QubitReg[4], QuantumVariableParams[step_diff * layer]);
    CNOT(QubitReg[3], QubitReg[4]);

    // Part of mixer Hamiltonian w/ X rotations on all qubits
    for (Index = 0; Index < N; Index++) {
      RX(QubitReg[Index], QuantumVariableParams[step_diff * layer + 1]);
    }
  }
}

///
///@brief runs the quantum kernel over number of layering circuits
///
///@param iqs_device - IQS device object
///@param params - Optimization Parameters
///@param symbop - SymbolicOperator object
///@param num_layers - Number of layering circuits
///@return double - expectation value
///
double run_qkernel(FullStateSimulator &iqs_device, const arma::mat &params,
                   SymbolicOperator &symbop, const int &num_layers) {
  double exp_val = 0.0;

  int step_diff = 2;
  for (int layer = 0; layer < num_layers; layer++) {
    QuantumVariableParams[step_diff * layer] = params[step_diff * layer];
    QuantumVariableParams[step_diff * layer + 1] =
        params[step_diff * layer + 1];
  }

  for (const auto &pstr : symbop.getOrderedPStringList()) {
    std::vector<double> probs;

    std::vector<std::reference_wrapper<qbit>> qids;
    for (int qubit = 0; qubit < N; ++qubit) {
      qids.push_back(std::ref(QubitReg[qubit]));
    }

    qaoaQ5();

    probs = iqs_device.getProbabilities(qids);

    double current_pstr_val =
        symbop.op_sum[pstr].real() *
        SymbolicOperatorUtils::getExpectValSglPauli(pstr, probs, N);

    exp_val += current_pstr_val;
  }

  return exp_val;
}

///
///@brief Object class to apply optimization step
///
///
class EnergyOfAnsatz {
public:
  EnergyOfAnsatz(SymbolicOperator &_symbop, arma::mat &_params,
                 FullStateSimulator &_iqs_device, const int &_num_layers)
      : symbop(_symbop), params(_params), iqs_device(_iqs_device),
        num_layers(_num_layers) {}

  double Evaluate(const arma::mat &theta) {
    steps_count++;

    return run_qkernel(iqs_device, params, symbop, num_layers);
  }

private:
  const arma::mat &params;
  SymbolicOperator symbop;
  FullStateSimulator &iqs_device;
  const int num_layers;
};

int main() {

  // Original Hamiltonian
  SymbolicOperator H_symbop;
  pstring inp_z1{{0, 'Z'}, {1, 'Z'}};
  H_symbop.addTerm(inp_z1, 0.5);

  pstring inp_z2{{0, 'Z'}, {3, 'Z'}};
  H_symbop.addTerm(inp_z2, 0.5);

  pstring inp_z3{{1, 'Z'}, {2, 'Z'}};
  H_symbop.addTerm(inp_z3, 0.5);

  pstring inp_z4{{1, 'Z'}, {3, 'Z'}};
  H_symbop.addTerm(inp_z4, 0.5);

  pstring inp_z5{{2, 'Z'}, {4, 'Z'}};
  H_symbop.addTerm(inp_z5, 0.5);

  pstring inp_z6{{3, 'Z'}, {4, 'Z'}};
  H_symbop.addTerm(inp_z6, 0.5);

  H_symbop.addIdentTerm(-3);

  // H_symbop = 3 * I - 0.5 * Z0Z1 - 0.5 * Z0Z3 - 0.5 * Z1Z2 - 0.5 * Z1Z3 - 0.5
  // * Z2Z4 - 0.5 * Z3Z4
  std::string H_symbop_charstring = H_symbop.getCharString();
  std::cout << "Original Hamiltonian:\n" << H_symbop_charstring << "\n";

  // Setup quantum device
  IqsConfig iqs_config(/*num_qubits*/ N,
                       /*simulation_type*/ "noiseless");
  FullStateSimulator iqs_device(iqs_config);
  if (QRT_ERROR_SUCCESS != iqs_device.ready()) {
    return -1;
  }

  // resetting optimization steps count
  steps_count = 0;

  // Using Simulated Annealing (SA) algorithm
  // Starting parameters
  arma::mat params(std::vector<double>(2 * num_layers, 0.01));

  // arbitrary function
  EnergyOfAnsatz eoa_sa_orig_hmtn(H_symbop, params, iqs_device, num_layers);

  // initialize the optimization algorithm
  ens::SA<> opt_sa_orig_hmtn(
      ens::ExponentialSchedule(), /* coolingSchedule - Instantiated cooling
                                    schedule (default ExponentialSchedule). */
      1000000, /* maxIterations - Maximum number of iterations allowed (0
                  indicates no limit). */
      1000.,   /* initT - Initial temperature. */
      1000,    /* initMoves - Number of initial iterations without changing
                  temperature. */
      100,     /* moveCtrlSweep - Sweeps per feedback move control. */
      1e-10,   /* tolerance - Tolerance to consider system frozen. */
      3,       /* maxToleranceSweep - Maximum sweeps below tolerance to consider
                  system frozen. */
      1.5,     /* maxMoveCoef - Maximum move size. */
      0.5,     /* initMoveCoef - Initial move size. */
      0.3      /* gain - Proportional control in feedback move control. */
  );

  // optimize the parameters
  double maximum_cut = -opt_sa_orig_hmtn.Optimize(eoa_sa_orig_hmtn, params);
  int sa_steps_count_orig_hmtn = steps_count;

  std::cout << "Quantum approximate optimization algorithm (QAOA) for the "
               "MaxCut Problem"
            << "\n";
  std::cout << "Optimized Parameters using Simulated Annealing (SA): "
            << "\n";
  params.print();
  std::cout << "Simulated Annealing (SA) execution count: "
            << sa_steps_count_orig_hmtn << "\n";
  std::cout << "Maximum Cut for the given graph is C = " << maximum_cut << "\n";

  return 0;
}
