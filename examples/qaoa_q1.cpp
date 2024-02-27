//===----------------------------------------------------------------------===//
// INTEL CONFIDENTIAL
//
// Copyright 2022-2024 Intel Corporation.
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

// 1-qubit example
// 0: ──RX(0.01)────RZ(0.01)───●─┤  <Z>

//
// Theoretical optimal value: -1
// H = Z0
// H = [[ 1  0]
//      [ 0 -1]]
//
// Using python API
// es,vs = la.eig(H)
// min_eig = min(es)
// Minimum eigenvalue: (-1+0j)
//

//
// Expected Output
// $ export HQ=<path to project directory>
// $ ./intel-quantum-compiler -I $HQ/build/include -L $HQ/build/lib -larmadillo -lhqcl $HQ/examples/qaoa_q1.cpp
// $ ./qaoa_q1
//
// Original Hamiltonian:
// 1.000000 [ Z0 ]
//
// Result for Original Hamiltonian:
// Optimized Parameters using Simulated Annealing (SA):
//    80.1106
//     7.1060
// Simulated Annealing (SA) execution count: 28186
// Total energy using Simulated Annealing (SA): -1
//

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
const int N = 1;
qbit QubitReg[N];
cbit CReg[N];

// Special global array to hold dynamic parameters for quantum algorithm
double QuantumVariableParams[2];

// Optimization steps
int steps_count = 0;

///
///@brief  defines the quantum kernel
///
///@return quantum_kernel - quantum kernel object
///
quantum_kernel void qaoaQ1() {
  // Index to loop over later
  int Index = 0;

  // Initialization of the qubits
  for (Index = 0; Index < N; Index++) {
    PrepZ(QubitReg[Index]);
  }

  // Part of mixer Hamiltonian w/ X rotations on all qubits
  for (Index = 0; Index < N; Index++) {
    RX(QubitReg[Index], QuantumVariableParams[0]);
  }

  // Part of original Hamiltonian
  for (Index = 0; Index < N; Index++) {
    RZ(QubitReg[Index], QuantumVariableParams[1]);
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
                   SymbolicOperator &symbop, int &num_layers) {
  double exp_val = 0.0;

  int step_diff = 2;
  for (int layer = 0; layer < num_layers; layer++) {
    QuantumVariableParams[0] = 2 * params[step_diff * layer];
    QuantumVariableParams[1] = 2 * params[step_diff * layer + 1];

    for (const auto &pstr : symbop.getOrderedPStringList()) {
      std::vector<double> probs;

      std::vector<std::reference_wrapper<qbit>> qids;
      for (int qubit = 0; qubit < N; ++qubit) {
        qids.push_back(std::ref(QubitReg[qubit]));
      }

      qaoaQ1();

      probs = iqs_device.getProbabilities(qids);

      double current_pstr_val =
          symbop.op_sum[pstr].real() *
          SymbolicOperatorUtils::getExpectValSglPauli(pstr, probs, N);

      exp_val += current_pstr_val;
    }
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
                 FullStateSimulator &_iqs_device, int &_num_layers)
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
  int num_layers;
};

int main() {

  int num_layers = 1;

  // Original Hamiltonian
  SymbolicOperator H_symbop;
  pstring inp_y1{{0, 'Z'}};
  H_symbop.addTerm(inp_y1, 1.0);

  // H_symbop = Z0
  std::string H_symbop_charstring = H_symbop.getCharString();
  std::cout << "Original Hamiltonian:\n" << H_symbop_charstring << "\n";

  /// Setup quantum device
  IqsConfig iqs_config(/*num_qubits*/ N,
                       /*simulation_type*/ "noiseless");
  FullStateSimulator iqs_device(iqs_config);
  if (QRT_ERROR_SUCCESS != iqs_device.ready()) {
    return -1;
  }

  // resetting steps count
  steps_count = 0;

  // Using Simulated Annealing (SA) algorithm
  // Starting parameters
  arma::mat params_sa_orig_hmtn(std::vector<double>{0.01, 0.01});

  // arbitrary function
  EnergyOfAnsatz eoa_sa_orig_hmtn(H_symbop, params_sa_orig_hmtn, iqs_device,
                                  num_layers);

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
  double opt_sa_val_orig_hmtn =
      opt_sa_orig_hmtn.Optimize(eoa_sa_orig_hmtn, params_sa_orig_hmtn);
  int sa_steps_count_orig_hmtn = steps_count;

  std::cout << "Result for Original Hamiltonian: " << std::endl;
  std::cout << "Optimized Parameters using Simulated Annealing (SA): "
            << "\n";
  params_sa_orig_hmtn.print();
  std::cout << "Simulated Annealing (SA) execution count: "
            << sa_steps_count_orig_hmtn << "\n";
  std::cout << "Total energy using Simulated Annealing (SA): "
            << opt_sa_val_orig_hmtn << "\n";

  return 0;
}

