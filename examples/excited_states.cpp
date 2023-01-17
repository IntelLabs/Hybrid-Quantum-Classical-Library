//===----------------------------------------------------------------------===//
// INTEL CONFIDENTIAL
//
// Copyright 2022-2023 Intel Corporation.
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
// VQE
// Variational Quantum Eigensolver
// Excited states - Folded Spectrum Method

//
// Expected Output
// $ export HQ=<path to project directory>
// $ ./intel-quantum-compiler -I $HQ/build/include -L $HQ/build/lib -larmadillo -lhqcl $HQ/examples/excited_states.cpp
// $ ./excited_states
// Original Hamiltonian:
// 1.500000 [ ]
// 1.000000 [ X0 ]
// -0.500000 [ Z0 ]
//
// Number of Qubitwise Commutation (QWC) Groups for original hamiltonian : 2
// Qubitwise Commutation (QWC) Groups:
// Group 0
// 	[ ]
// 	[ X0 ]
// Group 1
// 	[ Z0 ]
//
// Folded Hamiltonian:
// 1.610000 [ ]
// -1.200000 [ X0 ]
// 0.600000 [ Z0 ]
//
// Number of Qubitwise Commutation (QWC) Groups for folded hamiltonian: 2
// Qubitwise Commutation (QWC) Groups:
// Group 0
// 	[ ]
// 	[ X0 ]
// Group 1
// 	[ Z0 ]
//
// Result for Original Hamiltonian:
// Optimized Parameters using Simulated Annealing (SA):
//    2.6150e+01
//    2.2777e+02
// Simulated Annealing (SA) execution count: 30422
// Total energy using Simulated Annealing (SA): 0.381987
// Result for Folded Hamiltonian:
// Optimized Parameters using Simulated Annealing (SA):
//   -1.0423e+02
//    7.6969e+01
// Simulated Annealing (SA) execution count: 30696
// Total energy using Simulated Annealing (SA): 0.268385
//
// First excited state with respect to the original hamiltonian: 2.61801

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
double QuantumVariableParams[2 + N * 2];

// Optimization steps
int steps_count = 0;

///
///@brief  defines the quantum kernel
///
///@return quantum_kernel - quantum kernel object
///
quantum_kernel void vqeQ1() {

  // Index to loop over later
  int Index = 0;

  // Initialization of the qubits
  for (Index = 0; Index < N; Index++) {
    PrepZ(QubitReg[Index]);
  }

  // Ansatz
  RY(QubitReg[0], QuantumVariableParams[0]);
  RX(QubitReg[0], QuantumVariableParams[1]);

  // For changing Basis; not part of ansatz
  RY(QubitReg[0], QuantumVariableParams[2]);
  RX(QubitReg[0], QuantumVariableParams[3]);
}

///
///@brief executes the quantum kernel
///
///@param iqs_device - IQS device object
///@param params - Optimization Parameters
///@param symbop - SymbolicOperator object
///@param m_qwc_groups - Qubit-wise commutation groups
///@return double - expectation value
///
double run_qkernel(FullStateSimulator &iqs_device, const arma::mat &params,
                   SymbolicOperator &symbop, const QWCMap &m_qwc_groups) {
  int basis_change_variable_param_start_indx = 2;
  double exp_val = 0.0;

  QuantumVariableParams[0] = 2 * params[0];
  QuantumVariableParams[1] = 2 * params[1];

  for (auto &m_qwc_group : m_qwc_groups) {
    std::vector<double> probs;

    std::vector<double> variable_params(N * 2);

    QuantumVariableParams[2] = 0;
    QuantumVariableParams[3] = 0;

    SymbolicOperatorUtils::applyBasisChange(m_qwc_group.second, variable_params,
                                            N);

    for (auto indx = 0; indx < variable_params.size(); ++indx) {
      QuantumVariableParams[basis_change_variable_param_start_indx + indx] =
          variable_params[indx];
    }

    std::vector<std::reference_wrapper<qbit>> qids;
    for (int qubit = 0; qubit < N; ++qubit) {
      qids.push_back(std::ref(QubitReg[qubit]));
    }

    vqeQ1();

    probs = iqs_device.getProbabilities(qids);

    double current_pstr_val = SymbolicOperatorUtils::getExpectValSetOfPaulis(
        symbop, m_qwc_group.second, probs, N);

    exp_val += current_pstr_val;
  }

  return exp_val;
}

///
///@brief Object class to apply optimization step
///
class EnergyOfAnsatz {
public:
  EnergyOfAnsatz(SymbolicOperator &_symbop, QWCMap &_m_qwc_groups,
                 arma::mat &_params, FullStateSimulator &_iqs_device)
      : symbop(_symbop), m_qwc_groups(_m_qwc_groups), params(_params),
        iqs_device(_iqs_device) {}

  double Evaluate(const arma::mat &theta) {
    steps_count++;

    return run_qkernel(iqs_device, params, symbop, m_qwc_groups);
  }

private:
  const arma::mat &params;
  SymbolicOperator symbop;
  QWCMap &m_qwc_groups;
  FullStateSimulator &iqs_device;
};

int main() {

  SymbolicOperator H_symbop;
  pstring inp_y1{{0, 'X'}};
  H_symbop.addTerm(inp_y1, 1.0);

  pstring inp_y2{{0, 'Z'}};
  H_symbop.addTerm(inp_y2, -0.5);

  H_symbop.addIdentTerm(1.5);

  // H_symbop = X0 - 0.5 * Z0 + 1.5 * I0;
  std::string charstring = H_symbop.getCharString();
  std::cout << "Original Hamiltonian:\n" << charstring << "\n";

  // QWC
  QWCMap m_qwc_groups_orig_ham =
      SymbolicOperatorUtils::getQubitwiseCommutationGroups(H_symbop, N);
  std::cout << "Number of Qubitwise Commutation (QWC) Groups for original "
               "hamiltonian : "
            << m_qwc_groups_orig_ham.size() << std::endl;
  std::cout << "Qubitwise Commutation (QWC) Groups: " << m_qwc_groups_orig_ham
            << std::endl;

  // Setup quantum device
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
  arma::mat params_sa_orig_hmtn(std::vector<double>{0.1, 0.1});

  QuantumVariableParams[0] = params_sa_orig_hmtn[0];
  QuantumVariableParams[1] = params_sa_orig_hmtn[1];

  // arbitrary function
  EnergyOfAnsatz eoa_sa_orig_hmtn(H_symbop, m_qwc_groups_orig_ham,
                                  params_sa_orig_hmtn, iqs_device);

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

  // Folded Hamiltonian
  steps_count = 0; // resetting steps count

  double gamma = 2.1;
  SymbolicOperator H_fold_symbop =
      SymbolicOperatorUtils::getFoldedHamiltonian(H_symbop, gamma);
  charstring = H_fold_symbop.getCharString();
  std::cout << "Folded Hamiltonian:\n" << charstring << "\n";

  // QWC
  QWCMap m_qwc_groups_folded_ham =
      SymbolicOperatorUtils::getQubitwiseCommutationGroups(H_fold_symbop, N);
  std::cout
      << "Number of Qubitwise Commutation (QWC) Groups for folded hamiltonian: "
      << m_qwc_groups_folded_ham.size() << std::endl;
  std::cout << "Qubitwise Commutation (QWC) Groups: " << m_qwc_groups_folded_ham
            << std::endl;

  // Starting parameters
  arma::mat params_sa_folded_hmtn(std::vector<double>{0.1, 0.1});

  QuantumVariableParams[0] = params_sa_folded_hmtn[0];
  QuantumVariableParams[1] = params_sa_folded_hmtn[1];

  // arbitrary function
  EnergyOfAnsatz eoa_sa_folded_hmtn(H_fold_symbop, m_qwc_groups_folded_ham,
                                    params_sa_folded_hmtn, iqs_device);

  // initialize the optimization algorithm
  ens::SA<> opt_sa_folded_hmtn(
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
  double opt_sa_val_folded_hmtn =
      opt_sa_folded_hmtn.Optimize(eoa_sa_folded_hmtn, params_sa_folded_hmtn);
  int sa_steps_count_folded_hmtn = steps_count;

  std::cout << "Result for Original Hamiltonian: " << std::endl;
  std::cout << "Optimized Parameters using Simulated Annealing (SA): "
            << "\n";
  params_sa_orig_hmtn.print();
  std::cout << "Simulated Annealing (SA) execution count: "
            << sa_steps_count_orig_hmtn << "\n";
  std::cout << "Total energy using Simulated Annealing (SA): "
            << opt_sa_val_orig_hmtn << "\n";

  std::cout << "Result for Folded Hamiltonian: " << std::endl;
  std::cout << "Optimized Parameters using Simulated Annealing (SA): "
            << "\n";
  params_sa_folded_hmtn.print();
  std::cout << "Simulated Annealing (SA) execution count: "
            << sa_steps_count_folded_hmtn << "\n";
  std::cout << "Total energy using Simulated Annealing (SA): "
            << opt_sa_val_folded_hmtn << "\n";

  arma::mat h_fold_params(
      std::vector<double>{params_sa_folded_hmtn[0], params_sa_folded_hmtn[1]});
  std::cout
      << "\nFirst excited state with respect to the original hamiltonian: "
      << run_qkernel(iqs_device, h_fold_params, H_symbop, m_qwc_groups_orig_ham)
      << "\n";

  return 0;
}
