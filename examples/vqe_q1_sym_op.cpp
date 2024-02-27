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
// VQE
// Variational Quantum Eigensolver

//
// Expected Output
// $ export HQ=<path to project directory>
// $ ./intel-quantum-compiler -I $HQ/build/include -L $HQ/build/lib -larmadillo -lhqcl $HQ/examples/vqe_q1_sym_op.cpp
// $ ./vqe_q1_sym_op
// Hamiltonian:
// 0.250000 [ X0 ]
// 0.500000 [ Y0 ]

// Number of Qubitwise Commutation (QWC) Groups : 2
// Qubitwise Commutation (QWC) Groups:
// Group 0
// 	[ X0 ]
// Group 1
// 	[ Y0 ]

// Optimized Parameters using Simulated Annealing (SA):
//    7.0454e+01
//    2.0028e+02
// Simulated Annealing (SA) execution count: 32298
// Total energy using Simulated Annealing (SA): -0.559038
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

using namespace hybrid::quantum::core;
using namespace iqsdk;

// Initial setup
const int N = 1;
qbit QubitReg[N];
cbit CReg[N];

// Special global array to hold dynamic parameters for quantum algorithm
double QuantumVariableParams[2 + N * 2];

// Optimization steps
static int steps_count = 0;

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

  // not part of ansatz
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
                   SymbolicOperator &symbop,
                   std::map<int, std::set<pstring>> &m_qwc_groups) {
  int basis_change_variable_param_start_indx = 2;
  double total_energy = 0.0;

  QuantumVariableParams[0] = 2 * params[0];
  QuantumVariableParams[1] = 2 * params[1];

  for (auto &m_qwc_group : m_qwc_groups) {
    std::vector<double> ProbReg;

    std::vector<double> variable_params;
    variable_params.reserve(N * 2);

    SymbolicOperatorUtils::applyBasisChange((*m_qwc_group.second.begin()),
                                            variable_params, N);

    std::vector<std::reference_wrapper<qbit>> qids;
    for (int qubit = 0; qubit < N; ++qubit) {
      qids.push_back(std::ref(QubitReg[qubit]));
    }

    QuantumVariableParams[2] = 0;
    QuantumVariableParams[3] = 0;

    for (auto indx = 0; indx < variable_params.size(); ++indx) {
      QuantumVariableParams[basis_change_variable_param_start_indx + indx] =
          variable_params[indx];
    }

    vqeQ1();

    ProbReg = iqs_device.getProbabilities(qids);

    double current_pstr_val = SymbolicOperatorUtils::getExpectValSetOfPaulis(
        symbop, m_qwc_group.second, ProbReg, N);
    total_energy += current_pstr_val;
  }

  return total_energy;
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

    double total_energy = 0.0;

    return run_qkernel(iqs_device, params, symbop, m_qwc_groups);
  }

private:
  const arma::mat &params;
  SymbolicOperator &symbop;
  QWCMap &m_qwc_groups;
  FullStateSimulator &iqs_device;
};

int main() {

  SymbolicOperator so;
  pstring inp_y1{{0, 'Y'}};
  so.addTerm(inp_y1, 0.5);
  pstring inp_y2{{0, 'X'}};
  so.addTerm(inp_y2, 0.25);

  std::string charstring = so.getCharString();
  std::cout << "Hamiltonian:\n" << charstring << "\n";

  // QWC
  QWCMap m_qwc_groups =
      SymbolicOperatorUtils::getQubitwiseCommutationGroups(so);
  std::cout << "Number of Qubitwise Commutation (QWC) Groups : "
            << m_qwc_groups.size() << std::endl;
  std::cout << "Qubitwise Commutation (QWC) Groups: " << m_qwc_groups
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

  // Using SA algorithm
  // Starting parameters
  arma::mat params_sa(std::vector<double>{0.1, 0.1});

  QuantumVariableParams[0] = params_sa[0];
  QuantumVariableParams[1] = params_sa[1];

  // arbitrary function
  EnergyOfAnsatz eoa_sa(so, m_qwc_groups, params_sa, iqs_device);

  // initialize the optimization algorithm
  ens::SA<> opt_sa(
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
  double opt_sa_val = opt_sa.Optimize(eoa_sa, params_sa);
  int sa_steps_count = steps_count;

  std::cout << "Optimized Parameters using Simulated Annealing (SA): "
            << "\n";
  params_sa.print();
  std::cout << "Simulated Annealing (SA) execution count: " << sa_steps_count
            << "\n";
  std::cout << "Total energy using Simulated Annealing (SA): " << opt_sa_val
            << "\n";

  return 0;
}
