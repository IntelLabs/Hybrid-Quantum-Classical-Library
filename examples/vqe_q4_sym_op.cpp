//===----------------------------------------------------------------------===//
// INTEL CONFIDENTIAL
//
// Copyright 2022 Intel Corporation.
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

const int N = 4;
qbit QubitReg[N];
cbit CReg[N];

/* Special global array to hold dynamic parameters for quantum algorithm */
double QuantumVariableParams[16 + (N * 2)];

int steps_count = 0;

quantum_kernel void vqeQ4() {

  // Index to loop over later
  int Index = 0;

  // Initialization of the qubits
  for (Index = 0; Index < N; Index++) {
    PrepZ(QubitReg[Index]);
  }

  // Ansatz
  RX(QubitReg[0], QuantumVariableParams[0]);
  RY(QubitReg[1], QuantumVariableParams[1]);

  CNOT(QubitReg[0], QubitReg[1]);

  RX(QubitReg[0], QuantumVariableParams[2]);
  RY(QubitReg[1], QuantumVariableParams[3]);

  RY(QubitReg[0], QuantumVariableParams[4]);
  RX(QubitReg[1], QuantumVariableParams[5]);

  CNOT(QubitReg[0], QubitReg[1]);

  RY(QubitReg[0], QuantumVariableParams[6]);
  RX(QubitReg[1], QuantumVariableParams[7]);

  RX(QubitReg[2], QuantumVariableParams[8]);
  RY(QubitReg[3], QuantumVariableParams[9]);

  CNOT(QubitReg[2], QubitReg[3]);

  RX(QubitReg[2], QuantumVariableParams[10]);
  RY(QubitReg[3], QuantumVariableParams[11]);

  RY(QubitReg[2], QuantumVariableParams[12]);
  RX(QubitReg[3], QuantumVariableParams[13]);

  CNOT(QubitReg[2], QubitReg[3]);

  RY(QubitReg[2], QuantumVariableParams[14]);
  RX(QubitReg[3], QuantumVariableParams[15]);

  // X, Y -> Z basis change - not part of ansatz
  RY(QubitReg[0], QuantumVariableParams[16]);
  RX(QubitReg[0], QuantumVariableParams[17]);

  RY(QubitReg[1], QuantumVariableParams[18]);
  RX(QubitReg[1], QuantumVariableParams[19]);

  RY(QubitReg[2], QuantumVariableParams[20]);
  RX(QubitReg[2], QuantumVariableParams[21]);

  RY(QubitReg[3], QuantumVariableParams[22]);
  RX(QubitReg[3], QuantumVariableParams[23]);
}

double run_qkernel(FullStateSimulator &iqs_device, const arma::mat &params,
                   SymbolicOperator &symbop,
                   std::map<int, std::set<pstring>> &m_qwc_groups) {
  int basis_change_variable_param_start_indx = 16;
  double total_energy = 0.0;

  QuantumVariableParams[0] = 2 * params[0];
  QuantumVariableParams[1] = 2 * params[1];
  QuantumVariableParams[2] = 2 * params[2];
  QuantumVariableParams[3] = 2 * params[3];
  QuantumVariableParams[4] = 2 * params[4];
  QuantumVariableParams[5] = 2 * params[5];
  QuantumVariableParams[6] = 2 * params[6];
  QuantumVariableParams[7] = 2 * params[7];
  QuantumVariableParams[8] = 2 * params[8];
  QuantumVariableParams[9] = 2 * params[9];
  QuantumVariableParams[10] = 2 * params[10];
  QuantumVariableParams[11] = 2 * params[11];
  QuantumVariableParams[12] = 2 * params[12];
  QuantumVariableParams[13] = 2 * params[13];
  QuantumVariableParams[14] = 2 * params[14];
  QuantumVariableParams[15] = 2 * params[15];

  for (const auto &m_qwc_group : m_qwc_groups) {
    std::vector<double> ProbReg;

    std::vector<double> variable_params(N * 2);

    QuantumVariableParams[16] = 0;
    QuantumVariableParams[17] = 0;
    QuantumVariableParams[18] = 0;
    QuantumVariableParams[19] = 0;
    QuantumVariableParams[20] = 0;
    QuantumVariableParams[21] = 0;
    QuantumVariableParams[22] = 0;
    QuantumVariableParams[23] = 0;

    // Applying basis change to the first hamiltonian will be sufficient
    SymbolicOperatorUtils::applyBasisChange(m_qwc_group.second, variable_params,
                                            N);

    std::vector<std::reference_wrapper<qbit>> qids;
    for (int qubit = 0; qubit < N; ++qubit) {
      qids.push_back(std::ref(QubitReg[qubit]));
    }

    for (auto indx = 0; indx < variable_params.size(); ++indx) {
      QuantumVariableParams[basis_change_variable_param_start_indx + indx] =
          variable_params[indx];
    }

    vqeQ4();

    ProbReg = iqs_device.getProbabilities(qids);

    double current_pstr_val = SymbolicOperatorUtils::getExpectValSetOfPaulis(
        symbop, m_qwc_group.second, ProbReg, N);

    total_energy += current_pstr_val;
  }
  return total_energy;
}

class EnergyOfAnsatz {
public:
  EnergyOfAnsatz(SymbolicOperator &_symbop, QWCMap &_m_qwc_groups,
                 arma::mat &_params, FullStateSimulator &_iqs_device)
      : symbop(_symbop), m_qwc_groups(_m_qwc_groups), params(_params),
        iqs_device(_iqs_device) {}

  double Evaluate(const arma::mat &theta) {
    steps_count++;

    double total_energy = run_qkernel(iqs_device, params, symbop, m_qwc_groups);
    return total_energy;
  }

private:
  const arma::mat &params;
  SymbolicOperator &symbop;
  QWCMap &m_qwc_groups;
  FullStateSimulator &iqs_device;
};

int main() {

  SymbolicOperator H_symbop;
  pstring inp_y1{{0, 'Z'}};
  H_symbop.addTerm(inp_y1, 1.0);
  pstring inp_y2{{0, 'Z'}, {1, 'Z'}};
  H_symbop.addTerm(inp_y2, 1.0);
  pstring inp_y3{{0, 'Z'}, {1, 'Z'}, {2, 'Z'}};
  H_symbop.addTerm(inp_y3, 1.0);
  pstring inp_y4{{0, 'Z'}, {1, 'Z'}, {2, 'Z'}, {3, 'Z'}};
  H_symbop.addTerm(inp_y4, 1.0);
  pstring inp_y5{{2, 'X'}, {3, 'X'}};
  H_symbop.addTerm(inp_y5, 1.0);
  pstring inp_y6{{0, 'Y'}, {2, 'X'}, {3, 'X'}};
  H_symbop.addTerm(inp_y6, 1.0);
  pstring inp_y7{{0, 'Y'}, {1, 'Y'}, {2, 'X'}, {3, 'X'}};
  H_symbop.addTerm(inp_y7, 1.0);

  // Z0 + Z0Z1 + Z0Z1Z2 + Z0Z1Z2Z3 + X2X3 + Y0X2X3 + Y0Y1X2X3
  std::string charstring = H_symbop.getCharString();
  std::cout << "Hamiltonian:\n" << charstring << "\n";

  // QWC
  QWCMap m_qwc_groups =
      SymbolicOperatorUtils::getQubitwiseCommutationGroups(H_symbop, N);
  std::cout << "Number of Qubitwise Commutation (QWC) Groups : "
            << m_qwc_groups.size() << std::endl;
  std::cout << "Qubitwise Commutation (QWC) Groups: " << m_qwc_groups
            << std::endl;

  /// Setup quantum device
  IqsConfig iqs_config(/*num_qubits*/ N,
                        /*simulation_type*/ "noiseless");
  FullStateSimulator iqs_device(iqs_config);
  if (QRT_ERROR_SUCCESS != iqs_device.ready()) {
    return -1;
  }

  // Using Simulated Annealing algorithm
  // Starting parameters
  arma::mat params_sa(std::vector<double>{
      0.261799, 0.261799, 0.261799, 0.261799, 0.261799, 0.261799, 0.261799,
      0.261799, 0.261799, 0.261799, 0.261799, 0.261799, 0.261799, 0.261799,
      0.261799, 0.261799});

  QuantumVariableParams[0] = params_sa[0];
  QuantumVariableParams[1] = params_sa[1];
  QuantumVariableParams[2] = params_sa[2];
  QuantumVariableParams[3] = params_sa[3];
  QuantumVariableParams[4] = params_sa[4];
  QuantumVariableParams[5] = params_sa[5];
  QuantumVariableParams[6] = params_sa[6];
  QuantumVariableParams[7] = params_sa[7];
  QuantumVariableParams[8] = params_sa[8];
  QuantumVariableParams[9] = params_sa[9];
  QuantumVariableParams[10] = params_sa[10];
  QuantumVariableParams[11] = params_sa[11];
  QuantumVariableParams[12] = params_sa[12];
  QuantumVariableParams[13] = params_sa[13];
  QuantumVariableParams[14] = params_sa[14];
  QuantumVariableParams[15] = params_sa[15];

  // arbitrary function
  EnergyOfAnsatz eoa_sa(H_symbop, m_qwc_groups, params_sa, iqs_device);

  // initialize the optimization algorithm
  ens::SA<> opt_sa(
      ens::ExponentialSchedule(), /* coolingSchedule - Instantiated cooling
                                     schedule (default ExponentialSchedule).
                                     */
      1000000, /* maxIterations - Maximum number of iterations allowed (0
                  indicates no limit). */
      1000.,   /* initT - Initial temperature. */
      1000,    /* initMoves - Number of initial iterations without changing
                  temperature. */
      100,     /* moveCtrlSweep - Sweeps per feedback move control. */
      1e-10,   /* tolerance - Tolerance to consider system frozen. */
      3,       /* maxToleranceSweep - Maximum sweeps below tolerance to
      consider
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

  // Starting parameters
  params_sa = std::vector<double>{0.261799, 0.261799, 0.261799, 0.261799,
                                  0.261799, 0.261799, 0.261799, 0.261799,
                                  0.261799, 0.261799, 0.261799, 0.261799,
                                  0.261799, 0.261799, 0.261799, 0.261799};

  QuantumVariableParams[0] = params_sa[0];
  QuantumVariableParams[1] = params_sa[1];
  QuantumVariableParams[2] = params_sa[2];
  QuantumVariableParams[3] = params_sa[3];
  QuantumVariableParams[4] = params_sa[4];
  QuantumVariableParams[5] = params_sa[5];
  QuantumVariableParams[6] = params_sa[6];
  QuantumVariableParams[7] = params_sa[7];
  QuantumVariableParams[8] = params_sa[8];
  QuantumVariableParams[9] = params_sa[9];
  QuantumVariableParams[10] = params_sa[10];
  QuantumVariableParams[11] = params_sa[11];
  QuantumVariableParams[12] = params_sa[12];
  QuantumVariableParams[13] = params_sa[13];
  QuantumVariableParams[14] = params_sa[14];
  QuantumVariableParams[15] = params_sa[15];

  std::cout
      << "Expectation value of a single quantum kernel call w/ static params: ";
  std::cout << run_qkernel(iqs_device, params_sa, H_symbop, m_qwc_groups)
            << "\n";

  return 0;
}
