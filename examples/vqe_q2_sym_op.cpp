//===----------------------------------------------------------------------===//
// INTEL CONFIDENTIAL
//
// Copyright 2021 Intel Corporation.
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
// #include <clang/Quantum/quintrinsics.h>

/// Development mode
#include "../../clang/include/clang/Quantum/quintrinsics.h"
#include <cassert>
#include <iostream>
#include <qrt0.hpp>
#include <sstream>
#include <string>
#include <vector>

#include <armadillo>
#include <ensmallen.hpp>

#include "SymbolicOperator.hpp"
#include "SymbolicOperatorUtils.hpp"

using namespace hybrid::quantum::core;

const int N = 2;
qbit QubitReg[N];
cbit CReg[N];

const double FP_PI = 3.14159265359;
const double FP_2PI = 6.28318530718;
const double FP_PIby2 = 1.57079632679489661923;

/* Special global vector from QRT to get state probabilities */
extern std::vector<double> ProbabilityRegister;

/* Special global array to hold dynamic parameters for quantum algorithm */
quantum_shared_double_array QuantumVariableParams[4 + N * 2];

static int steps_count = 0;

quantum_kernel void vqeQ2() {

  // Index to loop over later
  int Index = 0;

  // Initialization of the qubits
  for (Index = 0; Index < N; Index++) {
    PrepZ(QubitReg[Index]);
  }

  // Ansatz
  RY(QubitReg[0], QuantumVariableParams[0]);
  RY(QubitReg[1], QuantumVariableParams[1]);

  CNOT(QubitReg[0], QubitReg[1]);

  RY(QubitReg[0], QuantumVariableParams[2]);
  RY(QubitReg[1], QuantumVariableParams[3]);

  // not part of ansatz
  RY(QubitReg[0], QuantumVariableParams[4]);
  RX(QubitReg[0], QuantumVariableParams[5]);

  RY(QubitReg[1], QuantumVariableParams[6]);
  RX(QubitReg[1], QuantumVariableParams[7]);

  // Measurements of all the qubits
  for (Index = 0; Index < N; Index++) {
    MeasZ(QubitReg[Index], CReg[Index]);
  }
}

double expectation_value(const pstring &pstr,
                         const std::vector<double> ProbReg) {
  arma::Row<double> I = {1, 1};
  arma::Row<double> Z = {1, -1};
  arma::rowvec rv(ProbReg);
  double exp_val;
  arma::Row<double> vR = {1};

  for (const auto &ps : pstr) {
    if (ps.second == 'X' || ps.second == 'Y' || ps.second == 'Z')
      vR = arma::kron(vR, Z).as_row();
    else
      vR = arma::kron(vR, I).as_row();
  }

  exp_val = arma::sum((rv % vR).as_row());

  return exp_val;
}

double run_qkernel(const arma::mat &params, SymbolicOperator &so) {
  int basis_change_start = 4;
  char first_basis_change;
  double total_energy = 0.0;

  QuantumVariableParams[0] = 2 * params[0];
  QuantumVariableParams[1] = 2 * params[1];
  QuantumVariableParams[2] = 2 * params[2];
  QuantumVariableParams[3] = 2 * params[3];

  for (const auto &o : so.getOrderedPStringList()) {
    // Find first basis change
    for (const auto &ps : o) {
      if (ps.second == 'X' || ps.second == 'Y' || ps.second == 'Z') {
        first_basis_change = ps.second;
        break;
      }
    }

    std::vector<double> ProbReg;

    QuantumVariableParams[4] = 0;
    QuantumVariableParams[5] = 0;
    QuantumVariableParams[6] = 0;
    QuantumVariableParams[7] = 0;
    for (const auto &ps : o) {
      if (ps.second == 'X') {
        QuantumVariableParams[basis_change_start + 2 * ps.first] = FP_PIby2;
        QuantumVariableParams[basis_change_start + 2 * ps.first + 1] = FP_PI;
      } else if (ps.second == 'Y') {
        QuantumVariableParams[basis_change_start + 2 * ps.first] = FP_PI;
        QuantumVariableParams[basis_change_start + 2 * ps.first + 1] = FP_PIby2;
      } else if (ps.second == 'Z') {
        QuantumVariableParams[basis_change_start + 2 * ps.first] = 0;
        QuantumVariableParams[basis_change_start + 2 * ps.first + 1] = 0;
      } else if ((ps.second != 'X' || ps.second != 'Y' || ps.second != 'Z')) {
        if (first_basis_change == 'X') {
          QuantumVariableParams[basis_change_start + 2 * ps.first] = FP_PIby2;
          QuantumVariableParams[basis_change_start + 2 * ps.first + 1] = FP_PI;
        } else if (first_basis_change == 'Y') {
          QuantumVariableParams[basis_change_start + 2 * ps.first] = FP_PI;
          QuantumVariableParams[basis_change_start + 2 * ps.first + 1] =
              FP_PIby2;
        } else if (first_basis_change == 'Z') {
          QuantumVariableParams[basis_change_start + 2 * ps.first] = 0;
          QuantumVariableParams[basis_change_start + 2 * ps.first + 1] = 0;
        }
      }
    }

    vqeQ2();
    for (auto j = 0; j < ProbabilityRegister.size(); j++) {
      auto p = ProbabilityRegister[j];
      if (p > 0) {
        ProbReg.push_back(ProbabilityRegister[j]);
        std::cout << j << " : " << ProbabilityRegister[j] << "\t";
        if ((j + 1) % 8 == 0)
          std::cout << "\n";
      }
    }
    std::cout << "\n";

    double current_pstr_val =
        so.op_sum[o].real() * SymbolicOperatorUtils::getExpectValSglPauli(o, ProbReg);
    total_energy += current_pstr_val;
  }

  return total_energy;
}

class EnergyOfAnsatz {
public:
  EnergyOfAnsatz(SymbolicOperator &so, arma::mat &params)
      : so(so), params(params) {}

  double Evaluate(const arma::mat &theta) {
    steps_count++;

    return run_qkernel(params, so);
  }

private:
  const arma::mat &params;
  SymbolicOperator so;
};

int main() {

  SymbolicOperator so;
  pstring inp_y1{{0, 'Z'}, {1, 'Z'}};
  so.addTerm(inp_y1, 0.5);
  pstring inp_y2{{0, 'X'}, {1, 'I'}};
  so.addTerm(inp_y2, 0.5);
  pstring inp_y3{{0, 'I'}, {1, 'X'}};
  so.addTerm(inp_y3, 0.25);
  std::string charstring = so.getCharString();
  std::cout << "Hamiltonian:\n" << charstring << "\n";

  // Using SPSA algorithm
  // Starting parameters
  arma::mat params_spsa(std::vector<double>{0.1, 0.1, 0.2, 0.1});

  QuantumVariableParams[0] = params_spsa[0];
  QuantumVariableParams[1] = params_spsa[1];
  QuantumVariableParams[2] = params_spsa[2];
  QuantumVariableParams[3] = params_spsa[3];

  // arbitrary function
  EnergyOfAnsatz eoa(so, params_spsa);

  // initialize the optimization algorithm
  ens::SPSA opt_spsa(
      0.2,   /* alpha - Scaling exponent for the step size. */
      0.101, /* gamma - Scaling exponent for evaluation step size. */
      0.16, /* stepSize - Scaling parameter for step size (named as ‘a’ in the
               paper). */
      0.3,    /* evaluationStepSize - Scaling parameter for evaluation step size
                 (named as ‘c’ in the paper). */
      100000, /* maxIterations - Maximum number of iterations allowed (0 means
                 no limit). */
      1e-10 /* tolerance - Maximum absolute tolerance to terminate algorithm. */
  );

  // optimize the parameters
  double opt_spsa_val = opt_spsa.Optimize(eoa, params_spsa);
  int spsa_steps_count = steps_count;

  // resetting steps count
  steps_count = 0;

  // Using SA algorithm
  // Starting parameters
  arma::mat params_sa(std::vector<double>{0.1, 0.1, 0.2, 0.1});

  QuantumVariableParams[0] = params_sa[0];
  QuantumVariableParams[1] = params_sa[1];
  QuantumVariableParams[2] = params_sa[2];
  QuantumVariableParams[3] = params_sa[3];

  // arbitrary function
  EnergyOfAnsatz eoa_sa(so, params_sa);

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

  std::cout << "Optimized Parameters using Simultaneous Perturbation "
               "Stochastic Approximation (SPSA): "
            << "\n";
  params_spsa.print();
  std::cout << "Simultaneous Perturbation Stochastic Approximation (SPSA) "
               "execution count: "
            << spsa_steps_count << "\n";
  std::cout << "Total energy using Simultaneous Perturbation Stochastic "
               "Approximation(SPSA): "
            << opt_spsa_val << "\n";

  std::cout << "Optimized Parameters using Simulated Annealing (SA): "
            << "\n";
  params_sa.print();
  std::cout << "Simulated Annealing (SA) execution count: " << sa_steps_count
            << "\n";
  std::cout << "Total energy using Simulated Annealing (SA): " << opt_sa_val
            << "\n";

  return 0;
}
