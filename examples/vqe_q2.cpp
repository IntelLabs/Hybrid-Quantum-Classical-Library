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
#include <vector>

#include <armadillo>
#include <ensmallen.hpp>

const int N = 2;
qbit QubitReg[N];

/* Special global vector from QRT to get measurement results */
extern std::vector<bool> ClassicalBitRegister;
/* Special global vector from QRT to get state probabilities */
extern std::vector<double> ProbabilityRegister;

/* Special global array to hold dynamic parameters for quantum algorithm */
quantum_shared_double_array QuantumVariableParams[4];

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

  // Measurements of all the qubits
  for (Index = 0; Index < N; Index++) {
    MeasZ(QubitReg[Index]);
  }
}

quantum_kernel void vqeQ2_X() {

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

  // Measurements of all the qubits
  for (Index = 0; Index < N; Index++) {
    H(QubitReg[Index]); // X-Basis to Z-Basis
    MeasZ(QubitReg[Index]);
  }
}


// Calculate expectation value on Z-basis
double Z0Z1_expectation(std::vector<double> ProbReg) {
  auto Z0Z1 = ProbReg[0] - ProbReg[1] - ProbReg[2] + ProbReg[3];
  return Z0Z1;
}

// Calculate expectation value on Z-basis
double Z0_expectation(std::vector<double> ProbReg) {
  auto Z0 = ProbReg[0] + ProbReg[1] - ProbReg[2] - ProbReg[3];
  return Z0;
}

// Calculate expectation value on Z-basis
double Z1_expectation(std::vector<double> ProbReg) {
  auto Z1 = ProbReg[0] - ProbReg[1] + ProbReg[2] - ProbReg[3];
  return Z1;
}

double run_qkernel(const arma::mat &params) {
    std::vector<double> ProbRegZ;
    std::vector<double> ProbRegX;

    QuantumVariableParams[0] = 2 * params[0];
    QuantumVariableParams[1] = 2 * params[1];
    QuantumVariableParams[2] = 2 * params[2];
    QuantumVariableParams[3] = 2 * params[3];

    vqeQ2();
    for (auto j = 0; j < ProbabilityRegister.size(); j++) {
      auto p = ProbabilityRegister[j];
      if (p > 0) {
        ProbRegZ.push_back(ProbabilityRegister[j]);
        std::cout << j << " : " << ProbabilityRegister[j] << "\t";
        if ((j + 1) % 8 == 0)
          std::cout << "\n";
      }
    }
    std::cout << "\n";

    vqeQ2_X();
    for (auto j = 0; j < ProbabilityRegister.size(); j++) {
      auto p = ProbabilityRegister[j];
      if (p > 0) {
        ProbRegX.push_back(ProbabilityRegister[j]);
        std::cout << j << " : " << ProbabilityRegister[j] << "\t";
        if ((j + 1) % 8 == 0)
          std::cout << "\n";
      }
    }
    std::cout << "\n";

    double energy_Z0Z1 = Z0Z1_expectation(ProbRegZ);

    double energy_Z0 = Z0_expectation(ProbRegX);

    double energy_Z1 = Z1_expectation(ProbRegX);

    double total_energy =
        0.5 * energy_Z0Z1 + 0.5 * energy_Z0 + 0.25 * energy_Z1;

    return total_energy;
}

class EnergyOfAnsatz {
public:
  EnergyOfAnsatz(arma::mat &params) : params(params) {}

  double Evaluate(const arma::mat &theta) {
    steps_count++;

    return run_qkernel(params);
  }

private:
  const arma::mat &params;
};

int main() {

  // Using SPSA algorithm
  // Starting parameters
  arma::mat params_spsa(std::vector<double>{0.1, 0.1, 0.2, 0.1});

  QuantumVariableParams[0] = params_spsa[0];
  QuantumVariableParams[1] = params_spsa[1];
  QuantumVariableParams[2] = params_spsa[2];
  QuantumVariableParams[3] = params_spsa[3];

  // arbitrary function
  EnergyOfAnsatz eoa(params_spsa);

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
  EnergyOfAnsatz eoa_sa(params_sa);

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
