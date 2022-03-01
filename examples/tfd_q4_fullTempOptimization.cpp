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
// TFD
// A single-step version of an algorithm to generate 4-qubit Thermofield Double
// States (TFD) which are purifications of Gibbs states.
//
// References:
//
// Theory:
//
// 1. Wu, J. & Hsieh, T. H. Variational thermal quantum simulation via
// thermofield double states. Phys. Rev. Lett. 123, 220502 (2019).
//
// 2. Ho, W. W. & Hsieh, T. H. Efficient variational simulation of non-trivial
// quantum states. SciPost Phys. 6, 29 (2019).
//
// 3. Premaratne, S. P. & Matsuura, A. Y. Engineering a cost function for
// real-world implementation of a variational quantum algorithm. In Proc. 2020
// IEEE Int. Conf. Quantum Comp. Eng., 278{285 (2020).
//
// Experiment:
//
// 1. Sagastizabal, R. et al. Variational preparation of  finite-temperature
// states on a quantum computer (2021). npj Quantum Inf. 7:130 (2021)
// https://doi.org/10.1038/s41534-021-00468-1.
//
// 2. Zhu, D. et al. Generation of thermofield double states and critical ground
// states with a quantum computer. P. Natl Acad. Sci. USA 117, 25402{25406
// (2020). https://doi.org/10.1073/pnas.2006337117
//
// 3. Premaratne, S. et al. Engineering a Cost Function for Real-world
// Implementation of a Variational Quantum Algorithm (2020).
// https://doi.org/10.1109/QCE49297.2020.00042.
//===----------------------------------------------------------------------===//

/// NOTE: Pick an include path depending on whether running in development
/// enviromnment or production environment
//===----------------------------------------------------------------------===//
/// Production mode
// #include <clang/Quantum/quintrinsics.h>

// AUTHOR : Shavindra Premaratne (August 26th 2021)
// Update : Comments added throughout to clarify the intentions behind code
// blocks (08/27/2021)

/// Development mode
#include "../../clang/include/clang/Quantum/quintrinsics.h"

#include <cassert> // to assert the size of the Probability Register
#include <fstream> // to write to files
#include <iostream>
#include <math.h> // to use pow()
#include <vector>

#include <qrt0.hpp>

// to perform the minimization of TFD generation cost funtion

#include <armadillo>
#include <ensmallen.hpp>

// ----------------------------------------------------------------------------------------
// ALL QUANTUM HELPER CODE
// ----------------------------------------------------------------------------------------

// Define the number of qubits needed for compilation
const int N = 4;
qbit QubitReg[N];
cbit CReg[N];

// Special global vector from QRT to get state probabilities
extern std::vector<double> ProbabilityRegister;

// Special global array to hold dynamic parameters for quantum algorithm
// alpha1, alpha2, gamma1, gamma2 is the order of the variational angles used in
// here
quantum_shared_double_array QuantumVariableParams[N];

static int steps_count = 0;

// When using the Intel Quantum Simulator (IQS) at the lowest level, it is
// necessary to decide how to calculate the cost function based on the results
// returned. If using the full state vector, then only one experiment is
// required per iteration.
//
// In this script, we will be using the probabilities corresponding to each of
// the basis states (e.g. 0000, 1010) These are returned simply by calculating
// the probability corresponding to each of the state vector element and stored
// in the special array ProbabilityRegister.

// First kind of experiment that needs to be run per iteration of optimization.
// Contains the basic quantum operations as described in Fig. 1 of Ref. 1
// with the appropriate decompositions such that the Intel Quantum Compiler can
// handle them.

quantum_kernel void tfdQ4_Z() {
  // Index to loop over later
  int Index = 0;

  // Initialization of the qubits
  for (Index = 0; Index < N; Index++)
    PrepZ(QubitReg[Index]);

  // preparation of Bell pairs (T -> Infinity)
  for (Index = 0; Index < 2; Index++)
    RY(QubitReg[Index], 1.57079632679);
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index], QubitReg[Index + 2]);

  // Single qubit variational terms
  for (Index = 0; Index < N; Index++)
    RX(QubitReg[Index], QuantumVariableParams[2]);

  // Two-qubit intra-system variational terms
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[2 * Index + 1], QubitReg[2 * Index]);
  for (Index = 0; Index < 2; Index++)
    RZ(QubitReg[2 * Index], QuantumVariableParams[3]);
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[2 * Index + 1], QubitReg[2 * Index]);

  // two-qubit inter-system XX variational terms
  for (Index = 0; Index < N; Index++)
    RY(QubitReg[Index], -1.57079632679);
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index + 2], QubitReg[Index]);
  for (Index = 0; Index < 2; Index++)
    RZ(QubitReg[Index], QuantumVariableParams[0]);
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index + 2], QubitReg[Index]);
  for (Index = 0; Index < N; Index++)
    RY(QubitReg[Index], 1.57079632679);

  // two-qubit inter-system ZZ variational terms
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index], QubitReg[Index + 2]);
  for (Index = 0; Index < 2; Index++)
    RZ(QubitReg[Index + 2], QuantumVariableParams[1]);
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index], QubitReg[Index + 2]);

  // Measurements of all the qubits
  for (Index = 0; Index < N; Index++)
    MeasZ(QubitReg[Index], CReg[Index]);
}

// Second kind of experiment that needs to be run per iteration of optimizaiton.
// Almost identical to the previous quantum kernel tfdQ4_Z, in this qkernel,
// four Hadamard operations are inserted just prior to measurement of the
// qubits. This enables mapping from X-basis to Z-basis, so that X related
// observables can be calculated.

quantum_kernel void tfdQ4_X() {
  // Index to loop over later
  int Index = 0;

  // Initialization of the qubits
  for (Index = 0; Index < N; Index++)
    PrepZ(QubitReg[Index]);

  // preparation of Bell pairs (T -> Infinity)
  for (Index = 0; Index < 2; Index++)
    RY(QubitReg[Index], 1.57079632679);
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index], QubitReg[Index + 2]);

  // Single qubit variational terms
  for (Index = 0; Index < N; Index++)
    RX(QubitReg[Index], QuantumVariableParams[2]);

  // Two-qubit intra-system variational terms
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[2 * Index + 1], QubitReg[2 * Index]);
  for (Index = 0; Index < 2; Index++)
    RZ(QubitReg[2 * Index], QuantumVariableParams[3]);
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[2 * Index + 1], QubitReg[2 * Index]);

  // two-qubit inter-system XX variational terms
  for (Index = 0; Index < N; Index++)
    RY(QubitReg[Index], -1.57079632679);
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index + 2], QubitReg[Index]);
  for (Index = 0; Index < 2; Index++)
    RZ(QubitReg[Index], QuantumVariableParams[0]);
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index + 2], QubitReg[Index]);
  for (Index = 0; Index < N; Index++)
    RY(QubitReg[Index], 1.57079632679);

  // two-qubit inter-system ZZ variational terms
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index], QubitReg[Index + 2]);
  for (Index = 0; Index < 2; Index++)
    RZ(QubitReg[Index + 2], QuantumVariableParams[1]);
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index], QubitReg[Index + 2]);

  // Mapping from X basis to Z basis
  for (Index = 0; Index < N; Index++)
    H(QubitReg[Index]);

  // Measurements of all the qubits
  for (Index = 0; Index < N; Index++)
    MeasZ(QubitReg[Index], CReg[Index]);
}

// The three expectation values of the observables have been written as
// functions below. The relationships between the expectation values and the
// probabilities obtained from the state vector amplitudes was derived
// independently with simple matrix calculations of Pauli matrix products.

double ZZApZZB_expectation(std::vector<double> ProbReg) {
  auto line0 = ProbReg[0] + ProbReg[3] - ProbReg[5] - ProbReg[6] - ProbReg[9] -
               ProbReg[10] + ProbReg[12] + ProbReg[15];
  return 2 * line0;
}

double ZZAB_expectation(std::vector<double> ProbReg) {
  auto line0 = ProbReg[0] - ProbReg[3] + ProbReg[5] - ProbReg[6] - ProbReg[9] +
               ProbReg[10] - ProbReg[12] + ProbReg[15];
  return 2 * line0;
}

double ZApZB_expectation(std::vector<double> ProbReg) {
  auto line0 = ProbReg[0] - ProbReg[15];
  auto line1 = ProbReg[1] + ProbReg[2] + ProbReg[4] - ProbReg[7] + ProbReg[8] -
               ProbReg[11] - ProbReg[13] - ProbReg[14];
  return (4 * line0 + 2 * line1);
}

// The total cost function is constructed here using the expectation values from
// before. There is freedom in the numeric coefficient and index used. Here we
// use the best engineered simple cost function from Ref. 3

double total_cost(double inv_temp, std::vector<double> ProbReg_Z,
                  std::vector<double> ProbReg_X) {
  double energy_X = ZApZB_expectation(ProbReg_X);
  double energy_ZZ = ZZApZZB_expectation(ProbReg_Z);
  double entropy = ZZAB_expectation(ProbReg_X) + ZZAB_expectation(ProbReg_Z);
  return energy_X + 1.60 * energy_ZZ - pow(inv_temp, -1.48) * entropy;
}

// Constructing a lambda function to be used for a single optimization iteration
// This function is directly called by the dlib optimization routine
double ansatz_run_lambda(const arma::mat &m, double inv_temp) {
  // loading the new variational angles into the special global array for hybrid
  // compilation
  QuantumVariableParams[0] = m(0);
  QuantumVariableParams[1] = m(1);
  QuantumVariableParams[2] = m(2);
  QuantumVariableParams[3] = m(3);

  // two local variables to store the results from the two different experiments
  // run in each iteration
  std::vector<double> ProbRegZ;
  std::vector<double> ProbRegX;

  // performing the Z-focused experiment, and storing the data in ProbRegZ
  tfdQ4_Z();
  std::cout << "Z - Probability Register Size: " << ProbabilityRegister.size()
            << "\n";
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

  // performing the X-focused experiment, and storing the data in ProbRegZ
  tfdQ4_X();
  std::cout << "X - Probability Register Size: " << ProbabilityRegister.size()
            << "\n";
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

  // calculating the cost function for this iteration
  double calc_cost = total_cost(inv_temp, ProbRegZ, ProbRegX);
  return calc_cost;
}

class EnergyOfAnsatz {
public:
  EnergyOfAnsatz(arma::mat &params, double beta) : params(params), beta(beta) {}

  double Evaluate(const arma::mat &theta) {
    steps_count++;

    double total_cost;
    total_cost = ansatz_run_lambda(params, beta);

    std::cout << "step: " << beta << ", ";
    for (int i = 0; i < 4; i++)
      std::cout << params[i] << ", ";
    std::cout << total_cost << "\n";

    return total_cost;
  }

private:
  const arma::mat &params;
  double beta;
};

int main() {
  // resetting steps count
  steps_count = 0;

  // resetting steps count
  steps_count = 0;

  // Clearing the file contents used to store all the results from TFD
  // minimization
  std::ofstream file_to_clear;
  file_to_clear.open("tfd_minimization_results.csv",
                     std::ofstream::out | std::ofstream::trunc);
  file_to_clear.close();

  std::ofstream results_file;
  results_file.open("tfd_minimization_results.csv", std::ios_base::app);
  for (int beta_index = -30; beta_index < 31; beta_index += 1) {
    // calculating the actual inverse temperature that is used during
    // calculations
    double beta = pow(10, (double)beta_index / 10);

    // Using SA algorithm
    // Starting parameters
    arma::mat params_sa(std::vector<double>{0, 0, 0, 0});

    // arbitrary function
    EnergyOfAnsatz eoa_sa(params_sa, beta);

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
        1e-7,    /* tolerance - Tolerance to consider system frozen. */
        3,   /* maxToleranceSweep - Maximum sweeps below tolerance to consider
                system frozen. */
        1.5, /* maxMoveCoef - Maximum move size. */
        0.5, /* initMoveCoef - Initial move size. */
        0.3  /* gain - Proportional control in feedback move control. */
    );

    // optimize the parameters
    double opt_sa_val = opt_sa.Optimize(eoa_sa, params_sa);
    int sa_steps_count = steps_count;
    std::cout << "SA Step Count: " << sa_steps_count << "\n";

    // writing the results from a single temperature to file
    results_file << beta << ", ";
    for (int i = 0; i < 4; i++)
      results_file << params_sa[i] << ", ";
    results_file << opt_sa_val << "\n";
  }

  results_file.close();

  return 0;
}

