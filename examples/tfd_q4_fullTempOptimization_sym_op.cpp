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

#include "SymbolicOperator.hpp"
#include "SymbolicOperatorUtils.hpp"

using namespace hybrid::quantum::core;

// ----------------------------------------------------------------------------------------
// ALL QUANTUM HELPER CODE
// ----------------------------------------------------------------------------------------

// Define the number of qubits needed for compilation
const int N = 4;
qbit QubitReg[N];
cbit CReg[N];

const double FP_PI = 3.14159265359;
const double FP_2PI = 6.28318530718;
const double FP_PIby2 = 1.57079632679;

/* Special global vector from QRT to get state probabilities */
extern std::vector<double> ProbabilityRegister;

/* Special global array to hold dynamic parameters for quantum algorithm */
/* alpha1, alpha2, gamma1, gamma2 is the order of the variational angles used in
 * here */
quantum_shared_double_array QuantumVariableParams[N + N * 2];

static int steps_count = 0;

/*
When using the Intel Quantum Simulator (IQS) at the lowest level, it is
necessary to decide how to calculate the cost function based on the results
returned. If using the full state vector, then only one experiment is required
per iteration.

In this script, we will be using the probabilities corresponding to each of the
basis states (e.g. 0000, 1010) These are returned simply by calculating the
probability corresponding to each of the state vector element and stored in the
special array ProbabilityRegister.
*/

//
// First kind of experiment that needs to be run per iteration of optimization.
// Contains the basic quantum operations as described in Fig. 1 of Ref. 1 with
// the appropriate decompositions such that the Intel Quantum Compiler can
// handle them. Further an additional 2 layers of rotations are added to help
// enable mapping from X-basis to Z-basis or Y-basis to Z-basis, so that X or Y
// related observables can be calculated.
//
quantum_kernel void tfdQ4() {
  // Index to loop over later
  int Index = 0;

  // Initialization of the qubits
  for (Index = 0; Index < N; Index++)
    PrepZ(QubitReg[Index]);

  // preparation of Bell pairs (T -> Infinity)
  for (Index = 0; Index < 2; Index++)
    RY(QubitReg[Index], FP_PIby2);
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
    RY(QubitReg[Index], -FP_PIby2);
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index + 2], QubitReg[Index]);
  for (Index = 0; Index < 2; Index++)
    RZ(QubitReg[Index], QuantumVariableParams[0]);
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index + 2], QubitReg[Index]);
  for (Index = 0; Index < N; Index++)
    RY(QubitReg[Index], FP_PIby2);

  // two-qubit inter-system ZZ variational terms
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index], QubitReg[Index + 2]);
  for (Index = 0; Index < 2; Index++)
    RZ(QubitReg[Index + 2], QuantumVariableParams[1]);
  for (Index = 0; Index < 2; Index++)
    CNOT(QubitReg[Index], QubitReg[Index + 2]);

  // not part of ansatz but additional layers to help
  // with X to Z or Y to Z basis mapping
  RY(QubitReg[0], QuantumVariableParams[4]);
  RX(QubitReg[0], QuantumVariableParams[5]);

  RY(QubitReg[1], QuantumVariableParams[6]);
  RX(QubitReg[1], QuantumVariableParams[7]);

  RY(QubitReg[2], QuantumVariableParams[8]);
  RX(QubitReg[2], QuantumVariableParams[9]);

  RY(QubitReg[3], QuantumVariableParams[10]);
  RX(QubitReg[3], QuantumVariableParams[11]);

  // Measurements of all the qubits
  for (Index = 0; Index < N; Index++)
    MeasZ(QubitReg[Index], CReg[Index]);
}

double run_qkernel(const arma::mat &params, SymbolicOperator &symbop) {
  // start point in the global parameter array corresponding to the basis change
  int basis_change_variable_param_start_indx = 4;
  // total cost
  double total_cost = 0.0;

  QuantumVariableParams[0] = 2 * params[0];
  QuantumVariableParams[1] = 2 * params[1];
  QuantumVariableParams[2] = 2 * params[2];
  QuantumVariableParams[3] = 2 * params[3];

  // loop through all the Pauli strings present in the SymbolicOperator object
  for (const auto &pstr : symbop.getOrderedPStringList()) {
    std::vector<double> ProbReg;

    QuantumVariableParams[4] = 0;
    QuantumVariableParams[5] = 0;
    QuantumVariableParams[6] = 0;
    QuantumVariableParams[7] = 0;
    QuantumVariableParams[8] = 0;
    QuantumVariableParams[9] = 0;
    QuantumVariableParams[10] = 0;
    QuantumVariableParams[11] = 0;

    // Update the global parameters array for the additional layers
    // used to enable mapping of X-basis to Z-basis or Y-basis to
    // Z-basis
    std::vector<double> variable_params;
    variable_params.reserve(N * 2);
    SymbolicOperatorUtils::applyBasisChange(pstr, variable_params, N);

    for (auto indx = 0; indx < variable_params.size(); ++indx) {
      QuantumVariableParams[basis_change_variable_param_start_indx + indx] = variable_params[indx];
    }

    // performing the experiment, and storing the data in ProbReg
    tfdQ4();
    std::cout << "Probability Register Size: " << ProbabilityRegister.size()
              << "\n";
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

    // calculate the expectation value
    double current_pstr_val =
        symbop.op_sum[pstr].real() * SymbolicOperatorUtils::getExpectValSglPauli(pstr, ProbReg, N);

    total_cost += current_pstr_val;
  }

  return total_cost;
}

// Constructing a lambda function to be used for a single optimization iteration
// This function is directly called by the dlib optimization routine
double ansatz_run_lambda(SymbolicOperator &symbop, const arma::mat &params,
                         double inv_temp) {
  // loading the new variational angles into the special global array for hybrid
  // compilation
  QuantumVariableParams[0] = params(0);
  QuantumVariableParams[1] = params(1);
  QuantumVariableParams[2] = params(2);
  QuantumVariableParams[3] = params(3);

  // runs the kernel to compute the total cost
  double total_cost = run_qkernel(params, symbop);

  return total_cost;
}

class EnergyOfAnsatz {
public:
  EnergyOfAnsatz(SymbolicOperator &symbop_, arma::mat &params_, double beta_)
      : symbop(symbop_), params(params_), beta(beta_) {}

  double Evaluate(const arma::mat &theta) {
    steps_count++;

    double total_cost;
    total_cost = ansatz_run_lambda(symbop, params, beta);

    std::cout << "step: " << beta << ", ";
    for (int i = 0; i < 4; i++)
      std::cout << params[i] << ", ";
    std::cout << total_cost << "\n";

    return total_cost;
  }

private:
  SymbolicOperator symbop;
  const arma::mat &params;
  double beta;
};

int main() {
  // H = X3 + X2 + X1 + X0 + 1.60 (Z3Z2 + Z1Z0) - (beta ^ -1.48((X1X3 + X0X2)+(Z1Z3 +
  // Z0Z2)))

  // resetting steps count
  steps_count = 0;

  // Clearing the file contents used to store all the results from TFD
  // minimization
  std::ofstream file_to_clear;
  file_to_clear.open("tfd_minimization_results_sym_op.csv",
                     std::ofstream::out | std::ofstream::trunc);
  file_to_clear.close();

  std::ofstream results_file;
  results_file.open("tfd_minimization_results_sym_op.csv", std::ios_base::app);
  for (int beta_index = -30; beta_index < 31; beta_index += 1) {
    // calculating the actual inverse temperature that is used during
    // calculations
    double beta = pow(10, (double)beta_index / 10);

    // construct object
    SymbolicOperator symbop;
    pstring Xa1{{0, 'X'}};
    symbop.addTerm(Xa1, 1);
    pstring Xa2{{1, 'X'}};
    symbop.addTerm(Xa2, 1);
    pstring Xb1{{2, 'X'}};
    symbop.addTerm(Xb1, 1);
    pstring Xb2{{3, 'X'}};
    symbop.addTerm(Xb2, 1);

    pstring Za{{0, 'Z'}, {1, 'Z'}};
    symbop.addTerm(Za, 1.60);
    pstring Zb{{2, 'Z'}, {3, 'Z'}};
    symbop.addTerm(Zb, 1.60);

    pstring Xa1b1{{0, 'X'}, {2, 'X'}};
    symbop.addTerm(Xa1b1, -pow(beta, -1.48));
    pstring Xa2b2{{1, 'X'}, {3, 'X'}};
    symbop.addTerm(Xa2b2, -pow(beta, -1.48));

    pstring Za1b1{{0, 'Z'}, {2, 'Z'}};
    symbop.addTerm(Za1b1, -pow(beta, -1.48));
    pstring Za2b2{{1, 'Z'}, {3, 'Z'}};
    symbop.addTerm(Za2b2, -pow(beta, -1.48));

    std::cout << "SymbolicOperator Object: " << symbop.getCharString() << "\n";

    // Using SA algorithm
    // Starting parameters
    arma::mat params_sa(std::vector<double>{0, 0, 0, 0});

    // arbitrary function
    EnergyOfAnsatz eoa_sa(symbop, params_sa, beta);

    // initialize the optimization algorithm
    ens::SA<> opt_sa(
        ens::ExponentialSchedule(), /* coolingSchedule - Instantiated cooling
                                       schedule (default ExponentialSchedule).
                                     */
        1000000, /* maxIterations - Maximum number of iterations allowed (0
                    indicates no limit). */
        1000.,   /* initT - Initial temperature. */
        3000,    /* initMoves - Number of initial iterations without changing
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
