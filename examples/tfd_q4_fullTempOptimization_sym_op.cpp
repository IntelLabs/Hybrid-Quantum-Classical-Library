//===----------------------------------------------------------------------===//
// INTEL CONFIDENTIAL
//
// Copyright 2021-2024 Intel Corporation.
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

//
// Expected Output
// $ export HQ=<path to project directory>
// $ ./intel-quantum-compiler -I $HQ/build/include -L $HQ/build/lib -larmadillo -lhqcl $HQ/examples/tfd_q4_fullTempOptimization_sym_op.cpp
// $ ./tfd_q4_fullTempOptimization_sym_op
// $ cat tfd_minimization_results_sym_op.csv
// 0.001, -196.266, 102.798, -6.28315, 4.71243, -110169
// 0.00125893, 64.1739, -97.9999, -1.57086, -3.14159, -78353.8
// 0.00158489, 109.566, 36.7368, -1.57079, -4.71229, -55726.3
// 0.00199526, -79.7349, 218.9, 6.28324, 1.5708, -39633.3
// 0.00251189, -51.463, -104.6, -6.28311, 1.57077, -28187.7
// 0.00316228, -31.0372, -76.2076, 1.57087, -4.71243, -20047.5
// 0.00398107, 17.653, 118.633, 4.71222, 3.14156, -14258
// 0.00501187, -89.9379, 233.341, -17.2789, -10.9954, -10140.5
// 0.00630957, -2.74985, 126.471, -6.28299, -0.000176139, -7212.07
// 0.00794328, -38.8793, -60.4591, -4.71281, -9.42504, -5129.33
// 0.01, 60.8685, -112.28, 26.7031, -1.57033, -3648.05
// 0.0125893, 89.1426, 62.0545, -7.85478, -15.7073, -2594.55
// 0.0158489, 14.5299, -5.49131, 17.2797, 6.2824, -1845.28
// 0.0199526, -144.121, -153.149, -12.5649, -21.9923, -1312.4
// 0.0251189, 4.31967, 71.4735, -10.9977, 37.7009, -933.411
// 0.0316228, -39.663, 5.49986, 23.5589, -40.8383, -663.874
// 0.0398107, -35.7359, -112.308, -84.8272, 1.56742, -472.184
// 0.0501187, -82.0739, 40.0544, -84.817, 7.85878, -335.862
// 0.0630957, -125.271, 82.4673, 15.7163, 160.215, -238.924
// 0.0794328, 88.3574, -55.7618, -78.5515, -70.6952, -170.002
// 0.1, 28.6671, -22.7765, -73.8439, -0.013147, -121.015
// 0.125893, -12.1736, -58.9045, 72.2797, 15.6895, -86.2173
// 0.158489, 45.9457, -184.568, -31.4482, -26.7295, -61.5284
// 0.199526, 57.7269, -2.35618, -25.0879, -125.627, -44.0508
// 0.251189, 31.8091, -253.684, 89.597, 72.2065, -31.7292
// 0.316228, -100.924, 137.445, 15.7915, -1.50195, -23.1047
// 0.398107, -14.5299, -58.9047, 12.4567, -58.0273, -17.135
// 0.501187, -45.9457, -2.35612, -29.7067, 138.349, -13.0633
// 0.630957, -89.928, -35.3428, -6.45084, 116.387, -10.3273
// 0.794328, -130.769, 8.63943, 99.1549, 58.2943, -8.50777
// 1, 125.271, 51.0509, -112.878, 41.0378, -7.30152
// 1.25893, -18.4568, 60.4758, 161.551, 45.3397, -6.49974
// 1.58489, -11.3883, 131.161, -33.2492, 91.3304, -5.9641
// 1.99526, 134.696, -47.9092, 144.23, 9.65503, -5.60513
// 2.51189, -19.2422, 11.781, 105.546, -42.1784, -5.36419
// 3.16228, 125.271, 129.591, 134.767, -73.5934, -5.20239
// 3.98107, -96.9966, -96.6041, -99.2976, -30.0796, -5.09373
// 5.01187, 13.7444, -43.1968, 84.4721, -35.8933, -5.02035
// 6.30957, -220.304, -3.92696, -48.3331, -101.866, -4.97062
// 7.94328, 75.0056, -74.6127, 103.303, 42.6478, -4.93664
// 10, 29.4525, 71.4713, -103.296, -211.82, -4.91324
// 12.5893, -68.7223, -126.449, -44.3632, 116.002, -4.89702
// 15.8489, -71.8639, 168.861, -25.5169, 150.559, -4.88571
// 19.9526, 185.747, 19.6351, -9.81109, 101.863, -4.87779
// 25.1189, 59.2976, -16.4933, 83.6405, 124.332, -4.87221
// 31.6228, -5.89021, -101.316, 40.4511, 9.18611, -4.86827
// 39.8107, 117.417, 109.17, 105.634, 48.9337, -4.86549
// 50.1187, -168.468, -58.9047, -138.621, -26.4645, -4.86352
// 63.0957, -86.0011, -131.161, 45.1617, 262.084, -4.86212
// 79.4328, 111.134, 57.3341, -137.838, -87.7253, -4.86113
// 100, 88.3574, -167.29, 137.051, 246.376, -4.86042
// 125.893, 65.5808, 18.0642, -9.81727, 50.5047, -4.85992
// 158.489, 156.687, 102.887, -16.8864, -64.1633, -4.85957
// 199.526, -218.733, 73.0421, -50.6579, 20.6597, -4.85932
// 251.189, -146.477, -137.445, -35.7359, 28.5137, -4.85914
// 316.228, 59.2976, 2.3563, -214.021, 28.5137, -4.85901
// 398.107, 50.6582, -38.4844, 7.4613, 31.1768, -4.85892
// 501.187, 55.3706, 237.976, -4.31966, -33.226, -4.85885
// 630.957, 21.5985, -96.604, -15.3153, 59.9297, -4.85881
// 794.328, -61.6537, -104.458, -120.559, 92.9164, -4.85877
// 1000, 173.18, -69.9003, -39.6626, -85.0623, -4.85875
//

/// NOTE: Pick an include path depending on whether running in development
/// enviromnment or production environment
//===----------------------------------------------------------------------===//
/// Production mode
#include <clang/Quantum/quintrinsics.h>

// AUTHOR : Shavindra Premaratne (August 26th 2021)
// Update : Comments added throughout to clarify the intentions behind code
// blocks (08/27/2021)

/// Development mode
// #include "../../clang/include/clang/Quantum/quintrinsics.h"

/// Quantum Runtime Library APIs
#include <quantum.hpp>

#include <cassert> // to assert the size of the Probability Register
#include <fstream> // to write to files
#include <iostream>
#include <math.h> // to use pow()
#include <vector>

// to perform the minimization of TFD generation cost function

#include <armadillo>
#include <ensmallen.hpp>

#include "SymbolicOperator.hpp"
#include "SymbolicOperatorUtils.hpp"

using namespace hybrid::quantum::core;
using namespace iqsdk;

// ----------------------------------------------------------------------------------------
// ALL QUANTUM HELPER CODE
// ----------------------------------------------------------------------------------------

// Define the number of qubits needed for compilation
const int N = 4;
qbit QubitReg[N];
cbit CReg[N];

// Special global array to hold dynamic parameters for quantum algorithm
// alpha1, alpha2, gamma1, gamma2 is the order of the variational angles used in
// here
double QuantumVariableParams[N + N * 2];

// Optimization steps
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
                   SymbolicOperator &symbop, QWCMap &m_qwc_groups) {
  // start point in the global parameter array corresponding to the basis change
  int basis_change_variable_param_start_indx = 4;
  // total cost
  double total_cost = 0.0;

  QuantumVariableParams[0] = 2 * params[0];
  QuantumVariableParams[1] = 2 * params[1];
  QuantumVariableParams[2] = 2 * params[2];
  QuantumVariableParams[3] = 2 * params[3];

  for (auto &m_qwc_group : m_qwc_groups) {
    std::vector<double> ProbReg;

    std::vector<double> variable_params;
    variable_params.reserve(N * 2);

    SymbolicOperatorUtils::applyBasisChange(m_qwc_group.second, variable_params,
                                            N);

    std::vector<std::reference_wrapper<qbit>> qids;
    for (int qubit = 0; qubit < N; ++qubit) {
      qids.push_back(std::ref(QubitReg[qubit]));
    }

    QuantumVariableParams[4] = 0;
    QuantumVariableParams[5] = 0;
    QuantumVariableParams[6] = 0;
    QuantumVariableParams[7] = 0;
    QuantumVariableParams[8] = 0;
    QuantumVariableParams[9] = 0;
    QuantumVariableParams[10] = 0;
    QuantumVariableParams[11] = 0;

    for (auto indx = 0; indx < variable_params.size(); ++indx) {
      QuantumVariableParams[basis_change_variable_param_start_indx + indx] =
          variable_params[indx];
    }

    // performing the experiment, and storing the data in ProbReg
    tfdQ4();

    ProbReg = iqs_device.getProbabilities(qids);

    double current_pstr_val = SymbolicOperatorUtils::getExpectValSetOfPaulis(
        symbop, m_qwc_group.second, ProbReg, N);
    total_cost += current_pstr_val;
  }

  return total_cost;
}

// Constructing a lambda function to be used for a single optimization iteration
// This function is directly called by the dlib optimization routine
double ansatz_run_lambda(SymbolicOperator &symbop, QWCMap &m_qwc_groups,
                         const arma::mat &params, double inv_temp,
                         FullStateSimulator &iqs_device) {
  // loading the new variational angles into the special global array for hybrid
  // compilation
  QuantumVariableParams[0] = params(0);
  QuantumVariableParams[1] = params(1);
  QuantumVariableParams[2] = params(2);
  QuantumVariableParams[3] = params(3);

  // runs the kernel to compute the total cost
  double total_cost = run_qkernel(iqs_device, params, symbop, m_qwc_groups);

  return total_cost;
}

///
///@brief Object class to apply optimization step
///
class EnergyOfAnsatz {
public:
  EnergyOfAnsatz(SymbolicOperator &symbop_, QWCMap &_m_qwc_groups,
                 arma::mat &params_, double beta_,
                 FullStateSimulator &_iqs_device)
      : symbop(symbop_), m_qwc_groups(_m_qwc_groups), params(params_),
        beta(beta_), iqs_device(_iqs_device) {}

  double Evaluate(const arma::mat &theta) {
    steps_count++;

    double total_cost;
    total_cost =
        ansatz_run_lambda(symbop, m_qwc_groups, params, beta, iqs_device);

    return total_cost;
  }

private:
  SymbolicOperator symbop;
  const arma::mat &params;
  double beta;
  QWCMap &m_qwc_groups;
  FullStateSimulator &iqs_device;
};

int main() {
  // H = X3 + X2 + X1 + X0 + 1.60 (Z3Z2 + Z1Z0) - (beta ^ -1.48((X1X3 +
  // X0X2)+(Z1Z3 + Z0Z2)))

  // resetting steps count
  steps_count = 0;

  // Clearing the file contents used to store all the results from TFD
  // minimization
  std::ofstream file_to_clear;
  file_to_clear.open("tfd_minimization_results_sym_op.csv",
                     std::ofstream::out | std::ofstream::trunc);
  file_to_clear.close();

  bool print_once_flag = false;
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

    std::cout << "Hamiltonian: " << symbop.getCharString() << "\n";

    // QWC
    QWCMap m_qwc_groups =
        SymbolicOperatorUtils::getQubitwiseCommutationGroups(symbop);
    if (!print_once_flag) {
      std::cout << "Number of Qubitwise Commutation (QWC) Groups : "
                << m_qwc_groups.size() << "\n";
      std::cout << "Qubitwise Commutation (QWC) Groups: " << m_qwc_groups
                << "\n";
      print_once_flag = true;
    }

    // Using SA algorithm
    // Starting parameters
    arma::mat params_sa(std::vector<double>{0, 0, 0, 0});

    // Setup quantum device
    IqsConfig iqs_config(/*num_qubits*/ N,
                          /*simulation_type*/ "noiseless");
    FullStateSimulator iqs_device(iqs_config);
    if (QRT_ERROR_SUCCESS != iqs_device.ready()) {
      return -1;
    }

    // arbitrary function
    EnergyOfAnsatz eoa_sa(symbop, m_qwc_groups, params_sa, beta, iqs_device);

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
