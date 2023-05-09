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
#include <clang/Quantum/quintrinsics.h>

// AUTHOR : Shavindra Premaratne (August 26th 2021)
// Update : Comments added throughout to clarify the intentions behind code
// blocks (08/27/2021)

/// Development mode
// #include "../../clang/include/clang/Quantum/quintrinsics.h"

#include <quantum.hpp>

#include <cassert> // to assert the size of the Probability Register
#include <fstream> // to write to files
#include <iostream>
#include <math.h> // to use pow()
#include <vector>

// to perform the minimization of TFD generation cost funtion
#include <dlib/global_optimization.h>
#include <dlib/optimization.h>

#include "SymbolicOperator.hpp"
#include "SymbolicOperatorUtils.hpp"
#include <armadillo>

using namespace iqsdk;
using namespace hybrid::quantum::core;

// ----------------------------------------------------------------------------------------
// ALL QUANTUM HELPER CODE
// ----------------------------------------------------------------------------------------

// Define the number of qubits needed for compilation
const int N_sub = 8; // Number of qubits in subsystem (thermal state size)
const int N_ss = 2;  // Number of subsystems (CAUTION : this is not a general
                     // parameter to be changed)
const int N = N_ss * N_sub; // Total number of qubits (TFD state size)

const int N_var_angles = 4;
const int N_map_angles = 2 * N;

qbit QubitReg[N];
cbit CReg[N];

/* Array to hold dynamic parameters for quantum algorithm */
double QuantumVariableParams[N_var_angles];
double QuantumMappingParams[N_map_angles];

/*
When using the Intel Quantum Simulator (IQS) at the lowest level, it is
necessary to decide how to calculate the cost function based on the results
returned. If using the full state vector, then only one experiment is required
per iteration. Here, we will be using the probabilities corresponding to each of
the basis states (e.g. 0000, 1010) These are returned simply by calculating the
probability corresponding to each of the state vector element and stored in the
special array ProbabilityRegister.
*/

quantum_kernel void TFD_terms() {

  int shifting_int =
      0; // this is to shift the variational parameter in case of more steps

  int index_intraX = 0;
  int index_intraZ = 0;
  int index_interX = 0;
  int index_interZ = 0;

  // ============================================================================================================
  // Single qubit variational terms
  for (index_intraX = 0; index_intraX < N; index_intraX++)
    RX(QubitReg[index_intraX], QuantumVariableParams[0 + shifting_int]);
  // ============================================================================================================

  // ============================================================================================================
  // Two-qubit intra-system variational terms (adjacent)
  for (int grand_intraZ = 0; grand_intraZ < N_sub - 1; grand_intraZ++) {
    for (index_intraZ = 0; index_intraZ < N_ss; index_intraZ++)
      CNOT(QubitReg[grand_intraZ + N_sub * index_intraZ + 1],
           QubitReg[grand_intraZ + N_sub * index_intraZ]);
    for (index_intraZ = 0; index_intraZ < N_ss; index_intraZ++)
      RZ(QubitReg[grand_intraZ + N_sub * index_intraZ],
         QuantumVariableParams[1 + shifting_int]);
    for (index_intraZ = 0; index_intraZ < N_ss; index_intraZ++)
      CNOT(QubitReg[grand_intraZ + N_sub * index_intraZ + 1],
           QubitReg[grand_intraZ + N_sub * index_intraZ]);
  }

  if (N_sub > 2) {
    // Two-qubit intra-system variational terms (boundary term)
    for (index_intraZ = 0; index_intraZ < N_ss; index_intraZ++)
      CNOT(QubitReg[N_sub * index_intraZ],
           QubitReg[N_sub * index_intraZ + (N_sub - 1)]);
    for (index_intraZ = 0; index_intraZ < N_ss; index_intraZ++)
      RZ(QubitReg[N_sub * index_intraZ + (N_sub - 1)],
         QuantumVariableParams[1 + shifting_int]);
    for (index_intraZ = 0; index_intraZ < N_ss; index_intraZ++)
      CNOT(QubitReg[N_sub * index_intraZ],
           QubitReg[N_sub * index_intraZ + (N_sub - 1)]);
  }
  // ============================================================================================================

  // ============================================================================================================
  // two-qubit inter-system XX variational terms
  for (index_interX = 0; index_interX < N_sub; index_interX++) {
    RY(QubitReg[index_interX + N_sub], -1.57079632679);
    RY(QubitReg[index_interX], -1.57079632679);
  }
  for (index_interX = 0; index_interX < N_sub; index_interX++)
    CNOT(QubitReg[index_interX + N_sub], QubitReg[index_interX]);
  for (index_interX = 0; index_interX < N_sub; index_interX++)
    RZ(QubitReg[index_interX], QuantumVariableParams[2 + shifting_int]);
  for (index_interX = 0; index_interX < N_sub; index_interX++)
    CNOT(QubitReg[index_interX + N_sub], QubitReg[index_interX]);
  for (index_interX = 0; index_interX < N_sub; index_interX++) {
    RY(QubitReg[index_interX + N_sub], 1.57079632679);
    RY(QubitReg[index_interX], 1.57079632679);
  }
  // ============================================================================================================

  // ============================================================================================================
  // two-qubit inter-system ZZ variational terms
  for (index_interZ = 0; index_interZ < N_sub; index_interZ++)
    CNOT(QubitReg[index_interZ], QubitReg[index_interZ + N_sub]);
  for (index_interZ = 0; index_interZ < N_sub; index_interZ++)
    RZ(QubitReg[index_interZ + N_sub], QuantumVariableParams[3 + shifting_int]);
  for (index_interZ = 0; index_interZ < N_sub; index_interZ++)
    CNOT(QubitReg[index_interZ], QubitReg[index_interZ + N_sub]);
  // ============================================================================================================
}

quantum_kernel void PrepZAll() {
  // Initialization of the qubits
  for (int Index = 0; Index < N; Index++)
    PrepZ(QubitReg[Index]);
}

quantum_kernel void BellPrep() {
  // Index to loop over later
  int Index = 0;

  // preparation of Bell pairs (T -> Infinity)
  for (Index = 0; Index < N_sub; Index++)
    RY(QubitReg[Index], 1.57079632679);
  for (Index = 0; Index < N_sub; Index++)
    CNOT(QubitReg[Index], QubitReg[Index + N_sub]);
}

quantum_kernel void DynamicMapping() {
  // not part of ansatz but additional layers to help
  // with X to Z or Y to Z basis mapping
  for (int qubit_index = 0; qubit_index < N; qubit_index++) {
    int map_index = 2 * qubit_index;
    RY(QubitReg[qubit_index], QuantumMappingParams[map_index]);
    RX(QubitReg[qubit_index], QuantumMappingParams[map_index + 1]);
  }
}

/*
Experiment using HQCL and dynamically setting up the measurement basis changes
*/
quantum_kernel void TFD_full() {
  PrepZAll();
  BellPrep();
  TFD_terms();
  DynamicMapping();
  // MeasZAll();
}

// Description below from official dlib documentation:
// In dlib, most of the general purpose solvers optimize functions that take a
// column vector as input and return a double.  So here we make a typedef for a
// variable length column vector of doubles.  This is the type we will use to
// represent the input to our objective functions which we will be minimizing.
// typedef dlib::matrix<double,0,1> column_vector;
typedef dlib::matrix<double, N_var_angles, 1> column_vector;

void clear_results_file(const std::string &csv_output_name);
void write_to_results_file(const std::string &csv_output_name, double beta,
                           column_vector starting_point, double result);

int main() {
  /// Setup quantum device
  iqsdk::IqsConfig iqs_config(N, "noiseless");
  iqsdk::FullStateSimulator iqs_device(iqs_config);
  if (iqsdk::QRT_ERROR_SUCCESS != iqs_device.ready()) {
    return 1;
  }

  const std::string &csv_output_name =
      "symbOp_20220902_TFD" + std::to_string(N) + "x1_O1.csv";

  // Clearing the file contents used to store all the results from TFD
  // minimization
  clear_results_file(csv_output_name);

  // initial starting point. Defining it here means I will reuse the best result
  // from from previous temperature when starting the next temperature run
  column_vector starting_point = {0, 0, 0, 0};

  // Looping over the index for range of inverse temperatures
  for (int beta_index = -30; beta_index < 31; beta_index += 1) {
    // calculating the actual inverse temperature that is used during
    // calculations
    double beta = pow(10, (double)beta_index / 10);

    // construct object
    SymbolicOperator symbop;
    pstring sym_term;

    // Single qubit variational terms
    for (int index_intraX = 0; index_intraX < N; index_intraX++) {
      sym_term = {std::make_pair(index_intraX, 'X')};
      symbop.addTerm(sym_term, 1);
    }

    // Two-qubit intra-system variational terms (adjacent)
    for (int grand_intraZ = 0; grand_intraZ < N_sub - 1; grand_intraZ++) {
      for (int index_intraZ = 0; index_intraZ < N_ss; index_intraZ++) {
        int qIndex0 = grand_intraZ + N_sub * index_intraZ;
        int qIndex1 = grand_intraZ + N_sub * index_intraZ + 1;
        // std::cout << "(" << qIndex0 << ", " << qIndex1 << ")\n";

        sym_term = {std::make_pair(qIndex0, 'Z'), std::make_pair(qIndex1, 'Z')};
        symbop.addTerm(sym_term, 1.00);
      }
    }
    // Two-qubit intra-system variational terms (boundary term)
    if (N_sub > 2) {
      for (int index_intraZ = 0; index_intraZ < N_ss; index_intraZ++) {
        int qIndex0 = N_sub * index_intraZ;
        int qIndex1 = N_sub * index_intraZ + (N_sub - 1);
        // std::cout << "(" << qIndex0 << ", " << qIndex1 << ")\n";

        sym_term = {std::make_pair(qIndex0, 'Z'), std::make_pair(qIndex1, 'Z')};
        symbop.addTerm(sym_term, 1.00);
      }
    }

    // two-qubit inter-system XX variational terms
    for (int index_interX = 0; index_interX < N_sub; index_interX++) {
      int qIndex0 = index_interX;
      int qIndex1 = index_interX + N_sub;
      // std::cout << "(" << qIndex0 << ", " << qIndex1 << ")\n";

      sym_term = {std::make_pair(qIndex0, 'X'), std::make_pair(qIndex1, 'X')};
      symbop.addTerm(sym_term, -pow(beta, -1.00));
    }

    // two-qubit inter-system XX variational terms
    for (int index_interZ = 0; index_interZ < N_sub; index_interZ++) {
      int qIndex0 = index_interZ;
      int qIndex1 = index_interZ + N_sub;

      sym_term = {std::make_pair(qIndex0, 'Z'), std::make_pair(qIndex1, 'Z')};
      symbop.addTerm(sym_term, -pow(beta, -1.00));
    }

    // Constructing a lambda function to be used for a single optimization
    // iteration This function is directly called by the dlib optimization
    // routine
    auto ansatz_run_lambda = [&](const column_vector &params) {
      // Setting all the variational angles to input values.
      for (int qubit_index = 0; qubit_index < N_map_angles; qubit_index++)
        QuantumVariableParams[qubit_index] = params(qubit_index);

      double total_cost = 0.0;

      // Setting all the variational angles to input values.
      for (int qubit_index = 0; qubit_index < N_map_angles; qubit_index++)
        QuantumVariableParams[qubit_index] = params(qubit_index);

      // loop through all the Pauli strings present in the SymbolicOperator
      // object
      for (const auto &pstr : symbop.getOrderedPStringList()) {
        // Setting all the mapping angles to 0.
        for (int map_index = 0; map_index < N_map_angles; map_index++)
          QuantumMappingParams[map_index] = 0;

        // used to enable mapping of X-basis to Z-basis or Y-basis to Z-basis
        std::vector<double> variable_params;
        variable_params.reserve(N_map_angles);
        SymbolicOperatorUtils::applyBasisChange(pstr, variable_params, N);

        for (auto indx = 0; indx < variable_params.size(); ++indx)
          QuantumMappingParams[indx] = variable_params[indx];

        // performing the experiment, and storing the data in ProbReg
        std::vector<std::reference_wrapper<qbit>> qids;
        for (int qubit = 0; qubit < N; ++qubit)
          qids.push_back(std::ref(QubitReg[qubit]));

        TFD_full();
        std::vector<double> ProbReg = iqs_device.getProbabilities(qids);

        // calculate the expectation value
        double current_pstr_val =
            symbop.op_sum[pstr].real() *
            SymbolicOperatorUtils::getExpectValSglPauli(pstr, ProbReg, N);

        total_cost += current_pstr_val;
      }

      return total_cost;
    };

    // running the full optimization for a given temperature
    auto result = dlib::find_min_bobyqa(
        ansatz_run_lambda, starting_point,
        2 * N_var_angles + 1, // number of interpolation points
        dlib::uniform_matrix<double>(N_var_angles, 1,
                                     -7.0), // lower bound constraint
        dlib::uniform_matrix<double>(N_var_angles, 1,
                                     7.0), // upper bound constraint
        1.5,                               // initial trust region radius
        1e-5,                              // stopping trust region radius
        10000 // max number of objective function evaluations
    );

    // writing the results from a single temperature to file
    write_to_results_file(csv_output_name, beta, starting_point, result);
  }

  return 0;
}

void clear_results_file(const std::string &csv_output_name) {
  // Clearing the file contents used to store all the results from TFD
  // minimization
  std::ofstream file_to_clear;
  file_to_clear.open(csv_output_name,
                     std::ofstream::out | std::ofstream::trunc);
  file_to_clear.close();
}

void write_to_results_file(const std::string &csv_output_name, double beta,
                           column_vector starting_point, double result) {
  // writing the results from a single temperature to file
  std::ofstream results_file;
  results_file.open(csv_output_name, std::ios_base::app);
  results_file << beta << ", ";
  for (int i = 0; i < N_var_angles; i++)
    results_file << starting_point(i) << ", ";
  results_file << result << "\n";
  results_file.close();
}