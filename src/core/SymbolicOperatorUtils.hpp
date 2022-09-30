//===----------------------------------------------------------------------===//
// INTEL CONFIDENTIAL
//
// Copyright 2021-2022 Intel Corporation.
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

#ifndef __HYBRID_QUANTUM_CORE_SYMBOLIC_OPERATOR_UTILS_DOT_HPP__
#define __HYBRID_QUANTUM_CORE_SYMBOLIC_OPERATOR_UTILS_DOT_HPP__

#include "SymbolicOperator.hpp"

#include <complex>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace hybrid {
namespace quantum {
namespace core {

enum class METHOD { DEFAULT = 0, QWC };

const double FP_PI = 3.14159265359;
const double FP_2PI = 6.28318530718;
const double FP_PIby2 = 1.57079632679489661923;

///
/// @brief Utility methods for Symbolic Operator Class
///
class SymbolicOperatorUtils {

public:
  ///
  /// @brief Get graph describing qubit-wise commutation
  ///
  /// @param symbop - SymbolicOperator object
  /// @return Vec2DMat
  ///
  static Vec2DMat qubitwiseCommutation(const SymbolicOperator &symbop);

  ///
  /// @brief Get groups that qubit-wise commute
  ///
  /// @param qwcgraph - Vec2DMat object
  /// @return vector<int>
  ///
  static std::vector<int> getGroupsQWC(const Vec2DMat &qwcgraph);

  ///
  /// @brief Get the string version of pstring object
  ///
  /// @param pstr - Pauli string
  /// @return std::string
  ///
  static std::string getCharString_pstring(const pstring &pstr);

  ///
  /// @brief Get the expectation value of the symbolic operator object
  ///
  /// @param symbop - SymbolicOperator object
  /// @param ProbReg - Probabilities obtained from the state vector amplitudes
  /// @param num_qubits - Number of qubits
  /// @param eps - Desired precision (epsilon)
  /// @param method - Method based on Qubit properties
  /// @return double
  ///
  static double getExpectVal(SymbolicOperator &symbop,
                             const std::vector<double> ProbReg, int num_qubits,
                             double eps = 0.0, METHOD method = METHOD::DEFAULT);

  ///
  /// @brief Get the expectation value for a set of pauli strings/objects
  ///
  /// @param pstr - Set of pauli operators
  /// @param ProbReg - Probabilities obtained from the state vector amplitudes
  /// @param num_qubits - Number of qubits
  /// @param eps - Desired precision (epsilon)
  /// @return double
  ///
  static double getExpectValSetOfPaulis(const std::vector<pstring> &v_pstr,
                                        const std::vector<double> ProbReg,
                                        int num_qubits, double eps = 0.0);

  ///
  /// @brief Get the expectation value for a set of pauli strings/objects
  ///
  /// @param symbop - Symbolic representation
  /// @param pstr - Set of pauli operators
  /// @param ProbReg - Probabilities obtained from the state vector amplitudes
  /// @param num_qubits - Number of qubits
  /// @param eps - Desired precision (epsilon)
  /// @return double
  ///
  static double getExpectValSetOfPaulis(SymbolicOperator &symbop,
                                        std::set<pstring> &s_pstr,
                                        const std::vector<double> ProbReg,
                                        int num_qubits, double eps = 0.0);

  ///
  /// @brief Get the expectation value of single pauli string/object
  ///
  /// @param pstr - Pauli String
  /// @param ProbReg - Probabilities obtained from the state vector amplitudes
  /// @param num_qubits - Number of qubits
  /// @param eps - Desired precision (epsilon)
  /// @return double
  ///
  static double getExpectValSglPauli(const pstring &pstr,
                                     const std::vector<double> ProbReg,
                                     int num_qubits, double eps = 0.0);

  /// @brief Applies basis change to a single Pauli string. For the missing
  /// qubit in the pauli term, first encountered
  ///        pauli operator's basis change is applied
  /// @param pstr - Pauli string
  /// @param variable_params - Variable Parameters
  /// @param num_qbits - Number of qubits
  static void applyBasisChange(const pstring &pstr,
                               std::vector<double> &variable_params,
                               int num_qbits);

  /// @brief Applies basis change to a set of pauli strings/objects.
  /// @param s_pstr - Set of pauli strings/objects
  /// @param variable_params - Variable Parameters
  /// @param num_qbits - Number of qubits
  static void applyBasisChange(std::set<pstring> &s_pstr,
                               std::vector<double> &variable_params,
                               int num_qbits, bool qwc_check = false);

  static QWCMap getQubitwiseCommutationGroups(SymbolicOperator &symbop,
                                              int num_qbits);

  /// @brief Finds the first pauli operator in the pauli string
  /// @param pstr - Pauli string
  /// @return
  static char findFirstPauliStringBasis(const pstring &pstr);

}; // end of class SymbolicOperator

std::ostream &operator<<(std::ostream &s, const QWCMap &qwc_groups_mapping);

} // namespace core
} // namespace quantum
} // namespace hybrid

#endif // __HYBRID_QUANTUM_CORE_SYMBOLIC_OPERATOR_UTILS_DOT_HPP__
