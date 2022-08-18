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
#include <set> // These are ordered
#include <string>
#include <utility> // pair
#include <vector>

namespace hybrid {
namespace quantum {
namespace core {

using op_pair = std::pair<int, char>;
using pstring = std::set<op_pair>;
using ComplexDP = std::complex<double>;
using MapPString = std::map<pstring, ComplexDP>;
using Vec2DMat = std::vector<std::vector<int>>;

enum class METHOD { DEFAULT = 0, QWC };

const double FP_PI = 3.14159265359;
const double FP_2PI = 6.28318530718;
const double FP_PIby2 = 1.57079632679489661923;

///
/// @brief Utility methods for Symbolic Operator Class
///
class SymbolicOperatorUtils {
  // Eventually need:
  // -- Commutation graph
  // -- Qubit-wise commutation graph
  // -- Anti-commutation graph
  // TODO:
  //- get the string addTerm to work. (do *not* do multi-term)
  //- get scalar mult to work

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
  /// @brief Get the Expect Val object
  ///
  /// @param symbop - SymbolicOperator object
  /// @param ProbReg - Probabilities obtained from the state vector amplitudes
  /// @param num_qubits - Number of qubits
  /// @param eps - Desired precision (epsilon)
  /// @param method - Method based on Qubit properties
  /// @return double
  ///
  static double
  getExpectVal(SymbolicOperator &symbop,
               const std::vector<double> ProbReg,
               int num_qubits, double eps = 0.0, METHOD method = METHOD::DEFAULT);

  ///
  /// @brief Get the Expect Val Set Of Paulis object
  ///
  /// @param pstr - List of Pauli Strings
  /// @param ProbReg - Probabilities obtained from the state vector amplitudes
  /// @param num_qubits - Number of qubits
  /// @param eps - Desired precision (epsilon)
  /// @return double
  ///
  static double getExpectValSetOfPaulis(
      const std::vector<pstring> &pstr,
      const std::vector<double> ProbReg,
      int num_qubits, double eps = 0.0);

  ///
  /// @brief Get the Expect Val Sgl Pauli object
  ///
  /// @param pstr - Pauli String
  /// @param ProbReg - Probabilities obtained from the state vector amplitudes
  /// @param num_qubits - Number of qubits
  /// @param eps - Desired precision (epsilon)
  /// @return double
  ///
  static double getExpectValSglPauli(
      const pstring &pstr,
      const std::vector<double> ProbReg,
      int num_qubits, double eps = 0.0);

  ///
  /// @brief Applies the basis change to the given Pauli string
  ///
  /// @param pstr - Pauli String
  /// @param variable_params - Probabilities obtained from the state vector amplitudes
  /// @param num_qubits - Desired precision (epsilon)
  ///
  static void applyBasisChange(const pstring &pstr,
                              std::vector<double> &variable_params, int
                              num_qubits);

  ///
  /// @brief Finds the basis of the first Pauli string
  ///
  /// @param pstr - Pauli String
  /// @return char
  ///
  static char findFirstPauliStringBasis(const pstring &pstr);

}; // end of class SymbolicOperator

} // namespace core
} // namespace quantum
} // namespace hybrid

#endif // __HYBRID_QUANTUM_CORE_SYMBOLIC_OPERATOR_UTILS_DOT_HPP__

