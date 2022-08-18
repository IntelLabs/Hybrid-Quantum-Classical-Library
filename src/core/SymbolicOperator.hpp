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

#ifndef __HYBRID_QUANTUM_CORE_SYMBOLIC_OPERATOR_DOT_HPP__
#define __HYBRID_QUANTUM_CORE_SYMBOLIC_OPERATOR_DOT_HPP__

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
// graph color -> set of commuting terms
using QWCMap = std::map<int, std::set<pstring>>;

///
/// @brief Symbolic Operator Class
///
class SymbolicOperator {

public:
  ///
  ///@brief Add pauli string (as set of pairs)
  ///
  ///@param inpp - Pauli string
  ///@param k - Arbitrary coefficient
  ///@param check_validity - For this constructor, we do not allow >1 op on same
  /// qubit id
  ///
  void addTerm(pstring &inpp, ComplexDP k = 1, bool check_validity = true);

  ///
  ///@brief Add pauli string (as character string)
  /// Constructor with vector of strings, e.g. {"X10","Z4"}
  /// This allows for inputs like [X0 X0 Y2], which addTerm(pstring&) does not
  /// allow for. Empty vector string mean identity is being added
  ///
  ///@param vecstr - List of Pauli strings
  ///@param k - Arbitrary coefficient
  ///
  void addTerm(std::vector<std::string> &vecstr, ComplexDP k = 1);

  ///
  ///@brief Add identity term with arbitrary coefficient
  ///
  ///@param k - Arbitrary coefficient
  ///
  void addIdentTerm(ComplexDP k = 1);

  ///
  ///@brief Process local character string, i.e. "Z10"
  ///
  ///@param inp - Local character string
  ///@return op_pair
  ///
  op_pair processLocCharString(std::string inp);

  ///
  ///@brief Get the Char String object
  ///
  ///@return std::string
  ///
  std::string getCharString();

  ///
  ///@brief Get the Num Terms object
  ///
  ///@return int
  ///
  int getNumTerms();

  ///
  ///@brief Get the Ordered P String List object
  ///
  ///@return std::vector<pstring>
  ///
  std::vector<pstring> getOrderedPStringList();

  ///
  ///@brief Remove all terms
  ///
  ///
  void removeAllTerms();

  ///
  ///@brief Construct Symbolic Operator object from a file
  ///
  ///
  int construct_hamiltonian_from_file(std::string filename);

  ///
  ///@brief Construct a new Symbolic Operator object
  ///
  ///
  SymbolicOperator() = default;

  ///
  ///@brief Construct a new Symbolic Operator object
  ///
  ///
  SymbolicOperator(const SymbolicOperator &) = default;

  ///
  ///@brief Overloading of assignment operator
  ///
  ///@param inp - SymbolicOperator object
  ///@return SymbolicOperator&
  ///
  SymbolicOperator &operator=(const SymbolicOperator &inp) {
    this->op_sum = inp.op_sum;
    return *this;
  }

  ///
  ///@brief Overloaded '+' to add two SymbolicOperator objects
  ///
  ///@param inp - SymbolicOperator object
  ///@return SymbolicOperator
  ///
  SymbolicOperator operator+(const SymbolicOperator &inp);

  ///
  ///@brief Overloaded '*' to multiply two SymbolicOperator objects
  ///
  ///@param p_right - SymbolicOperator
  ///@return SymbolicOperator
  ///
  SymbolicOperator operator*(const SymbolicOperator &p_right);

#ifndef DOXYGEN_SKIP
  ///@brief Overloaded '*' to multiply SymbolicOperator with arbitrary
  /// coefficient
  /// TODO
  ///@param inp_right
  ///@return SymbolicOperator
  ///
  // SymbolicOperator operator*(const ComplexDP &inp_right);
#endif /* DOXYGEN_SKIP */

  ///
  ///@brief Check for equality between two SymbolicOperator objects
  ///
  ///@param rhs - SymbolicOperator
  ///@return true
  ///@return false
  ///
  bool operator==(const SymbolicOperator &rhs);

  /// @brief Strips whitespace from LHS
  /// @param s
  /// @param matches
  /// @return stripped whitespace string from LHS
  std::string lstrip(std::string s, std::string matches);

  /// @brief Strips whitespace from RHS
  /// @param s
  /// @param matches
  /// @return stripped whitespace string from RHS
  std::string rstrip(std::string s, std::string matches);

  /// @brief Strips whitespace throughout the string
  /// @param s
  /// @return string without the whitespace
  std::string stripws(std::string s);

  /// @brief Splits the string into tokens as per the given delimiter
  /// @param str
  /// @param delim
  /// @return vector of string tokens
  std::vector<std::string> split(const std::string &str, const char delimiter);

  ///
  /// @brief Zero Threshold Value
  ///
  double zero_thresh = 1.0e-11;

  ///
  /// @brief Stores the pauli term and it's coefficient
  ///
  MapPString op_sum; // The data structure

  ///
  /// @brief Stores qubit-wise commutation graph
  ///
  Vec2DMat adj_matrix;
}; // end of class SymbolicOperator

} // namespace core
} // namespace quantum
} // namespace hybrid
#endif // __HYBRID_QUANTUM_CORE_SYMBOLIC_OPERATOR_DOT_HPP__
