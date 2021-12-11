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

#include "SymbolicOperatorUtils.hpp"

#include <armadillo>

#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
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

using namespace std;
using namespace arma;

const double FP_PI = 3.14159265359;
const double FP_2PI = 6.28318530718;
const double FP_PIby2 = 1.57079632679489661923;

// Two Pauli strings
// If for *every* qubit id, either (at least one local op is I) or (both local
// ops are same Pauli), then the two pauli strings are QWC.

// X0 Y1 Z2 and I0 Y1 Z2 --> yes
// X0 Y1 Z2 and X0 --> yes
// X0 Y1 Z2 and Z0 I1 I2 --> no

// X0 Y1 and Y0 X1 --> *not* QWC, though they *do* commute.

// X0*Y0 = iZ0
// Y1*X1 = -iZ1
// (X0 Y1)*(Y0 X1) = -(-1) Z0 Z1

Vec2DMat
SymbolicOperatorUtils::qubitwiseCommutation(const SymbolicOperator &symbop) {

  int M = symbop.op_sum.size();

  Vec2DMat vmat = Vec2DMat(M, vector<int>(M, 0));
  int row = 0, col = 0;
  for (auto iIter = symbop.op_sum.begin(); iIter != symbop.op_sum.end();
       ++iIter, ++row) {
    col = row + 1;

    for (auto jIter = next(iIter); jIter != symbop.op_sum.end();
         ++jIter, ++col) {

      pstring::const_iterator loc_i_end = iIter->first.end();
      pstring::const_iterator loc_j_end = jIter->first.end();

      // Loop over the two individual Pauli strings
      for (pstring::const_iterator ipstrIter = iIter->first.begin();
           ipstrIter != loc_i_end; ++ipstrIter) {

        bool flag_non_commuting_found = false;

        int a = ipstrIter->first;
        char U = ipstrIter->second;

        for (pstring::const_iterator jpstrIter = jIter->first.begin();
             jpstrIter != loc_j_end; ++jpstrIter) {

          int b = jpstrIter->first;
          char V = jpstrIter->second;

          if (U != V && a == b) {
            // Only one mismatch is requried for the two pstrings to *not* be
            // QWC Add edge for non-commutation
            vmat[row][col] = vmat[col][row] = 1;
            flag_non_commuting_found = true;
            break;
          }
        }

        if (flag_non_commuting_found)
          break;
      }
    }
  }

  return vmat;
}

string
SymbolicOperatorUtils::getCharString_pstring(const pstring &inp_pstring) {

  string charstring = "";

  auto it = inp_pstring.begin();
  while (it != inp_pstring.end()) {
    charstring += it->second;
    charstring += to_string(it->first);
    charstring += " ";

    it++;
  }

  return charstring;
}

double
SymbolicOperatorUtils::getExpectValSglPauli(const pstring &pstr,
                                            const std::vector<double> ProbReg,
                                            int num_qbits, double eps) {
  arma::Row<double> I = {1, 1};
  arma::Row<double> Z = {1, -1};
  arma::Row<double> vR = {1};
  double exp_val;

  arma::rowvec rv(ProbReg);

  for (auto i = 0; i < num_qbits; ++i) {
    auto foundX = pstr.find(op_pair(i, 'X'));
    auto foundY = pstr.find(op_pair(i, 'Y'));
    auto foundZ = pstr.find(op_pair(i, 'Z'));

    if (foundX != pstr.end() || foundY != pstr.end() || foundZ != pstr.end()) {
      vR = arma::kron(vR, Z).as_row();
    } else {
      vR = arma::kron(vR, I).as_row();
    }
  }

  exp_val = arma::sum((rv % vR).as_row());

  return exp_val;
}

double SymbolicOperatorUtils::getExpectValSetOfPaulis(
    const std::vector<pstring> &v_pstr, const std::vector<double> ProbReg,
    int num_qbits, double eps) {
  arma::Row<double> I = {1, 1};
  arma::Row<double> Z = {1, -1};
  arma::Row<double> vR = {1};
  arma::rowvec rv(ProbReg);
  double exp_val;

  for (const auto &pstr : v_pstr) {
    exp_val += getExpectValSglPauli(pstr, ProbReg, num_qbits);
  }

  return exp_val;
}

double SymbolicOperatorUtils::getExpectVal(SymbolicOperator &symbop,
                                           const std::vector<double> ProbReg,
                                           int num_qbits, double eps,
                                           METHOD method) {
  arma::Row<double> I = {1, 1};
  arma::Row<double> Z = {1, -1};
  arma::Row<double> vR = {1};
  arma::rowvec rv(ProbReg);
  double exp_val;

  for (const auto &pstr : symbop.getOrderedPStringList()) {
    exp_val += symbop.op_sum[pstr].real() *
               getExpectValSglPauli(pstr, ProbReg, num_qbits);
  }

  return exp_val;
}

char SymbolicOperatorUtils::findFirstPauliStringBasis(const pstring &pstr) {
  char first_pstr_basis;
  // Find first basis change
  for (const auto &ps : pstr) {
    if (ps.second == 'X' || ps.second == 'Y' || ps.second == 'Z') {
      first_pstr_basis = ps.second;
      break;
    }
  }
  return first_pstr_basis;
}

void SymbolicOperatorUtils::applyBasisChange(
    const pstring &pstr, std::vector<double> &variable_params, int num_qbits) {

  char first_pstr_basis = findFirstPauliStringBasis(pstr);
  for (auto i = 0; i < num_qbits; ++i) {
    auto foundX = pstr.find(op_pair(i, 'X'));
    auto foundY = pstr.find(op_pair(i, 'Y'));
    auto foundZ = pstr.find(op_pair(i, 'Z'));

    if (foundX != pstr.end()) {
      variable_params.push_back(FP_PIby2);
      variable_params.push_back(FP_PI);
    } else if (foundY != pstr.end()) {
      variable_params.push_back(FP_PI);
      variable_params.push_back(FP_PIby2);
    } else if (foundZ != pstr.end()) {
      variable_params.push_back(0);
      variable_params.push_back(0);
    } else {
      if (first_pstr_basis == 'X') {
        variable_params.push_back(FP_PIby2);
        variable_params.push_back(FP_PI);
      } else if (first_pstr_basis == 'Y') {
        variable_params.push_back(FP_PI);
        variable_params.push_back(FP_PIby2);
      } else if (first_pstr_basis == 'Z') {
        variable_params.push_back(0);
        variable_params.push_back(0);
      }
    }
  }
}

} // namespace core
} // namespace quantum
} // namespace hybrid
