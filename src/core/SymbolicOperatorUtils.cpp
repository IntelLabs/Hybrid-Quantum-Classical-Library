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

#include "SymbolicOperatorUtils.hpp"

#include <armadillo>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>

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

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>
    Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
typedef int vertices_size_type;
typedef boost::property_map<Graph, boost::vertex_index_t>::const_type
    vertex_index_map;
typedef std::pair<int, int> Edge;

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

vector<int> SymbolicOperatorUtils::getGroupsQWC(const Vec2DMat &qwcgraph) {
  // Based on
  // https://valelab4.ucsf.edu/svn/3rdpartypublic/boost/libs/graph/doc/sequential_vertex_coloring.html

  int numNodes = qwcgraph.size();

  // Transfer Vec2DMat to an edge list
  vector<Edge> edges;
  for (int i = 0; i < qwcgraph.size(); i++) {
    for (int j = 0; j < i; j++) {
      if (qwcgraph[i][j] == 1) {
        edges.push_back(Edge(i, j));
      }
    }
  }

  // Define the graph
  Graph g(edges.begin(), edges.end(), numNodes);

  std::vector<vertices_size_type> color_vec(boost::num_vertices(g));
  boost::iterator_property_map<vertices_size_type *, vertex_index_map> color(
      &color_vec.front(), get(boost::vertex_index, g));
  vertices_size_type num_colors = sequential_vertex_coloring(g, color);

  return color_vec;
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
  arma::Row<double> vR = {1}; // To maintain the shape of the matrix
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
  double exp_val = 0.0;

  for (const auto &pstr : v_pstr) {
    exp_val += getExpectValSglPauli(pstr, ProbReg, num_qbits);
  }

  return exp_val;
}

double SymbolicOperatorUtils::getExpectValSetOfPaulis(
    SymbolicOperator &symbop, const std::set<pstring> &s_pstr,
    const std::vector<double> ProbReg, int num_qbits, double eps) {
  double exp_val = 0.0;

  for (const auto &pstr : s_pstr) {
    exp_val += symbop.op_sum[pstr].real() *
               getExpectValSglPauli(pstr, ProbReg, num_qbits);
  }

  return exp_val;
}

double SymbolicOperatorUtils::getExpectVal(SymbolicOperator &symbop,
                                           const std::vector<double> ProbReg,
                                           int num_qbits, double eps,
                                           METHOD method) {
  double exp_val = 0.0;

  for (const auto &pstr : symbop.getOrderedPStringList()) {
    exp_val += symbop.op_sum[pstr].real() *
               getExpectValSglPauli(pstr, ProbReg, num_qbits);
  }

  return exp_val;
}

void SymbolicOperatorUtils::applyBasisChange(
    const pstring &pstr, std::vector<double> &variable_params, int num_qbits) {
  if (variable_params.empty()) {
    variable_params = std::vector<double>(num_qbits * 2, 0);
  }

  for (auto idx = 0; idx < num_qbits; ++idx) {
    auto foundX = pstr.find(op_pair(idx, 'X'));
    auto foundY = pstr.find(op_pair(idx, 'Y'));
    auto foundZ = pstr.find(op_pair(idx, 'Z'));

    int basis_step_diff = 2;
    if (foundX != pstr.end()) {
      variable_params[basis_step_diff * idx] = FP_PIby2;
      variable_params[basis_step_diff * idx + 1] = 0;
    } else if (foundY != pstr.end()) {
      variable_params[basis_step_diff * idx] = 0;
      variable_params[basis_step_diff * idx + 1] = FP_PIby2;
    } else if (foundZ != pstr.end()) {
      variable_params[basis_step_diff * idx] = 0;
      variable_params[basis_step_diff * idx + 1] = 0;
    }
  }
}

QWCMap SymbolicOperatorUtils::getQubitwiseCommutationGroups(
    const SymbolicOperator &symbop, int num_qbits) {

  Vec2DMat qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(symbop);
  std::vector<int> coloring = SymbolicOperatorUtils::getGroupsQWC(qwcmat);
  QWCMap qbtwise_comm_groups;
  for (auto pos = 0; pos < coloring.size(); pos++) {
    const auto ord_list = symbop.getOrderedPStringList();
    qbtwise_comm_groups[coloring[pos]].emplace(ord_list[pos]);
  }

  return qbtwise_comm_groups;
}

void SymbolicOperatorUtils::applyBasisChange(
    const std::set<pstring> &s_pstr, std::vector<double> &variable_params,
    int num_qbits, bool qwc_check) {

  if (qwc_check) {
    // Check whether the set of pauli strings commutes with each other
    SymbolicOperator symbop;
    for (const auto &pstr : s_pstr) {
      symbop.addTerm(pstr);
    }
    QWCMap qwc_group_mapping = getQubitwiseCommutationGroups(symbop, num_qbits);
    if (qwc_group_mapping.size() > 1) {
      throw std::logic_error("Provided pauli strings does not qubitwise "
                             "commute. Basis change cannot be applied.");
    }
  }

  if (variable_params.empty()) {
    variable_params = std::vector<double>(num_qbits * 2, 0);
  }

  for (const auto &pstr : s_pstr) {
    applyBasisChange(pstr, variable_params, num_qbits);
  }
}

SymbolicOperator
SymbolicOperatorUtils::getFoldedHamiltonian(SymbolicOperator &symbop,
                                                  double gamma) {
  // (H - γI)^2 = H^2 - 2*γI*H + γI^I
  // Create a SymbolicOperator object for gamma
  SymbolicOperator gamma_symbop;
  gamma_symbop.addIdentTerm(gamma);

  // Create a SymbolicOperator object for coefficient 2
  SymbolicOperator coeff;
  coeff.addIdentTerm(-2);

  // Create a SymbolicOperator object for negative identity to propagate sign
  SymbolicOperator H_folded =
      symbop * symbop +
      (coeff * symbop * gamma_symbop) +
      gamma_symbop * gamma_symbop;
  return H_folded;
}

std::ostream &operator<<(std::ostream &s, const QWCMap &qwc_groups) {
  s << "\n";
  for (const auto &qwc_group : qwc_groups) {
    s << "Group " << qwc_group.first << "\n";
    for (const auto &s_pstr : qwc_group.second) {
      s << "\t[ ";
      for (const auto &pstr : s_pstr) {
        s << pstr.second << std::to_string(pstr.first) << " ";
      }
      s << "]\n";
    }
  }
  return s;
}

std::ostream &operator<<(std::ostream &s, const std::set<pstring> &qwc_group) {
  s << "\n";
  for (const auto &s_pstr : qwc_group) {
    s << "\t[ ";
    for (const auto &pstr : s_pstr) {
      s << pstr.second << std::to_string(pstr.first) << " ";
    }
    s << "]\n";
  }
  return s;
}

} // namespace core
} // namespace quantum
} // namespace hybrid
