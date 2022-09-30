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

#include <gtest/gtest.h>

#include <armadillo>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using op_pair = std::pair<int, char>;
using pstring = std::set<op_pair>;
using ComplexDP = std::complex<double>;
using MapPString = std::map<pstring, ComplexDP>;
using Vec2DMat = std::vector<std::vector<int>>;

using namespace std;
using namespace arma;

using namespace hybrid::quantum::core;

TEST(QWCTests, Complete) {

  // Test less-than & greater-than ops for pstring
  /*
  pstring p1 = {{0, 'Z'}};
  pstring p2 = {{1, 'X'},{2,'Y'},{3,'Z'}};
  pstring p3 = {{1, 'X'},{0,'Y'},{3,'Z'}};
  pstring p4 = {{1,'Y'},        {3,'X'},{4,'X'}};
  pstring p5 = {{1,'Y'},{2,'Y'},{3,'X'},{4,'X'}};
  cout << "(p1 > p2)? " << (p1 > p2) << endl;
  cout << "(p2 > p1)? " << (p2 > p1) << endl;
  cout << "---" << endl;
  cout << "(p3 < p1)? " << (p3 < p1) << endl;
  cout << "(p4<p5)? " << (p4 < p5) << endl;
  */

  Vec2DMat qwcmat;
  pstring inp;
  int num_qubits = 1;

  // Empty op
  SymbolicOperator op;
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ(qwcmat, Vec2DMat({}));

  // Identity op
  SymbolicOperator op_id;
  op_id.addIdentTerm(42.0);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op_id);
  EXPECT_EQ(qwcmat, Vec2DMat({{0}}));

  // Op with single pauli string of length-1
  op.removeAllTerms();
  inp = {{0, 'Z'}};
  op.addTerm(inp, 1.);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ(qwcmat, Vec2DMat({{0}}));

  // Op with single longer pauli string
  op.removeAllTerms();
  inp = {{0, 'Z'}, {2, 'X'}, {3, 'X'}};
  op.addTerm(inp, 1.);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ(qwcmat, Vec2DMat({{0}}));

  // op = 2*Z0 + 3*Z1 (another sum of Paulis)
  op.removeAllTerms();
  inp = {{0, 'Z'}};
  op.addTerm(inp, 2);
  inp = {{1, 'Z'}};
  op.addTerm(inp, 3);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ(qwcmat, Vec2DMat({{0, 0}, {0, 0}}));

  // op = 2*Z0 + 3*Y0 (another sum of Paulis)
  num_qubits = 1;
  op.removeAllTerms();
  inp = {{0, 'Z'}};
  op.addTerm(inp, 2);
  inp = {{0, 'Y'}};
  op.addTerm(inp, 3);

  QWCMap qwc_groups_mapping;
  std::vector<std::vector<double>> v_variable_params;
  double actual_cost = 0.0;

  std::set<pstring> y0{{make_pair<int, char>(0, 'Y')}};
  std::set<pstring> z0{{make_pair<int, char>(0, 'Z')}};
  QWCMap gold_qwc_groups_mapping{{0, y0}, {1, z0}};
  vector<int> gold_coloring_small = {0, 1};
  std::vector<std::vector<double>> gold_v_variable_params{
      {0, 1.5707963267948966}, {0, 0}};
  double expected_cost = 1.9978159519295624;

  std::vector<double> y0_prob_reg{0.0000000004756236, 0.0000000005184095};
  std::vector<double> z0_prob_reg{0.9994539875174638, 0.0005460114885037};
  std::map<int, std::vector<double>> precomputed_prob_reg_map{{0, y0_prob_reg},
                                                              {1, z0_prob_reg}};

  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ(qwcmat, Vec2DMat({{0, 1}, {1, 0}}));

  vector<int> coloring_small = SymbolicOperatorUtils::getGroupsQWC(qwcmat);
  EXPECT_EQ(coloring_small, gold_coloring_small);

  qwc_groups_mapping =
      SymbolicOperatorUtils::getQubitwiseCommutationGroups(op, num_qubits);
  ASSERT_TRUE(qwc_groups_mapping.size() == gold_qwc_groups_mapping.size());
  EXPECT_EQ(qwc_groups_mapping, gold_qwc_groups_mapping);

  for (auto &qwc_group_mapping : qwc_groups_mapping) {
    std::vector<double> variable_params;
    SymbolicOperatorUtils::applyBasisChange(qwc_group_mapping.second,
                                            variable_params, num_qubits);
    v_variable_params.push_back(variable_params);
    actual_cost += SymbolicOperatorUtils::getExpectValSetOfPaulis(
        op, qwc_group_mapping.second,
        precomputed_prob_reg_map[qwc_group_mapping.first], num_qubits);
  }

  ASSERT_TRUE(v_variable_params.size() == gold_v_variable_params.size());
  EXPECT_EQ(v_variable_params, gold_v_variable_params);
  EXPECT_DOUBLE_EQ(actual_cost, expected_cost);

  // op_4 = 2*Z0 + 3*Y0 + 3*Z1   (another sum of Paulis)
  op.removeAllTerms();
  inp = {{0, 'Z'}};
  op.addTerm(inp, 2);
  inp = {{0, 'Y'}};
  op.addTerm(inp, 3);
  inp = {{1, 'Z'}};
  op.addTerm(inp, 3);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ(qwcmat, Vec2DMat({{0, 1, 0}, {1, 0, 0}, {0, 0, 0}}));

  // op_4 = 2*Z0 + 3*Y0 + 3*Z1   (another sum of Paulis)
  // Changing order in construction (ought to be same result)
  op.removeAllTerms();
  inp = {{1, 'Z'}};
  op.addTerm(inp, 3);
  inp = {{0, 'Z'}};
  op.addTerm(inp, 2);
  inp = {{0, 'Y'}};
  op.addTerm(inp, 3);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ(qwcmat, Vec2DMat({{0, 1, 0}, {1, 0, 0}, {0, 0, 0}}));

  // Example from arXiv:1907.03358 ([...] Izmaylov)
  // Deliberately placing in random order
  // Alphabetical order is:
  // [Y1 X3 X4]
  // [Y1 Y2 X3 X4]
  // [Z1]
  // [Z1 Z2]
  // [Z1 Z2 Z3]
  // [Z1 Z2 Z3 Z4]
  // [X3 X4]
  // cout << "*******" << endl;
  op.removeAllTerms();
  vector<string> inp_char;
  inp_char = {"X3", "X4"};
  op.addTerm(inp_char, 1.0);
  inp_char = {"Z1", "Z2"};
  op.addTerm(inp_char, 1.0);
  inp_char = {"Y1", "X3", "X4"};
  op.addTerm(inp_char, 1.0);
  inp_char = {"Z1", "Z2", "Z3"};
  op.addTerm(inp_char, 1.0);
  inp_char = {"Z1", "Z2", "Z3", "Z4"};
  op.addTerm(inp_char, 1.0);
  inp_char = {"Z1"};
  op.addTerm(inp_char, 1.0);
  inp_char = {"Y1", "Y2", "X3", "X4"};
  op.addTerm(inp_char, 1.0);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);

  Vec2DMat goldMat = {{0, 0, 1, 1, 1, 1, 0}, {0, 0, 1, 1, 1, 1, 0},
                      {1, 1, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 0, 0, 0},
                      {1, 1, 0, 0, 0, 0, 1}, {1, 1, 0, 0, 0, 0, 1},
                      {0, 0, 0, 0, 1, 1, 0}};
  EXPECT_EQ(qwcmat, goldMat);

  vector<int> coloring = SymbolicOperatorUtils::getGroupsQWC(qwcmat);
  vector<int> gold_coloring = {0, 0, 1, 1, 1, 1, 0};
  EXPECT_EQ(coloring, gold_coloring);

  num_qubits = 4;
  gold_qwc_groups_mapping = {{0,
                              {{{1, 'Y'}, {2, 'Y'}, {3, 'X'}, {4, 'X'}},
                               {{1, 'Y'}, {3, 'X'}, {4, 'X'}},
                               {{3, 'X'}, {4, 'X'}}}},
                             {1,
                              {{{1, 'Z'}},
                               {{1, 'Z'}, {2, 'Z'}},
                               {{1, 'Z'}, {2, 'Z'}, {3, 'Z'}},
                               {{1, 'Z'}, {2, 'Z'}, {3, 'Z'}, {4, 'Z'}}}}};
  gold_v_variable_params = {
      {0, 1.5707963267948966, 0, 1.5707963267948966, 1.5707963267948966, 0, 0,
       1.5707963267948966, 1.5707963267948966, 0, 1.5707963267948966, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  std::vector<double> grp1_prob_reg{
      0.0624998969931823, 0.0625058893829035, 0.0625058888083628,
      0.0624879133625461, 0.0625027854931024, 0.0625087781597683,
      0.0625087775852011, 0.0624908013086298, 0.0624972149424357,
      0.0625032070750061, 0.0625032065004901, 0.0624852318260516,
      0.0625001019970131, 0.0625060944063897, 0.0625060938318472,
      0.0624881183270698};
  std::vector<double> grp2_prob_reg{
      0.0000000022587437, 0.0000000000000000, 0.0000000000000000,
      0.0000000000000000, 0.0170722599322310, 0.0000000000000000,
      0.0000000000000000, 0.0000000000392311, 0.0000000000392026,
      0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
      0.9829277354718768, 0.0000000000000000, 0.0000000000000000,
      0.0000000022587152};
  precomputed_prob_reg_map = {{0, grp1_prob_reg}, {1, grp2_prob_reg}};
  expected_cost = -3.9999520436097669;

  v_variable_params.clear();
  qwc_groups_mapping.clear();
  actual_cost = 0.0;

  qwc_groups_mapping =
      SymbolicOperatorUtils::getQubitwiseCommutationGroups(op, num_qubits);
  ASSERT_TRUE(qwc_groups_mapping.size() == gold_qwc_groups_mapping.size());
  EXPECT_EQ(qwc_groups_mapping, gold_qwc_groups_mapping);

  for (auto &qwc_group_mapping : qwc_groups_mapping) {
    std::vector<double> variable_params;
    SymbolicOperatorUtils::applyBasisChange(qwc_group_mapping.second,
                                            variable_params, num_qubits);
    v_variable_params.push_back(variable_params);
    actual_cost += SymbolicOperatorUtils::getExpectValSetOfPaulis(
        op, qwc_group_mapping.second,
        precomputed_prob_reg_map[qwc_group_mapping.first], num_qubits);
  }

  EXPECT_EQ(v_variable_params.size(), gold_v_variable_params.size());
  EXPECT_EQ(v_variable_params, gold_v_variable_params);
  EXPECT_DOUBLE_EQ(actual_cost, expected_cost);

  // op = 0.5 * Z0 * Z1 + 0.5 * X0 + 0.25 * X1
  num_qubits = 2;
  op.removeAllTerms();
  inp = {{0, 'Z'}, {1, 'Z'}};
  op.addTerm(inp, 0.5);
  inp = {{0, 'X'}};
  op.addTerm(inp, 0.5);
  inp = {{1, 'X'}};
  op.addTerm(inp, 0.25);

  v_variable_params.clear();
  qwc_groups_mapping.clear();
  actual_cost = 0.0;

  set<pstring> x0x1 = {{make_pair<int, char>(0, 'X')},
                       {make_pair<int, char>(1, 'X')}};
  set<pstring> z0z1 = {
      {make_pair<int, char>(0, 'Z'), make_pair<int, char>(1, 'Z')}};
  gold_qwc_groups_mapping = {{0, x0x1}, {1, z0z1}};
  gold_coloring = {0, 1, 0};
  gold_v_variable_params = {{1.5707963267948966, 0, 1.5707963267948966, 0},
                            {0, 0, 0, 0}};
  expected_cost = -0.90138667910371084;

  std::vector<double> x0x1_prob_reg{0.0839878007634251, 0.0000000758196198,
                                    0.0000010586065377, 0.9160110648104176};
  std::vector<double> z0z1_prob_reg{0.1108608098371496, 0.3891550046005872,
                                    0.3882139807711494, 0.1117702047911133};
  precomputed_prob_reg_map = {{0, x0x1_prob_reg}, {1, z0z1_prob_reg}};

  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ(qwcmat, Vec2DMat({{0, 1, 0}, {1, 0, 1}, {0, 1, 0}}));

  coloring = SymbolicOperatorUtils::getGroupsQWC(qwcmat);
  EXPECT_EQ(coloring, gold_coloring);

  qwc_groups_mapping =
      SymbolicOperatorUtils::getQubitwiseCommutationGroups(op, num_qubits);
  std::cout << "Qubitwise Commutation(QWC) Groups:"
            << qwc_groups_mapping << std::endl;
  ASSERT_TRUE(qwc_groups_mapping.size() == gold_qwc_groups_mapping.size());
  EXPECT_EQ(qwc_groups_mapping, gold_qwc_groups_mapping);

  for (auto &qwc_group_mapping : qwc_groups_mapping) {
    std::vector<double> variable_params;
    SymbolicOperatorUtils::applyBasisChange(qwc_group_mapping.second,
                                            variable_params, num_qubits);
    v_variable_params.push_back(variable_params);
    actual_cost += SymbolicOperatorUtils::getExpectValSetOfPaulis(
        op, qwc_group_mapping.second,
        precomputed_prob_reg_map[qwc_group_mapping.first], num_qubits);
  }
  ASSERT_TRUE(v_variable_params.size() == gold_v_variable_params.size());
  EXPECT_EQ(v_variable_params, gold_v_variable_params);
  EXPECT_DOUBLE_EQ(actual_cost, expected_cost);

  // pauli strings does not commute
  num_qubits = 2;
  std::vector<double> variable_params;
  std::set<pstring> s_pstr;
  inp = {{0, 'Z'}, {1, 'Z'}};
  s_pstr.insert(inp);
  inp = {{0, 'X'}};
  s_pstr.insert(inp);
  EXPECT_THROW(SymbolicOperatorUtils::applyBasisChange(s_pstr, variable_params,
                                                       num_qubits, true),
               std::logic_error);
}

TEST(UtilTests, Complete) {
  string charstring;

  // Z0 Y1
  pstring inp_Z0Y1 = {{0, 'Z'}, {1, 'Y'}};
  charstring = SymbolicOperatorUtils::getCharString_pstring(inp_Z0Y1);
  EXPECT_EQ(charstring, "Z0 Y1 ");

  // X0 Y1 Z2
  pstring inp_X0Y1Z2 = {{0, 'X'}, {1, 'Y'}, {2, 'Z'}};
  EXPECT_EQ(SymbolicOperatorUtils::findFirstPauliStringBasis(inp_X0Y1Z2), 'X');

  // Y0 X1 Z2
  pstring inp_Y0X1Z2 = {{0, 'Y'}, {1, 'X'}, {2, 'Z'}};
  EXPECT_EQ(SymbolicOperatorUtils::findFirstPauliStringBasis(inp_Y0X1Z2), 'Y');

  // Z0 Y1 X2
  pstring inp_Z0Y1X2 = {{0, 'Z'}, {1, 'Y'}, {2, 'X'}};
  EXPECT_EQ(SymbolicOperatorUtils::findFirstPauliStringBasis(inp_Z0Y1X2), 'Z');

  // Empty Pauli string
  pstring inp_empty = {};
  EXPECT_THROW(SymbolicOperatorUtils::findFirstPauliStringBasis(inp_empty),
               std::logic_error);

  // Incorrect Pauli string
  pstring inp_incorrect = {{0, 'A'}};
  EXPECT_THROW(SymbolicOperatorUtils::findFirstPauliStringBasis(inp_incorrect),
               std::invalid_argument);

  std::vector<double> variable_params;
  std::vector<double> e_variable_params{0.0, 0.0, 0.0, FP_PIby2};
  int num_qbits = 2;

  SymbolicOperatorUtils::applyBasisChange(inp_Z0Y1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);

  variable_params = {};
  e_variable_params = {FP_PIby2, 0.0, 0.0, FP_PIby2};
  pstring inp_X0Y1 = {{0, 'X'}, {1, 'Y'}};
  SymbolicOperatorUtils::applyBasisChange(inp_X0Y1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);

  variable_params = {};
  e_variable_params = {0.0, FP_PIby2, FP_PIby2, 0.0};
  pstring inp_Y0X1 = {{0, 'Y'}, {1, 'X'}};
  SymbolicOperatorUtils::applyBasisChange(inp_Y0X1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);

  variable_params = {};
  e_variable_params = {FP_PIby2, 0.0, FP_PIby2, 0.0};
  pstring inp_X1 = {{1, 'X'}};
  SymbolicOperatorUtils::applyBasisChange(inp_X1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);

  variable_params = {};
  e_variable_params = {0.0, FP_PIby2, 0.0, FP_PIby2};
  pstring inp_Y1 = {{1, 'Y'}};
  SymbolicOperatorUtils::applyBasisChange(inp_Y1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);

  variable_params = {};
  e_variable_params = {0.0, 0.0, 0.0, 0.0};
  pstring inp_Z1 = {{1, 'Z'}};
  SymbolicOperatorUtils::applyBasisChange(inp_Z1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);
}
