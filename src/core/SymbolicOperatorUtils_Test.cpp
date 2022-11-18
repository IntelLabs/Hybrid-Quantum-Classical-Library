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
  QWCMap qwc_groups_mapping;
  QWCMap gold_qwc_groups_mapping;
  std::vector<std::vector<double>> v_variable_params;
  double actual_cost = 0.0;
  vector<int> coloring;
  vector<int> gold_coloring;

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

  std::set<pstring> y0{{make_pair<int, char>(0, 'Y')}};
  std::set<pstring> z0{{make_pair<int, char>(0, 'Z')}};
  gold_qwc_groups_mapping = {{0, y0}, {1, z0}};
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

  // Single Qubit QWC
  // H = 0.5 * Y0 + 0.5 * X0
  num_qubits = 1;
  SymbolicOperator orig_H;
  pstring inp_y1{{0, 'Y'}};
  orig_H.addTerm(inp_y1, 0.5);

  pstring inp_y2{{0, 'X'}};
  orig_H.addTerm(inp_y2, 0.25);

  v_variable_params.clear();
  qwc_groups_mapping.clear();
  actual_cost = 0.0;

  set<pstring> x = {{make_pair<int, char>(0, 'X')}};
  set<pstring> y = {{make_pair<int, char>(0, 'Y')}};
  gold_qwc_groups_mapping = {{0, x}, {1, y}};
  gold_coloring = {0, 1};
  gold_v_variable_params = {{1.5707963267948966, 0}, {0, 1.5707963267948966}};
  expected_cost = -0.55903842387228087;

  std::vector<double> x_prob_reg{0.27643693541442487804,
                                 0.72356306458557562156};
  std::vector<double> y_prob_reg{0.052743108420506820688,
                                 0.94725689157949322095};
  precomputed_prob_reg_map = {{0, x_prob_reg}, {1, y_prob_reg}};

  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(orig_H);
  EXPECT_EQ(qwcmat, Vec2DMat({{0, 1}, {1, 0}}));

  coloring = SymbolicOperatorUtils::getGroupsQWC(qwcmat);
  EXPECT_EQ(coloring, gold_coloring);

  qwc_groups_mapping =
      SymbolicOperatorUtils::getQubitwiseCommutationGroups(orig_H, num_qubits);
  std::cout << "Qubitwise Commutation(QWC) Groups:" << qwc_groups_mapping
            << std::endl;
  ASSERT_TRUE(qwc_groups_mapping.size() == gold_qwc_groups_mapping.size());
  EXPECT_EQ(qwc_groups_mapping, gold_qwc_groups_mapping);

  for (auto &qwc_group_mapping : qwc_groups_mapping) {
    std::vector<double> variable_params;
    SymbolicOperatorUtils::applyBasisChange(qwc_group_mapping.second,
                                            variable_params, num_qubits);
    v_variable_params.push_back(variable_params);
    actual_cost += SymbolicOperatorUtils::getExpectValSetOfPaulis(
        orig_H, qwc_group_mapping.second,
        precomputed_prob_reg_map[qwc_group_mapping.first], num_qubits);
  }
  ASSERT_TRUE(v_variable_params.size() == gold_v_variable_params.size());
  EXPECT_EQ(v_variable_params, gold_v_variable_params);
  EXPECT_DOUBLE_EQ(actual_cost, expected_cost);

  // Example from arXiv:1907.03358 ([...] Izmaylov)
  // Deliberately placing in random order
  // Alphabetical order is:
  // [Y0 X2 X3]
  // [Y0 Y1 X2 X3]
  // [Z0]
  // [Z0 Z1]
  // [Z0 Z1 Z2]
  // [Z0 Z1 Z2 Z3]
  // [X2 X3]
  // cout << "*******" << endl;
  op.removeAllTerms();
  vector<string> inp_char;
  inp_char = {"X2", "X3"};
  op.addTerm(inp_char, 1.0);
  inp_char = {"Z0", "Z1"};
  op.addTerm(inp_char, 1.0);
  inp_char = {"Y0", "X2", "X3"};
  op.addTerm(inp_char, 1.0);
  inp_char = {"Z0", "Z1", "Z2"};
  op.addTerm(inp_char, 1.0);
  inp_char = {"Z0", "Z1", "Z2", "Z3"};
  op.addTerm(inp_char, 1.0);
  inp_char = {"Z0"};
  op.addTerm(inp_char, 1.0);
  inp_char = {"Y0", "Y1", "X2", "X3"};
  op.addTerm(inp_char, 1.0);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);

  Vec2DMat goldMat = {{0, 0, 1, 1, 1, 1, 0}, {0, 0, 1, 1, 1, 1, 0},
                      {1, 1, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 0, 0, 0},
                      {1, 1, 0, 0, 0, 0, 1}, {1, 1, 0, 0, 0, 0, 1},
                      {0, 0, 0, 0, 1, 1, 0}};
  EXPECT_EQ(qwcmat, goldMat);

  coloring = SymbolicOperatorUtils::getGroupsQWC(qwcmat);
  gold_coloring = {0, 0, 1, 1, 1, 1, 0};
  EXPECT_EQ(coloring, gold_coloring);

  num_qubits = 4;
  gold_qwc_groups_mapping = {{0,
                              {{{0, 'Y'}, {1, 'Y'}, {2, 'X'}, {3, 'X'}},
                               {{0, 'Y'}, {2, 'X'}, {3, 'X'}},
                               {{2, 'X'}, {3, 'X'}}}},
                             {1,
                              {{{0, 'Z'}},
                               {{0, 'Z'}, {1, 'Z'}},
                               {{0, 'Z'}, {1, 'Z'}, {2, 'Z'}},
                               {{0, 'Z'}, {1, 'Z'}, {2, 'Z'}, {3, 'Z'}}}}};
  gold_v_variable_params = {{0, 1.5707963267948965580, 0, 1.5707963267948965580,
                             1.5707963267948965580, 0, 1.5707963267948965580,
                             0},
                            {0, 0, 0, 0, 0, 0, 0, 0}};
  std::vector<double> grp1_prob_reg{
      1.3845591418918789389e-02, 2.6297890660156664211e-01,
      2.6308482979159647508e-01, 1.3837970684570516378e-02,
      2.2232302313607353421e-03, 4.2227351485173952872e-02,
      4.2244359905482704864e-02, 2.2220065460391343753e-03,
      1.6569726233180367709e-03, 3.1472028573188495781e-02,
      3.1484704942203262101e-02, 1.6560606111256648276e-03,
      7.2776530687329223784e-03, 1.3822950488238314182e-01,
      1.3828518124936861611e-01, 7.2736473849713805667e-03};
  std::vector<double> grp2_prob_reg{
      1.0813809754668755186e-02, 6.0852570919786738219e-11,
      1.1042754839469855444e-10, 4.2486590906362607656e-03,
      9.5129950703424948077e-02, 5.3532494126565204755e-10,
      9.7143998955154058212e-10, 3.7375794379347855589e-02,
      6.0717052170840069003e-01, 3.4167317597493791952e-09,
      6.2002526113279708641e-09, 2.3855242649418809120e-01,
      4.8164729282900208343e-03, 2.7103746699951031567e-11,
      4.9184451127491999486e-11, 1.8923535697257763688e-03};
  precomputed_prob_reg_map = {{0, grp1_prob_reg}, {1, grp2_prob_reg}};
  expected_cost = -4.8121480259189262085;

  v_variable_params.clear();
  qwc_groups_mapping.clear();
  actual_cost = 0.0;

  qwc_groups_mapping =
      SymbolicOperatorUtils::getQubitwiseCommutationGroups(op, num_qubits);
  std::cout << "Qubitwise Commutation(QWC) Groups:" << qwc_groups_mapping
            << std::endl;
  ASSERT_TRUE(qwc_groups_mapping.size() == gold_qwc_groups_mapping.size());
  EXPECT_EQ(qwc_groups_mapping, gold_qwc_groups_mapping);

  for (auto &qwc_group_mapping : qwc_groups_mapping) {
    std::vector<double> variable_params;
    SymbolicOperatorUtils::applyBasisChange(qwc_group_mapping.second,
                                            variable_params, num_qubits);
    v_variable_params.push_back(variable_params);

    double current_cost = SymbolicOperatorUtils::getExpectValSetOfPaulis(
        op, qwc_group_mapping.second,
        precomputed_prob_reg_map[qwc_group_mapping.first], num_qubits);
    actual_cost += current_cost;
  }

  EXPECT_EQ(v_variable_params.size(), gold_v_variable_params.size());
  EXPECT_EQ(v_variable_params, gold_v_variable_params);
  EXPECT_DOUBLE_EQ(actual_cost, expected_cost);

  // Single Quantum Kernel call w/ static parameters
  grp1_prob_reg = {8.0360546671606790728e-04, 1.5224091907570652671e-02,
                   4.4994059027517307320e-03, 2.7466122575537224582e-02,
                   8.0374917510785310899e-04, 1.5226814424844437382e-02,
                   4.5002105294162644475e-03, 2.7471034329460592566e-02,
                   3.5504575814976936654e-03, 6.7262474900198779282e-02,
                   1.9879095478210224296e-02, 1.2134972591857257074e-01,
                   1.1586331401959174339e-02, 2.1949996788327180708e-01,
                   6.4872141940806551941e-02, 3.9600477058407790310e-01};
  grp2_prob_reg = {2.4609250508941676872e-01, 1.1149352807346514516e-01,
                   8.8246630779790516397e-02, 5.0244452128424550719e-02,
                   1.1149352807346521455e-01, 5.0512740311826523354e-02,
                   3.9980608928583002970e-02, 2.2763518262699379557e-02,
                   8.8246630779790558030e-02, 3.9980608928582968276e-02,
                   3.1644473858133674582e-02, 1.8017223296170257335e-02,
                   5.0244452128424564596e-02, 2.2763518262699379557e-02,
                   1.8017223296170295499e-02, 1.0258357802356793051e-02};
  precomputed_prob_reg_map = {{0, grp1_prob_reg}, {1, grp2_prob_reg}};
  expected_cost = 8.7993924085391994616e-01;

  v_variable_params.clear();
  actual_cost = 0.0;

  for (auto &qwc_group_mapping : qwc_groups_mapping) {
    std::vector<double> variable_params;
    SymbolicOperatorUtils::applyBasisChange(qwc_group_mapping.second,
                                            variable_params, num_qubits);
    v_variable_params.push_back(variable_params);

    double current_cost = SymbolicOperatorUtils::getExpectValSetOfPaulis(
        op, qwc_group_mapping.second,
        precomputed_prob_reg_map[qwc_group_mapping.first], num_qubits);
    actual_cost += current_cost;
  }

  EXPECT_EQ(v_variable_params.size(), gold_v_variable_params.size());
  EXPECT_EQ(v_variable_params, gold_v_variable_params);
  EXPECT_DOUBLE_EQ(actual_cost, expected_cost);

  // Individual costs for each Pauli string in Z0 + Z0Z1 + Z0Z1Z2 + Z0Z1Z2Z3 +
  // X2X3 + Y0X2X3 + Y0Y1X2X3
  std::map<pstring, double> expected_pauli_str_cost{
      {{{0, 'Y'}, {1, 'Y'}, {2, 'X'}, {3, 'X'}}, -6.2067338507351621502e-01},
      {{{0, 'Y'}, {2, 'X'}, {3, 'X'}}, -2.5679956294393513350e-01},
      {{{2, 'X'}, {3, 'X'}}, -9.0001373486192592921e-01},
      {{{0, 'Z'}}, -7.0486356878775413559e-01},
      {{{0, 'Z'}, {1, 'Z'}}, -9.5645740881822216561e-01},
      {{{0, 'Z'}, {1, 'Z'}, {2, 'Z'}}, -4.1688297836771287530e-01},
      {{{0, 'Z'}, {1, 'Z'}, {2, 'Z'}, {3, 'Z'}}, -9.5645738706585958777e-01}};

  std::map<pstring, std::vector<double>> precomputed_pauli_str_probs = {
      {{{0, 'Y'}, {1, 'Y'}, {2, 'X'}, {3, 'X'}},
       {1.3845591418918763368e-02, 2.6297890660156653109e-01,
        2.6308482979159647508e-01, 1.3837970684570512908e-02,
        2.2232302313607314390e-03, 4.2227351485173932055e-02,
        4.2244359905482718742e-02, 2.2220065460391309059e-03,
        1.6569726233180350362e-03, 3.1472028573188495781e-02,
        3.1484704942203255162e-02, 1.6560606111256665623e-03,
        7.2776530687329111027e-03, 1.3822950488238311406e-01,
        1.3828518124936853284e-01, 7.2736473849713831688e-03}},
      {{{0, 'Y'}, {2, 'X'}, {3, 'X'}},
       {1.3845591418918778981e-02, 2.6297890660156653109e-01,
        2.6308482979159647508e-01, 1.3837970684570516378e-02,
        2.2232302313607240664e-03, 4.2227351485173918177e-02,
        4.2244359905482725681e-02, 2.2220065460391213649e-03,
        1.6569726233180283141e-03, 3.1472028573188447209e-02,
        3.1484704942203227407e-02, 1.6560606111256552866e-03,
        7.2776530687328989597e-03, 1.3822950488238303079e-01,
        1.3828518124936861611e-01, 7.2736473849713640868e-03}},
      {{{2, 'X'}, {3, 'X'}},
       {1.3845591418918778981e-02, 2.6297890660156653109e-01,
        2.6308482979159647508e-01, 1.3837970684570512908e-02,
        2.2232302313607149591e-03, 4.2227351485173938994e-02,
        4.2244359905482718742e-02, 2.2220065460391131250e-03,
        1.6569726233180265794e-03, 3.1472028573188516598e-02,
        3.1484704942203234346e-02, 1.6560606111256559372e-03,
        7.2776530687329015618e-03, 1.3822950488238311406e-01,
        1.3828518124936853284e-01, 7.2736473849713614848e-03}},
      {{{0, 'Z'}},
       {1.0813809754668784677e-02, 6.0852570919554171220e-11,
        1.1042754839434974272e-10, 4.2486590906362503572e-03,
        9.5129950703424975833e-02, 5.3532494126360580951e-10,
        9.7143998954847070807e-10, 3.7375794379347827834e-02,
        6.0717052170840069003e-01, 3.4167317597363159454e-09,
        6.2002526113083782641e-09, 2.3855242649418798018e-01,
        4.8164729282900147628e-03, 2.7103746699847401346e-11,
        4.9184451127336689864e-11, 1.8923535697257776698e-03}},
      {{{0, 'Z'}, {1, 'Z'}},
       {1.0813809754668770799e-02, 6.0852570919902918321e-11,
        1.1042754839492054904e-10, 4.2486590906362590309e-03,
        9.5129950703425114611e-02, 5.3532494126667475298e-10,
        9.7143998955349210798e-10, 3.7375794379347820895e-02,
        6.0717052170840102310e-01, 3.4167317597559073046e-09,
        6.2002526113404348216e-09, 2.3855242649418803569e-01,
        4.8164729282900243038e-03, 2.7103746700002759436e-11,
        4.9184451127591022053e-11, 1.8923535697257763688e-03}},
      {{{0, 'Z'}, {1, 'Z'}, {2, 'Z'}},
       {1.0813809754668770799e-02, 6.0852570919770246306e-11,
        1.1042754839485190598e-10, 4.2486590906362607656e-03,
        9.5129950703425045222e-02, 5.3532494126550708415e-10,
        9.7143998955288847292e-10, 3.7375794379347807017e-02,
        6.0717052170840069003e-01, 3.4167317597484531665e-09,
        6.2002526113365743697e-09, 2.3855242649418798018e-01,
        4.8164729282900164975e-03, 2.7103746699943622484e-11,
        4.9184451127560274199e-11, 1.8923535697257709478e-03}},
      {{{0, 'Z'}, {1, 'Z'}, {2, 'Z'}, {3, 'Z'}},
       {1.0813809754668746513e-02, 6.0852570919857862827e-11,
        1.1042754839480294722e-10, 4.2486590906362477552e-03,
        9.5129950703424975833e-02, 5.3532494126627801648e-10,
        9.7143998955245854580e-10, 3.7375794379347807017e-02,
        6.0717052170840069003e-01, 3.4167317597533699281e-09,
        6.2002526113338273029e-09, 2.3855242649418792467e-01,
        4.8164729282900243038e-03, 2.7103746699982622758e-11,
        4.9184451127538476697e-11, 1.8923535697257672615e-03}}};

  for (auto &qwc_group_mapping : qwc_groups_mapping) {
    for (const auto &pstr : qwc_group_mapping.second) {
      double pauli_str_cost =
          op.op_sum[pstr].real() *
          SymbolicOperatorUtils::getExpectValSglPauli(
              pstr, precomputed_pauli_str_probs[pstr], num_qubits);
      EXPECT_EQ(pauli_str_cost, expected_pauli_str_cost[pstr]);
    }
  }

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
  std::cout << "Qubitwise Commutation(QWC) Groups:" << qwc_groups_mapping
            << std::endl;
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

  // Excited states - Single Qubit
  // Original Hamiltonian
  // X0 - 0.5 * Z0 + 1.5 * I0
  num_qubits = 1;
  orig_H.removeAllTerms();
  inp_y1 = {{0, 'X'}};
  orig_H.addTerm(inp_y1, 1.0);

  inp_y2 = {{0, 'Z'}};
  orig_H.addTerm(inp_y2, -0.5);

  orig_H.addIdentTerm(1.5);

  v_variable_params.clear();
  qwc_groups_mapping.clear();
  actual_cost = 0.0;

  x = {{}, {make_pair<int, char>(0, 'X')}};
  set<pstring> z = {{make_pair<int, char>(0, 'Z')}};
  gold_qwc_groups_mapping = {{0, x}, {1, z}};
  gold_coloring = {0, 0, 1};
  gold_v_variable_params = {{1.5707963267948966, 0}, {0, 0}};
  expected_cost = 0.38198754909266114;

  x_prob_reg = {0.052882488498675203625, 0.94711751150132517107};
  std::vector<double> z_prob_reg{0.72377742790468968526,
                                 0.27622257209531064781};
  precomputed_prob_reg_map = {{0, x_prob_reg}, {1, z_prob_reg}};

  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(orig_H);
  EXPECT_EQ(qwcmat, Vec2DMat({{0, 0, 0}, {0, 0, 1}, {0, 1, 0}}));

  coloring = SymbolicOperatorUtils::getGroupsQWC(qwcmat);
  EXPECT_EQ(coloring, gold_coloring);

  qwc_groups_mapping =
      SymbolicOperatorUtils::getQubitwiseCommutationGroups(orig_H, num_qubits);
  std::cout << "Qubitwise Commutation(QWC) Groups:" << qwc_groups_mapping
            << std::endl;
  ASSERT_TRUE(qwc_groups_mapping.size() == gold_qwc_groups_mapping.size());
  EXPECT_EQ(qwc_groups_mapping, gold_qwc_groups_mapping);

  for (auto &qwc_group_mapping : qwc_groups_mapping) {
    std::vector<double> variable_params;
    SymbolicOperatorUtils::applyBasisChange(qwc_group_mapping.second,
                                            variable_params, num_qubits);
    v_variable_params.push_back(variable_params);
    actual_cost += SymbolicOperatorUtils::getExpectValSetOfPaulis(
        orig_H, qwc_group_mapping.second,
        precomputed_prob_reg_map[qwc_group_mapping.first], num_qubits);
  }
  ASSERT_TRUE(v_variable_params.size() == gold_v_variable_params.size());
  EXPECT_EQ(v_variable_params, gold_v_variable_params);
  EXPECT_DOUBLE_EQ(actual_cost, expected_cost);

  // Folded Hamiltonian - Single qubit
  // -1.2 * X0 + 0.6 * Z0 + 1.61 * I0
  SymbolicOperator folded_H;
  pstring inp_y1_fold{{0, 'X'}};
  folded_H.addTerm(inp_y1_fold, -1.2);

  pstring inp_y2_fold{{0, 'Z'}};
  folded_H.addTerm(inp_y2_fold, 0.6);

  folded_H.addIdentTerm(1.61);

  v_variable_params.clear();
  qwc_groups_mapping.clear();
  actual_cost = 0.0;

  x = {{}, {make_pair<int, char>(0, 'X')}};
  z = {{make_pair<int, char>(0, 'Z')}};
  gold_qwc_groups_mapping = {{0, x}, {1, z}};
  gold_coloring = {0, 0, 1};
  gold_v_variable_params = {{1.5707963267948966, 0}, {0, 0}};
  expected_cost = 0.26838503344570969;

  x_prob_reg = {0.94712823921120792292, 0.052871760788792666885};
  z_prob_reg = {0.27624400629383893957, 0.72375599370616083839};
  precomputed_prob_reg_map = {{0, x_prob_reg}, {1, z_prob_reg}};

  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(folded_H);
  EXPECT_EQ(qwcmat, Vec2DMat({{0, 0, 0}, {0, 0, 1}, {0, 1, 0}}));

  coloring = SymbolicOperatorUtils::getGroupsQWC(qwcmat);
  EXPECT_EQ(coloring, gold_coloring);

  qwc_groups_mapping = SymbolicOperatorUtils::getQubitwiseCommutationGroups(
      folded_H, num_qubits);
  std::cout << "Qubitwise Commutation(QWC) Groups:" << qwc_groups_mapping
            << std::endl;
  ASSERT_TRUE(qwc_groups_mapping.size() == gold_qwc_groups_mapping.size());
  EXPECT_EQ(qwc_groups_mapping, gold_qwc_groups_mapping);

  for (auto &qwc_group_mapping : qwc_groups_mapping) {
    std::cout << "Calculating cost for group: " << qwc_group_mapping.second
              << std::endl;
    std::vector<double> variable_params;
    SymbolicOperatorUtils::applyBasisChange(qwc_group_mapping.second,
                                            variable_params, num_qubits);
    v_variable_params.push_back(variable_params);
    actual_cost += SymbolicOperatorUtils::getExpectValSetOfPaulis(
        folded_H, qwc_group_mapping.second,
        precomputed_prob_reg_map[qwc_group_mapping.first], num_qubits);
  }
  ASSERT_TRUE(v_variable_params.size() == gold_v_variable_params.size());
  EXPECT_EQ(v_variable_params, gold_v_variable_params);
  EXPECT_DOUBLE_EQ(actual_cost, expected_cost);

  // First Excited State Check
  v_variable_params.clear();
  actual_cost = 0.0;
  x_prob_reg = {0.94712823921120792292, 0.052871760788792666885};
  z_prob_reg = {0.27624400629383893957, 0.72375599370616083839};
  precomputed_prob_reg_map = {{0, x_prob_reg}, {1, z_prob_reg}};
  expected_cost = 2.6180124721285774;

  for (auto &qwc_group_mapping : qwc_groups_mapping) {
    std::vector<double> variable_params;
    SymbolicOperatorUtils::applyBasisChange(qwc_group_mapping.second,
                                            variable_params, num_qubits);
    v_variable_params.push_back(variable_params);
    actual_cost += SymbolicOperatorUtils::getExpectValSetOfPaulis(
        orig_H, qwc_group_mapping.second,
        precomputed_prob_reg_map[qwc_group_mapping.first], num_qubits);
  }
  ASSERT_TRUE(v_variable_params.size() == gold_v_variable_params.size());
  EXPECT_EQ(v_variable_params, gold_v_variable_params);
  EXPECT_DOUBLE_EQ(actual_cost, expected_cost);

  // Use getExcitedStateHamiltonian method w/ expected values from the above
  // test
  double gamma = 0.5;
  folded_H = SymbolicOperatorUtils::getFoldedHamiltonian(orig_H, gamma);

  EXPECT_EQ(qwcmat, Vec2DMat({{0, 0, 0}, {0, 0, 1}, {0, 1, 0}}));

  EXPECT_EQ(coloring, gold_coloring);

  ASSERT_TRUE(qwc_groups_mapping.size() == gold_qwc_groups_mapping.size());
  EXPECT_EQ(qwc_groups_mapping, gold_qwc_groups_mapping);

  ASSERT_TRUE(v_variable_params.size() == gold_v_variable_params.size());
  EXPECT_EQ(v_variable_params, gold_v_variable_params);
  EXPECT_DOUBLE_EQ(actual_cost, expected_cost);

  // First Excited State Check
  ASSERT_TRUE(v_variable_params.size() == gold_v_variable_params.size());
  EXPECT_EQ(v_variable_params, gold_v_variable_params);
  EXPECT_DOUBLE_EQ(actual_cost, expected_cost);
}

TEST(UtilTests, Complete) {
  string charstring;
  int num_qbits = 2;
  std::vector<double> variable_params;
  std::vector<double> e_variable_params;

  // Z0 Y1
  pstring inp_Z0Y1 = {{0, 'Z'}, {1, 'Y'}};
  charstring = SymbolicOperatorUtils::getCharString_pstring(inp_Z0Y1);
  EXPECT_EQ(charstring, "Z0 Y1 ");

  e_variable_params = {0.0, 0.0, 0.0, FP_PIby2};
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
  e_variable_params = {0, 0.0, FP_PIby2, 0.0};
  pstring inp_X1 = {{1, 'X'}};
  SymbolicOperatorUtils::applyBasisChange(inp_X1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);

  variable_params = {};
  e_variable_params = {0.0, 0.0, 0.0, FP_PIby2};
  pstring inp_Y1 = {{1, 'Y'}};
  SymbolicOperatorUtils::applyBasisChange(inp_Y1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);

  variable_params = {};
  e_variable_params = {0.0, 0.0, 0.0, 0.0};
  pstring inp_Z1 = {{1, 'Z'}};
  SymbolicOperatorUtils::applyBasisChange(inp_Z1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);
}
