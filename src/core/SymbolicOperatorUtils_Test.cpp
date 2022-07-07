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

// Still need to change name to "SymbolicOperator"

#include <gtest/gtest.h>

#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility> // pair
#include <vector>
#include <armadillo>


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

  // Empty op
  SymbolicOperator op;
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ( qwcmat , Vec2DMat({}) );

  // Identity op
  SymbolicOperator op_id;
  op_id.addIdentTerm(42.0);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op_id);
  EXPECT_EQ( qwcmat , Vec2DMat({{0}}) );
  
  // Op with single pauli string of length-1
  op.removeAllTerms();
  inp = {{0,'Z'}};
  op.addTerm( inp , 1.);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ( qwcmat , Vec2DMat({{0}}) );

  // Op with single longer pauli string
  op.removeAllTerms();
  inp = {{0,'Z'},{2,'X'},{3,'X'}};
  op.addTerm( inp , 1.);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ( qwcmat , Vec2DMat({{0}}) );

  // op = 2*Z0 + 3*Z1 (another sum of Paulis)
  op.removeAllTerms();
  inp = {{0, 'Z'}};
  op.addTerm(inp, 2);
  inp = {{1, 'Z'}};
  op.addTerm(inp, 3);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ( qwcmat , Vec2DMat({{0,0},{0,0}}) );

  // op = 2*Z0 + 3*Y0 (another sum of Paulis)
  op.removeAllTerms();
  inp = {{0, 'Z'}};
  op.addTerm(inp, 2);
  inp = {{0, 'Y'}};
  op.addTerm(inp, 3);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ( qwcmat , Vec2DMat({{0,1},{1,0}}) );

  // op_4 = 2*Z0 + 3*Y0 + 3*Z1   (another sum of Paulis)
  //SymbolicOperator op_5;
  op.removeAllTerms();
  inp = {{0, 'Z'}};
  op.addTerm(inp, 2);
  inp = {{0, 'Y'}};
  op.addTerm(inp, 3);
  inp = {{1, 'Z'}};
  op.addTerm(inp, 3);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ( qwcmat , Vec2DMat({{0,1,0},{1,0,0},{0,0,0}}) );

  // op_4 = 2*Z0 + 3*Y0 + 3*Z1   (another sum of Paulis)
  // Changing order in construction (ought to be same result)
  //SymbolicOperator op_6;
  op.removeAllTerms();
  inp = {{1, 'Z'}};
  op.addTerm(inp, 3);
  inp = {{0, 'Z'}};
  op.addTerm(inp, 2);
  inp = {{0, 'Y'}};
  op.addTerm(inp, 3);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
  EXPECT_EQ( qwcmat , Vec2DMat({{0,1,0},{1,0,0},{0,0,0}}) );

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
  inp_char = {"X3","X4"};
  op.addTerm(inp_char,1.0);
  inp_char = {"Z1","Z2"};
  op.addTerm(inp_char,1.0);
  inp_char = {"Y1","X3","X4"};
  op.addTerm(inp_char,1.0);
  inp_char = {"Z1","Z2","Z3"};
  op.addTerm(inp_char,1.0);
  inp_char = {"Z1","Z2","Z3","Z4"};
  op.addTerm(inp_char,1.0);
  inp_char = {"Z1"};
  op.addTerm(inp_char,1.0);
  inp_char = {"Y1","Y2","X3","X4"};
  op.addTerm(inp_char,1.0);
  qwcmat = SymbolicOperatorUtils::qubitwiseCommutation(op);
//  cout << qwcmat << endl;
  Vec2DMat goldMat = 
    {{0,0,1,1,1,1,0},
     {0,0,1,1,1,1,0},
     {1,1,0,0,0,0,0},
     {1,1,0,0,0,0,0},
     {1,1,0,0,0,0,1},
     {1,1,0,0,0,0,1},
     {0,0,0,0,1,1,0}};
  EXPECT_EQ( qwcmat , goldMat );
  
  vector<int> coloring = SymbolicOperatorUtils::getGroupsQWC(qwcmat);
  vector<int> goldColoring = {0,0,1,1,1,1,0};
  EXPECT_EQ( coloring, goldColoring );
  /*
   * Would next use op.getOrderedPStringList() to be able to relate rows/cols to pstrings
   * */
  

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
  EXPECT_THROW(SymbolicOperatorUtils::findFirstPauliStringBasis(inp_empty), std::logic_error);

  // Incorrect Pauli string
  pstring inp_incorrect = {{0, 'A'}};
  EXPECT_THROW(SymbolicOperatorUtils::findFirstPauliStringBasis(inp_incorrect), std::invalid_argument);

  std::vector<double> variable_params;
  std::vector<double> e_variable_params{0.0, 0.0, FP_PI, FP_PIby2};
  int num_qbits = 2;

  SymbolicOperatorUtils::applyBasisChange(inp_Z0Y1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);

  variable_params = {};
  e_variable_params = {FP_PIby2, FP_PI, FP_PI, FP_PIby2};
  pstring inp_X0Y1 = {{0, 'X'}, {1, 'Y'}};
  SymbolicOperatorUtils::applyBasisChange(inp_X0Y1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);

  variable_params = {};
  e_variable_params = {FP_PI, FP_PIby2, FP_PIby2, FP_PI};
  pstring inp_Y0X1 = {{0, 'Y'}, {1, 'X'}};
  SymbolicOperatorUtils::applyBasisChange(inp_Y0X1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);

  variable_params = {};
  e_variable_params = {FP_PIby2, FP_PI, FP_PIby2, FP_PI};
  pstring inp_X1 = {{1, 'X'}};
  SymbolicOperatorUtils::applyBasisChange(inp_X1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);

  variable_params = {};
  e_variable_params = {FP_PI, FP_PIby2, FP_PI, FP_PIby2};
  pstring inp_Y1 = {{1, 'Y'}};
  SymbolicOperatorUtils::applyBasisChange(inp_Y1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);

  variable_params = {};
  e_variable_params = {0.0, 0.0, 0.0, 0.0};
  pstring inp_Z1 = {{1, 'Z'}};
  SymbolicOperatorUtils::applyBasisChange(inp_Z1, variable_params, num_qbits);
  EXPECT_EQ(variable_params, e_variable_params);

}





