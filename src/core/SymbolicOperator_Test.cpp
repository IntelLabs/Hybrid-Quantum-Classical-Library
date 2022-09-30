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

#include "SymbolicOperator.hpp"
#include "SymbolicOperatorUtils.hpp"

#include <gtest/gtest.h>

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace std;
using namespace hybrid::quantum::core;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

string tmpfilefortest(const string &prefix, const string &dir) {
  string s = dir + "/" + prefix + "XXXXXX";
  char *cs = new char[s.size() + 1];
  cs[s.copy(cs, s.size())] = '\0';

  // create tmpfile and get file descriptor
  int fd = mkstemp(cs);
  if (fd == -1 || close(fd) == -1) {
    std::cout << "Error creating temporary file (" << errno
              << "): " << strerror(errno);
  }

  string filename(cs);

  delete[] cs;
  return filename;
}

TEST(IdentityTests, Complete) {

  /*
   * TESTS
   *
   * X*Y = iZ
   * Y*Z = iX
   * Z*X = iY
   *
   */

  string charstring;

  // Zero
  SymbolicOperator op_zero;
  charstring = op_zero.getCharString();
  EXPECT_EQ(charstring, "0.0");

  // Identity
  SymbolicOperator op_id;
  pstring inp_id{};
  op_id.addTerm(inp_id);
  charstring = op_id.getCharString();
  EXPECT_EQ(charstring, "1.000000 [ ]\n");

  // Y1
  SymbolicOperator op_1;
  pstring inp_y1{{1, 'Y'}};
  op_1.addTerm(inp_y1);
  charstring = op_1.getCharString();
  EXPECT_EQ(charstring, "1.000000 [ Y1 ]\n");

  // 1.5 [ X0 Y2 Z3 ]
  SymbolicOperator op_2;
  pstring inp{{0, 'X'}, {2, 'Y'}, {3, 'Z'}};
  op_2.addTerm(inp, 1.5);
  charstring = op_2.getCharString();
  EXPECT_EQ(charstring, "1.500000 [ X0 Y2 Z3 ]\n");

  // Y1 + 1.5 [ X0 Y2 Z3 ]
  SymbolicOperator op_3;
  inp = {{1, 'Y'}};
  op_3.addTerm(inp, 1);
  inp = {{0, 'X'}, {2, 'Y'}, {3, 'Z'}};
  op_3.addTerm(inp, 1.5);
  charstring = op_3.getCharString();
  EXPECT_EQ(charstring, "1.500000 [ X0 Y2 Z3 ]\n1.000000 [ Y1 ]\n");

  // Same, but constructed through adding
  // Y1 + 1.5 [ X0 Y2 Z3 ]
  SymbolicOperator op;
  op = op_1 + op_2;
  charstring = op.getCharString();
  EXPECT_EQ(charstring, "1.500000 [ X0 Y2 Z3 ]\n1.000000 [ Y1 ]\n");

  // op_3+op_3
  op = op_3 + op_3;
  charstring = op.getCharString();
  EXPECT_EQ(charstring, "3.000000 [ X0 Y2 Z3 ]\n2.000000 [ Y1 ]\n");

  // op_1 * op_2 = 1.5 [ X0 Y1 Y2 Z3 ]
  op = op_1 * op_2;
  charstring = op.getCharString();
  EXPECT_EQ(charstring, "1.500000 [ X0 Y1 Y2 Z3 ]\n");

  // Z0 X2
  SymbolicOperator op_Z0X2;
  inp = {{0, 'Z'}, {2, 'X'}};
  op_Z0X2.addTerm(inp);
  charstring = op_Z0X2.getCharString();
  EXPECT_EQ(charstring, "1.000000 [ Z0 X2 ]\n");

  // X2
  SymbolicOperator op_X2;
  inp = {{2, 'X'}};
  op_X2.addTerm(inp);

  // Y2
  SymbolicOperator op_Y2;
  inp = {{2, 'Y'}};
  op_Y2.addTerm(inp);

  // X2 * (Z0 X2) --> Z0
  op = op_Z0X2 * op_X2;
  charstring = op.getCharString();
  EXPECT_EQ(charstring, "1.000000 [ Z0 ]\n");

  // Y2 * (Z0 X2)
  op = op_Y2 * op_Z0X2;
  charstring = op.getCharString();
  EXPECT_EQ(charstring, "0.000000 + -1.000000 i [ Z0 Z2 ]\n");

  // (Z0 X2) * Y2
  op = op_Z0X2 * op_Y2;
  charstring = op.getCharString();
  EXPECT_EQ(charstring, "0.000000 + 1.000000 i [ Z0 Z2 ]\n");

  // op_4 = 2*Z0 + 3*Z1 (another sum of Paulis)
  SymbolicOperator op_4;
  inp = {{0, 'Z'}};
  op_4.addTerm(inp, 2);
  inp = {{1, 'Z'}};
  op_4.addTerm(inp, 3);

  // op_3 * op_4 -->
  // ( Y1 + 1.5 [ X0 Y2 Z3 ] ) * ( 2 Z0 + 3 Z1 )
  // 2 Z0 Y1 + (i)3 X1 + (-i)3 Y0 Y2 Z3 + 4.5 X0 Z1 Y2 Z3
  op = op_3 * op_4;
  charstring = op.getCharString();
  EXPECT_EQ(charstring, "4.500000 [ X0 Z1 Y2 Z3 ]\n"
                        "0.000000 + -3.000000 i [ Y0 Y2 Z3 ]\n"
                        "2.000000 [ Z0 Y1 ]\n"
                        "0.000000 + 3.000000 i [ X1 ]\n");

  // Make sure it works with e.g. Y2*Y2=Id
  op = op_Y2 * op_Y2;
  charstring = op.getCharString();
  EXPECT_EQ(charstring, "1.000000 [ ]\n");

  // Make sure multiply by scaled Identity works
  // 3*I*op_4 = 3*I*(2*Z0 + 3*Z1) -> 6*Z0 + 9*Z1
  SymbolicOperator op_3I;
  op_3I.addIdentTerm(3.0);
  op = op_3I * op_4;
  charstring = op.getCharString();
  EXPECT_EQ(charstring, "6.000000 [ Z0 ]\n9.000000 [ Z1 ]\n");

  // Ensure *char string* input (overridden addTerm()) parses correctly
  // Z0 Z0 X2 --> X2
  op = SymbolicOperator();
  vector<string> inp_vecstr;
  inp_vecstr = {"Z0", "Z0", "X2"};
  op.addTerm(inp_vecstr);
  charstring = op.getCharString();
  EXPECT_EQ(charstring, "1.000000 [ X2 ]\n");

  // Equality checks
  // Zero
  SymbolicOperator op_zero_1;
  SymbolicOperator op_zero_2;
  ASSERT_TRUE(op_zero_1 == op_zero_2);

  // Identity
  SymbolicOperator op_id_1, op_id_2;
  pstring inp_id_1, inp_id_2;
  op_id_1.addTerm(inp_id_1);
  op_id_2.addTerm(inp_id_2);
  ASSERT_TRUE(op_id_1 == op_id_2);

  // Test getOrderedPStringList
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
  vector<pstring> ordered_pstrings = op.getOrderedPStringList();
  EXPECT_EQ(ordered_pstrings,
            vector<pstring>({{{1, 'Y'}, {2, 'Y'}, {3, 'X'}, {4, 'X'}},
                             {{1, 'Y'}, {3, 'X'}, {4, 'X'}},
                             {{1, 'Z'}},
                             {{1, 'Z'}, {2, 'Z'}},
                             {{1, 'Z'}, {2, 'Z'}, {3, 'Z'}},
                             {{1, 'Z'}, {2, 'Z'}, {3, 'Z'}, {4, 'Z'}},
                             {{3, 'X'}, {4, 'X'}}}));

  // Incorrect Pauli String - A0
  SymbolicOperator op_A0;
  vector<string> v_inp = {"A0"};
  EXPECT_THROW(op_A0.addTerm(v_inp, 2.0), invalid_argument);

  // Z0 * X0
  SymbolicOperator op_Z0;
  pstring inp_Z0 = {{0, 'Z'}};
  op_Z0.addTerm(inp_Z0, 1.0, false);
  SymbolicOperator op_X0;
  pstring inp_X0 = {{0, 'X'}};
  op_X0.addTerm(inp_X0);
  SymbolicOperator op_Z0X0 = op_Z0 * op_X0;
  charstring = op_Z0X0.getCharString();
  EXPECT_EQ(charstring, "0.000000 + 1.000000 i [ Y0 ]\n");
  // InEquality Test
  ASSERT_FALSE(op_Z0 == op_X0);
  EXPECT_EQ(op_Z0X0.getNumTerms(), 1);

  // Z0 * Y0
  SymbolicOperator op_Y0;
  pstring inp_Y0 = {{0, 'Y'}};
  op_Y0.addTerm(inp_Y0, 1.0, false);
  SymbolicOperator op_Z0Y0 = op_Z0 * op_Y0;
  charstring = op_Z0Y0.getCharString();
  EXPECT_EQ(charstring, "0.000000 + -1.000000 i [ X0 ]\n");

  // Z1 * X1
  SymbolicOperator op_Z1;
  pstring inp_Z1 = {{1, 'Z'}};
  op_Z1.addTerm(inp_Z1, 0.0, false);
  SymbolicOperator op_X1;
  pstring inp_X1 = {{1, 'X'}};
  op_X1.addTerm(inp_X1);
  SymbolicOperator op_Z1X1 = op_Z1 * op_X1;
  charstring = op_Z1X1.getCharString();
  EXPECT_EQ(charstring, "0.0");

  /*
  Y1 Y2 X3 X4
  Y1 X3 X4
  Z1
  Z1 Z2
  Z1 Z2 Z3
  Z1 Z2 Z3 Z4
  X3 X4
  */
}

TEST(ExpectationValueTests, Complete) {

  // test 1
  SymbolicOperator so;
  pstring inp_y1{{0, 'Z'}, {1, 'Z'}};
  so.addTerm(inp_y1, 0.5);

  std::vector<double> ProbReg{0.1, 0.2, 0.1, 0.1};
  std::vector<double> expected{-0.050000000000000003};
  std::vector<double> actual;
  for (const auto &pstr : so.getOrderedPStringList()) {
    actual.push_back(
        so.op_sum[pstr].real() *
        SymbolicOperatorUtils::getExpectValSglPauli(pstr, ProbReg, 2));
  }
  ASSERT_TRUE(expected.size() == actual.size());
  for (int i = 0; i < expected.size(); ++i) {
    EXPECT_DOUBLE_EQ(expected[i], actual[i]);
  }

  // test 2
  SymbolicOperator so_1;
  pstring inp_y2{{0, 'X'}};
  so_1.addTerm(inp_y2, 0.5);

  std::vector<double> ProbReg_1{0.1, 0.2, 0.1, 0.1};
  std::vector<double> expected_1{0.050000000000000003};
  std::vector<double> actual_1;
  for (const auto &pstr : so_1.getOrderedPStringList()) {
    actual_1.push_back(
        so_1.op_sum[pstr].real() *
        SymbolicOperatorUtils::getExpectValSglPauli(pstr, ProbReg_1, 2));
  }
  ASSERT_TRUE(expected_1.size() == actual_1.size());
  for (int i = 0; i < expected_1.size(); ++i) {
    EXPECT_DOUBLE_EQ(expected_1[i], actual_1[i]);
  }

  // test 3
  SymbolicOperator so_2;
  pstring inp_y3{{1, 'X'}};
  so_2.addTerm(inp_y3, 0.25);

  std::vector<double> ProbReg_2{0.1, 0.2, 0.1, 0.1};
  std::vector<double> expected_2{-0.025000000000000008};
  std::vector<double> actual_2;
  for (const auto &pstr : so_2.getOrderedPStringList()) {
    actual_2.push_back(
        so_2.op_sum[pstr].real() *
        SymbolicOperatorUtils::getExpectValSglPauli(pstr, ProbReg_2, 2));
  }
  ASSERT_TRUE(expected_2.size() == actual_2.size());
  for (int i = 0; i < expected_2.size(); ++i) {
    EXPECT_DOUBLE_EQ(expected_2[i], actual_2[i]);
  }

  // test 4
  std::vector<pstring> v_pstr{{{0, 'Z'}, {1, 'Z'}}, {{0, 'X'}}, {{1, 'X'}}};

  std::vector<double> ProbReg_3{0.1, 0.2, 0.1, 0.1};
  std::vector<double> expected_3{-0.10000000000000003};
  std::vector<double> actual_3;
  actual_3.push_back(
      SymbolicOperatorUtils::getExpectValSetOfPaulis(v_pstr, ProbReg_3, 2));
  ASSERT_TRUE(expected_3.size() == actual_3.size());
  for (int i = 0; i < expected_3.size(); ++i) {
    EXPECT_DOUBLE_EQ(expected_3[i], actual_3[i]);
  }

  // test 5
  SymbolicOperator symbop;
  symbop.addTerm(inp_y1, 0.5);
  symbop.addTerm(inp_y2, 0.5);
  symbop.addTerm(inp_y3, 0.25);

  std::vector<double> ProbReg_4{0.1, 0.2, 0.1, 0.1};
  std::vector<double> expected_4{-0.025000000000000008};
  std::vector<double> actual_4;
  actual_4.push_back(SymbolicOperatorUtils::getExpectVal(symbop, ProbReg_4, 2));

  ASSERT_TRUE(expected_4.size() == actual_4.size());
  for (int i = 0; i < expected_4.size(); ++i) {
    EXPECT_DOUBLE_EQ(expected_4[i], actual_4[i]);
  }
}

TEST(ConstructHamiltonianFromFileTests, Complete) {
  SymbolicOperator symbop;
  symbop.construct_hamiltonian_from_file(
      "../src/core/bosehubb_gray_d-4_Nx-002.ham");

  std::vector<double> ProbReg{0.1, 0.2, 0.1, 0.1};
  std::vector<double> expected{1.8000000000000007};
  std::vector<double> actual;

  actual.push_back(SymbolicOperatorUtils::getExpectVal(symbop, ProbReg, 2));

  ASSERT_TRUE(expected.size() == actual.size());
  for (int i = 0; i < expected.size(); ++i) {
    EXPECT_DOUBLE_EQ(expected[i], actual[i]);
  }

  // test 2
  string s1 = "2 qubits\n"
              "  Z0 Z1 : 0.5;  \n"
              "X0 : 0.5;\n"
              "X1 : 0.25;\n";

  string ham_file = tmpfilefortest("ham_file_for_test", "/tmp");
  std::ofstream out(ham_file);
  out << s1;
  out.close();

  ASSERT_TRUE(symbop.construct_hamiltonian_from_file(ham_file) == EXIT_SUCCESS);

  expected.clear();
  expected.push_back(1.7750000000000004);
  actual.clear();
  actual.push_back(SymbolicOperatorUtils::getExpectVal(symbop, ProbReg, 2));

  ASSERT_TRUE(expected.size() == actual.size());
  for (int i = 0; i < expected.size(); ++i) {
    EXPECT_DOUBLE_EQ(expected[i], actual[i]);
  }

  remove(ham_file.c_str());

  // test 3
  // Z is not assigned a qubit number which makes it invalid
  // ignores empty newline
  string s2 = "2 qubits\n"
              "\n"
              "I : 0.5\n"
              "Z Z : 0.5;\n"
              "X0 : 0.5;\n"
              "X1 : 0.25;\n";

  ham_file = tmpfilefortest("ham_file_for_test", "/tmp");
  std::ofstream out2(ham_file);
  out2 << s2;
  out2.close();

  ASSERT_TRUE(symbop.construct_hamiltonian_from_file(ham_file) == EXIT_FAILURE);

  remove(ham_file.c_str());

  // test 4
  // invalid header with extra text
  string s3 = "2 qubits example\n"
              "Z Z : 0.5;\n"
              "X0 : 0.5;\n"
              "X1 : 0.25;\n";

  ham_file = tmpfilefortest("ham_file_for_test", "/tmp");
  std::ofstream out3(ham_file);
  out3 << s3;
  out3.close();

  ASSERT_TRUE(symbop.construct_hamiltonian_from_file(ham_file) == EXIT_FAILURE);

  remove(ham_file.c_str());

  // test 5
  // invalid header
  string s4 = "2 qubts\n"
              "Z Z : 0.5;\n"
              "X0 : 0.5;\n"
              "X1 : 0.25;\n";

  ham_file = tmpfilefortest("ham_file_for_test", "/tmp");
  std::ofstream out4(ham_file);
  out4 << s4;
  out4.close();

  ASSERT_TRUE(symbop.construct_hamiltonian_from_file(ham_file) == EXIT_FAILURE);

  remove(ham_file.c_str());

  // test 6
  // invalid header
  string s5 = "A qubits\n"
              "Z Z : 0.5;\n"
              "X0 : 0.5;\n"
              "X1 : 0.25;\n";

  ham_file = tmpfilefortest("ham_file_for_test", "/tmp");
  std::ofstream out5(ham_file);
  out5 << s5;
  out5.close();

  ASSERT_TRUE(symbop.construct_hamiltonian_from_file(ham_file) == EXIT_FAILURE);

  remove(ham_file.c_str());

  // test 7
  // invalid header
  string s6 = "0 qubits\n"
              "Z Z : 0.5;\n"
              "X0 : 0.5;\n"
              "X1 : 0.25;\n";

  ham_file = tmpfilefortest("ham_file_for_test", "/tmp");
  std::ofstream out6(ham_file);
  out6 << s6;
  out6.close();

  ASSERT_TRUE(symbop.construct_hamiltonian_from_file(ham_file) == EXIT_FAILURE);

  remove(ham_file.c_str());

  // test 8
  // invalid header
  string s7 = "1 qubits\n"
              "A0 : 0.5;\n";

  ham_file = tmpfilefortest("ham_file_for_test", "/tmp");
  std::ofstream out7(ham_file);
  out7 << s7;
  out7.close();

  ASSERT_TRUE(symbop.construct_hamiltonian_from_file(ham_file) == EXIT_FAILURE);

  remove(ham_file.c_str());

  // test 9
  // invalid filename
  ham_file = "/tmp/non-existent-file";
  ASSERT_TRUE(symbop.construct_hamiltonian_from_file(ham_file) == EXIT_FAILURE);

  // test 10
  // test lstrip functionality
  EXPECT_EQ(symbop.lstrip("  TestString", " "), "TestString");

  // test 11
  // test lstrip functionality
  EXPECT_EQ(symbop.lstrip("|", "|"), "");

}
