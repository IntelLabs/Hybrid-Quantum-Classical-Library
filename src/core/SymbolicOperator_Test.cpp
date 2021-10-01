#include "SymbolicOperator.hpp"

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

/*
 *
 * g++ -std=c++0x SymbolicOperator_Test.cpp
 *
 */

using namespace std;

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
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
  SymbolicOperator op_izmaylov;
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
  vector<pstring> ordered_pstrings = op.getOrderedPStringList();
  EXPECT_EQ( ordered_pstrings, vector<pstring>({
      {{1,'Y'},{2,'Y'},{3,'X'},{4,'X'}},
      {{1,'Y'},{3,'X'},{4,'X'}},
      {{1,'Z'}},
      {{1,'Z'},{2,'Z'}},
      {{1,'Z'},{2,'Z'},{3,'Z'}},
      {{1,'Z'},{2,'Z'},{3,'Z'},{4,'Z'}},
      {{3,'X'},{4,'X'}} 
      }));
      
      
    

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



