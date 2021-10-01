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

using op_pair = std::pair<int, char>;
using pstring = std::set<op_pair>;
using ComplexDP = std::complex<double>;
using MapPString = std::map<pstring, ComplexDP>;
using Vec2DMat = std::vector<std::vector<int>>;

using namespace std;
using namespace arma;

  
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

Vec2DMat SymbolicOperatorUtils::qubitwiseCommutation( const SymbolicOperator &symbop ) {
  
  int M = symbop.op_sum.size();

  Vec2DMat vmat =
      Vec2DMat(M, vector<int>(M, 0));
  int row = 0, col = 0;
  for (auto iIter = symbop.op_sum.begin(); iIter != symbop.op_sum.end();
       ++iIter, ++row) {
    col = row+1;
    

    for (auto jIter = next(iIter); jIter != symbop.op_sum.end();
         ++jIter, ++col) {

      pstring::const_iterator loc_i_end = iIter->first.end();
      pstring::const_iterator loc_j_end = jIter->first.end();
      
      // Loop over the two individual Pauli strings
      for (pstring::const_iterator ipstrIter = iIter->first.begin();
           ipstrIter != loc_i_end ; ++ipstrIter) {
        
        bool flag_non_commuting_found = false;

        int a = ipstrIter->first;
        char U = ipstrIter->second;

        for (pstring::const_iterator jpstrIter = jIter->first.begin();
             jpstrIter != loc_j_end; ++jpstrIter) {

          int b = jpstrIter->first;
          char V = jpstrIter->second;
          //cout << "& " << V << b << endl;
          
          if (U!=V && a==b) {
            // Only one mismatch is requried for the two pstrings to *not* be QWC
            // Add edge for non-commutation
            vmat[row][col] = vmat[col][row] = 1;
            flag_non_commuting_found = true;
            break;
          }
          
        }
        
        if (flag_non_commuting_found) break;

      }
    }
  }


  return vmat;
 

}




string SymbolicOperatorUtils::getCharString_pstring( const pstring& inp_pstring ) {

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


 







