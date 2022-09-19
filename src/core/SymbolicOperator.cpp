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

#include <fstream>
#include <iostream>
#include <istream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace hybrid {

namespace quantum {

namespace core {

using std::string;
using std::vector;
using namespace std;

// Add pauli string (as set of pairs)
void SymbolicOperator::addTerm(pstring &inpp, ComplexDP k,
                               bool check_validity) {

  if (check_validity) {
    // TODO: Add more validation
    // Iterate through set
    pstring::iterator it = inpp.begin();

    // For this constructor, we do not allow >1 op on same qubit id (wouldn't
    // even know what order.) Instead, may use the other constructor
  }

  // Add to map (works properly whether or not key already existed)
  op_sum[inpp] += k;

  // Remove if zero
  if (abs(op_sum[inpp]) < zero_thresh) {
    op_sum.erase(inpp);
  }
}

// Add pauli string (as character string)
// Constructor with vector of strings, e.g. {"X10","Z4"}
// This allows for inputs like [X0 X0 Y2], which addTerm(pstring&) does not
// allow for. Empty vector string mean identity is being added
void SymbolicOperator::addTerm(vector<string> &vecstr, ComplexDP k) {

  // Create new object, that will ensure algebra is correct
  SymbolicOperator newop;
  newop.addIdentTerm(); // Will be multiplied

  for (int i = 0; i < vecstr.size(); i++) {

    string rawstr(vecstr[i]);
    pstring p = {processLocCharString(rawstr)};
    SymbolicOperator locop;
    locop.addTerm(p);

    newop = newop * locop;
  }

  auto iter = newop.op_sum.begin();
  if (iter != newop.op_sum.end()) {
    pstring locstring = newop.op_sum.begin()->first;
    this->op_sum[locstring] += k;
  }
}

// Add identity term with arbitrary coefficient
void SymbolicOperator::addIdentTerm(ComplexDP k) {

  pstring inp_id{};
  this->addTerm(inp_id, k);
}

// Process local character string, i.e. "Z10"
op_pair SymbolicOperator::processLocCharString(std::string inp) {

  // Must be 'X', 'Y', or 'Z'
  char sglP(inp[0]);
  if (sglP != 'X' && sglP != 'Y' && sglP != 'Z') {
    throw std::invalid_argument("Char must be in {X,Y,Z}");
  }

  // Qubit ID (arbitrarily large int)
  int qid = stoi(inp.substr(1, inp.length()));

  // Add local Pauli string to data structure
  return make_pair(qid, sglP);
}

// Addition
SymbolicOperator SymbolicOperator::operator+(const SymbolicOperator &inp) {

  // New object, copy from this
  SymbolicOperator newop;
  newop.op_sum = this->op_sum;

  // Iterate through inp
  for (MapPString::const_iterator it = inp.op_sum.begin();
       it != inp.op_sum.end(); ++it) {

    newop.op_sum[it->first] += it->second;
  }

  return newop;
}

// Multiplication between two SymbolicOperators
SymbolicOperator SymbolicOperator::operator*(const SymbolicOperator &p_right) {

  // New object
  SymbolicOperator newop;

  for (MapPString::iterator Pi = this->op_sum.begin(); Pi != this->op_sum.end();
       ++Pi) {

    for (MapPString::const_iterator Pj = p_right.op_sum.begin();
         Pj != p_right.op_sum.end(); ++Pj) {

      pstring::iterator loc_i = Pi->first.begin();
      pstring::iterator loc_j = Pj->first.begin();

      pstring new_pstring;

      // Used to keep track of the i & -i that accumulate
      ComplexDP quaternion_coeff = {1, 0};

      // For a given pstring-pstring pair:
      while (loc_i != Pi->first.end() || loc_j != Pj->first.end()) {

        // When only one is at end of iterator
        if (loc_i == Pi->first.end()) {

          while (loc_j != Pj->first.end()) {
            int b = loc_j->first;
            char V = loc_j->second;
            new_pstring.insert({b, V});
            ++loc_j;
          }
          break;

        } else if (loc_j == Pj->first.end()) {

          while (loc_i != Pi->first.end()) {
            int a = loc_i->first;
            char U = loc_i->second;
            new_pstring.insert({a, U});
            ++loc_i;
          }
          break;
        }

        int a = loc_i->first;
        int b = loc_j->first;
        char U = loc_i->second;
        char V = loc_j->second;

        if (a == b) {
          // Same qubit id

          if (U == V) {
            // e.g. X2*X2=Id. Cancels. Add no op.
          } else {
            // X*Y = iZ  |  Y*Z = iX  |  Z*X = iY
            if (U == 'X' && V == 'Y') {
              new_pstring.insert({a, 'Z'});
              quaternion_coeff *= ComplexDP{0, 1};
            } else if (U == 'Y' && V == 'Z') {
              new_pstring.insert({a, 'X'});
              quaternion_coeff *= ComplexDP{0, 1};
            } else if (U == 'Z' && V == 'X') {
              new_pstring.insert({a, 'Y'});
              quaternion_coeff *= ComplexDP{0, 1};
            } else if (U == 'Y' && V == 'X') {
              new_pstring.insert({a, 'Z'});
              quaternion_coeff *= ComplexDP{0, -1};
            } else if (U == 'Z' && V == 'Y') {
              new_pstring.insert({a, 'X'});
              quaternion_coeff *= ComplexDP{0, -1};
            } else if (U == 'X' && V == 'Z') {
              new_pstring.insert({a, 'Y'});
              quaternion_coeff *= ComplexDP{0, -1};
            }
          }

          // Update both iterators
          loc_i++;
          loc_j++;
        } else if (a < b) {
          // They are on different qubits. Incorporate only the lower one.
          new_pstring.insert({a, U});

          // Update only first local Pauli
          loc_i++;
        } else if (b < a) {
          // They are on different qubits. Incorporate the lower one.
          new_pstring.insert({b, V});

          // Update only second local Pauli
          loc_j++;
        }
      }

      // Add the new term
      ComplexDP new_k = quaternion_coeff * Pi->second * Pj->second;
      newop.addTerm(new_pstring, new_k);
    }
  }

  return newop;
}

// // Multiplication by a scalar
// SymbolicOperator SymbolicOperator::operator*(const ComplexDP &inp_right) {
//   // New object, copy from this
//   SymbolicOperator newop;
//   newop.op_sum = this->op_sum;

//   // Iterate through inp
//   for (MapPString::const_iterator it = newop.op_sum.begin();
//        it != newop.op_sum.end(); ++it) {

//     newop.op_sum[it->first] *= inp_right * it->second;
//   }

//   return newop;
// }

bool SymbolicOperator::operator==(const SymbolicOperator &rhs) {
  if (this->op_sum != rhs.op_sum)
    return false;
  return true;
}

// Get char string repr
string SymbolicOperator::getCharString() {

  string strP = "";

  if (op_sum.size() == 0) {
    return "0.0";
  }

  for (MapPString::iterator it = op_sum.begin(); it != op_sum.end(); it++) {

    pstring opstring = it->first;
    ComplexDP k = it->second;

    strP += to_string(k.real()) +
            (imag(k) == 0 ? "" : " + " + to_string(imag(k)) + " i") + " [";

    for (pstring::iterator it_loc = opstring.begin(); it_loc != opstring.end();
         it_loc++) {

      strP += " ";
      strP += it_loc->second + to_string(it_loc->first);
    }

    strP += " ]\n";
  }

  return strP;
}

// Get number of terms
int SymbolicOperator::getNumTerms() { return this->op_sum.size(); }

// Get ordered list of pstrings. Needed to interpret e.g. QWC graph.
vector<pstring> SymbolicOperator::getOrderedPStringList() {

  vector<pstring> listPStrings;
  for (auto it = this->op_sum.begin(); it != this->op_sum.end(); ++it) {
    listPStrings.push_back(it->first);
  }

  return listPStrings;
}

// Remove all terms
void SymbolicOperator::removeAllTerms() { this->op_sum.clear(); }

std::string SymbolicOperator::lstrip(std::string s, std::string matches) {
  const size_t b = s.find_first_not_of(matches);
  if (b == std::string::npos)
    return "";
  return s.substr(b, string::npos);
}

std::string SymbolicOperator::rstrip(std::string s, std::string matches) {
  const size_t e = s.find_last_not_of(matches);
  const size_t range = e + 1;
  return s.substr(0, range);
}

std::string SymbolicOperator::stripws(std::string s) {
  const string &ws = " \t\n\r";
  const size_t b = s.find_first_not_of(ws);
  if (b == string::npos)
    return "";
  const size_t e = s.find_last_not_of(ws);
  const size_t range = e - b + 1;
  return s.substr(b, range);
}

std::vector<std::string> SymbolicOperator::split(const std::string &str,
                                                 const char delimiter) {
  std::string s1;

  string elem;
  vector<string> columns;
  stringstream ss(str);
  while (getline(ss, elem, delimiter)) {
    columns.push_back(elem);
  }

  return columns;
}

// Construct Hamiltonian from file
int SymbolicOperator::construct_hamiltonian_from_file(std::string filename) {
  std::ifstream ifs;
  try {
    ifs.open(filename);
    if (!ifs.is_open()) {
      throw runtime_error("invalid file -> " + filename);
    }
  } catch (std::exception &e) {
    std::cout << e.what() << '\n';
    return EXIT_FAILURE;
  }

  std::string line;
  bool first_line = false;
  unsigned long numqbits;
  try {
    while (std::getline(ifs, line)) {

      if (!first_line) {
        first_line = true;
        vector<string> header = split(stripws(line), ' ');
        if (header.size() > 2)
          throw runtime_error("invalid header -> " + line);
        if (header.size() > 2 && header[1] != "qubits")
          throw runtime_error("invalid data -> " + line);
        std::size_t pos{};
        numqbits = stoul(header[0], &pos);
        if (header.size() > 1 && (pos != header[0].size() || numqbits < 1)) {
          throw runtime_error("invalid input for the number of qubits");
        }
        continue;
      }

      pstring inp_id{};

      std::istringstream iss(line);
      std::string s_line = rstrip(stripws(line), ";");
      std::vector<std::string> pauliterm = split(s_line, ':');

      // read pauli term
      if (pauliterm.size() > 1) {
        std::vector<std::string> terms = split(pauliterm[0], ' ');
        for (auto &t : terms) {
          if (t.size() > 2 || (t.size() > 1 && ((int)t[1] - '0') < 0) ||
              (t.size() > 1 && (((int)t[1] - '0') > numqbits - 1)) ||
              (t.size() == 1 && t[0] != 'I') ||
              (t.size() > 1 && t[0] != 'X' && t[0] != 'Y' && t[0] != 'Z' &&
               t[0] != 'I')) {
            throw runtime_error("invalid term -> " + t);
          }

          // handle pauli operators X, Y & Z
          if (t.size() > 1) {
            inp_id.insert({(int)t[1] - '0', t[0]});
          } else {
            // handle identity
            for (auto qno = 0; qno < numqbits; qno++) {
              inp_id.insert({qno, t[0]});
            }
          }
        }

        // Add the coefficient to the pauli term
        this->addTerm(inp_id, stod(stripws(pauliterm[1])));
      }
    }
  } catch (std::exception &e) {
    std::cout << "error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  ifs.close();
  return EXIT_SUCCESS;
}

} // namespace core
} // namespace quantum
} // namespace hybrid
