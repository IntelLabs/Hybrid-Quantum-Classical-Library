#include <complex>
#include <map>
#include <set> // These are ordered
#include <string>
#include <utility> // pair
#include <vector>

using op_pair = std::pair<int, char>;
using pstring = std::set<op_pair>;
using ComplexDP = std::complex<double>;
using MapPString = std::map<pstring, ComplexDP>;
using Vec2DMat = std::vector<std::vector<int>>;

class PauliOp {

  /*
   *
   * To compile: g++ -std=c++0x PauliOp_Test.cpp
   *
   * Eventually need:
   * -- Commutation graph
   * -- Qubit-wise commutation graph
   *

      Methods to add:
      - equals (assignment)
      - equality (comparison)
      - Qubit-wise commutation (QWC) Graph
      - Commutation graph
      - Anti-commutation graph

      External methods to add:
      -

      Done:
      - addition
      - multiplication
      - string-based input

  to do:
  - get the string addTerm to work. (do *not* do multi-term)
  - get scalar mult to work
  - do assertions
  - do the other tests above


   */

public:
  void addTerm(pstring &inpp, ComplexDP k = 1, bool check_validity = true);
  void addTerm(std::vector<std::string> &vecstr, ComplexDP k = 1);
  void addIdentTerm(ComplexDP k = 1);
  Vec2DMat qubitwiseCommutation();
  op_pair processLocCharString(std::string inp);
  std::string getCharString();

  PauliOp() = default;
  PauliOp(const PauliOp &) = default;
  PauliOp &operator=(const PauliOp &inp) {
    this->op_sum = inp.op_sum;
    return *this;
  }
  PauliOp operator+(const PauliOp &inp);
  PauliOp operator*(const PauliOp &p_right);
  PauliOp operator*(const ComplexDP &inp_right);
  bool operator==(const PauliOp &rhs);

  double zero_thresh = 1.0e-11;

  MapPString op_sum; // The data structure
  Vec2DMat adj_matrix;
}; // end of class PauliOp
