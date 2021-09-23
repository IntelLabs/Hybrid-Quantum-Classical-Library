#include "SymbolicOperator.hpp"
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

class SymbolicOperatorUtils {

  /*
   *
   * Eventually need:
   * -- Commutation graph
   * -- Qubit-wise commutation graph
   * -- Anti-commutation graph

  to do:
  - get the string addTerm to work. (do *not* do multi-term)
  - get scalar mult to work
  - do assertions
  - do the other tests above

   */


public:

  static Vec2DMat qubitwiseCommutation( const SymbolicOperator &symbop );
  static std::string getCharString_pstring( const pstring& inp_pstring );

}; // end of class SymbolicOperator


