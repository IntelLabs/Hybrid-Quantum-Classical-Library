#include "PauliOp.hpp"

#include <iostream>
#include <stdexcept>

using std::vector;
using std::string;
using namespace std;

// Add pauli string (as set of pairs)
void PauliOp::addTerm( pstring& inpp , ComplexDP k, bool check_validity) {

    if (check_validity) {

        // Iterate through set
        pstring::iterator it = inpp.begin();

        // For this constructor, we do not allow >1 op on same qubit id (wouldn't even know what order.)
        // Instead, may use the other constructor

    }

    // Add to map (works properly whether or not key already existed)
    op_sum[inpp] += k;

    // Remove if zero
    if (abs(op_sum[inpp])<zero_thresh) {
        op_sum.erase(inpp);
    }

}

// Add pauli string (as character string)
// Constructor with vector of strings, e.g. {"X10","Z4"}
// This allows for inputs like [X0 X0 Y2], which addTerm(pstring&) does not allow for.
// Empty vector string mean identity is being added
void PauliOp::addTerm( vector<string>& vecstr, ComplexDP k ) {

    // Create new object, that will ensure algebra is correct
    PauliOp newop;
    newop.addIdentTerm(); // Will be multiplied

    for( int i=0; i<vecstr.size(); i++ ) {

        string rawstr = vecstr[i];
        pstring p = {processLocCharString(rawstr)};
        PauliOp locop;
        locop.addTerm(p);

        newop = newop * locop;

    }

    pstring locstring = newop.op_sum.begin()->first;
    this->op_sum[ locstring ] += k;

}

// Add identity term with arbitrary coefficient
void PauliOp::addIdentTerm(ComplexDP k) {

    pstring inp_id{};
    this->addTerm(inp_id,k);

}

// Process local character string, i.e. "Z10"
op_pair PauliOp::processLocCharString( std::string inp ) {

    // Must be 'X', 'Y', or 'Z'
    char sglP = inp[0];
    if(sglP!='X' && sglP!='Y' && sglP!='Z') {
        throw std::invalid_argument("Char must be in {X,Y,Z}");
    }

    // Qubit ID (arbitrarily large int)
    int qid = stoi( inp.substr(1,inp.length()) );

    // Add local Pauli string to data structure
    return make_pair(qid,sglP);

}

// Addition
PauliOp PauliOp::operator+(const PauliOp& inp) {

    // New object, copy from this
    PauliOp newop;
    newop.op_sum = this->op_sum;

    // Iterate through inp
    for(MapPString::const_iterator it = inp.op_sum.begin(); it != inp.op_sum.end(); ++it) {

        newop.op_sum[ it->first ] += it->second;

    }

    return newop;

}

// Multiplication between two PauliOps
PauliOp PauliOp::operator*(const PauliOp& p_right) {

    // New object
    PauliOp newop;

    for(MapPString::iterator Pi = this->op_sum.begin(); Pi != this->op_sum.end(); ++Pi) {

        for(MapPString::const_iterator Pj = p_right.op_sum.begin(); Pj != p_right.op_sum.end(); ++Pj) {

            pstring::iterator loc_i = Pi->first.begin();
            pstring::iterator loc_j = Pj->first.begin();

            pstring new_pstring;

            // Used to keep track of the i & -i that accumulate
            ComplexDP quaternion_coeff = {1,0};

            // For a given pstring-pstring pair:
            while( loc_i != Pi->first.end() || loc_j != Pj->first.end() ) {

                // When only one is at end of iterator
                if ( loc_i == Pi->first.end() ) {

                    while ( loc_j != Pj->first.end() ) {
                        int b = loc_j->first;
                        char V = loc_j->second;
                        new_pstring.insert( {b,V} );
                        ++loc_j;
                    }
                    break;

                } else if ( loc_j == Pj->first.end() ) {

                    while ( loc_i != Pi->first.end() ) {
                        int a = loc_i->first;
                        char U = loc_i->second;
                        new_pstring.insert( {a,U} );
                        ++loc_i;
                    }
                    break;

                }

                int a = loc_i->first;
                int b = loc_j->first;
                char U = loc_i->second;
                char V = loc_j->second;

                if (a == b) {
                    //Same qubit id

                    if(U == V) {
                        // e.g. X2*X2=Id. Cancels. Add no op.
                    } else {
                        // X*Y = iZ  |  Y*Z = iX  |  Z*X = iY
                        if( U == 'X' && V == 'Y' ) {
                            new_pstring.insert( { a, 'Z' } );
                            quaternion_coeff *= ComplexDP{0,1};
                        } else if ( U == 'Y' && V == 'Z' ) {
                            new_pstring.insert( { a, 'X' } );
                            quaternion_coeff *= ComplexDP{0,1};
                        } else if ( U == 'Z' && V == 'X' ) {
                            new_pstring.insert( { a, 'Y' } );
                            quaternion_coeff *= ComplexDP{0,1};
                        } else if ( U == 'Y' && V == 'X' ) {
                            new_pstring.insert( { a, 'Z' } );
                            quaternion_coeff *= ComplexDP{0,-1};
                        } else if ( U == 'Z' && V == 'Y' ) {
                            new_pstring.insert( { a, 'X' } );
                            quaternion_coeff *= ComplexDP{0,-1};
                        } else if ( U == 'X' && V == 'Z' ) {
                            new_pstring.insert( { a, 'Y' } );
                            quaternion_coeff *= ComplexDP{0,-1};
                        }
                    }

                    // Update both iterators
                    loc_i++;
                    loc_j++;
                } else if (a < b) {
                    // They are on different qubits. Incorporate only the lower one.
                    new_pstring.insert( {a,U} );

                    // Update only first local Pauli
                    loc_i++;
                } else if (b < a) {
                    // They are on different qubits. Incorporate the lower one.
                    new_pstring.insert( {b,V} );

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
// PauliOp operator* (ComplexDP inp_right) {

//     // TO DO

// }

bool PauliOp::operator== (const PauliOp& rhs) {
    if(this->op_sum != rhs.op_sum) return false;
    return true;
}

// Get char string repr
string PauliOp::getCharString() {

    string strP = "";

    if (op_sum.size()==0) {
        return "0.0";
    }

    for (MapPString::iterator it=op_sum.begin(); it!=op_sum.end(); it++) {

        pstring opstring = it->first;
        ComplexDP k = it->second;

        strP += to_string(k.real()) +  ( imag(k)==0 ? "" : " + "+to_string(imag(k))+" i" ) + " [";

        for( pstring::iterator it_loc=opstring.begin(); it_loc!=opstring.end(); it_loc++) {

            strP += " ";
            strP += it_loc->second + to_string(it_loc->first);

        }

        strP += " ]\n";

    }

    return strP;

}

// Two Pauli strings
// If for *every* qubit id, either (at least one local op is I) or (both local ops are same Pauli), then the two pauli strings are QWC.

// X0 Y1 Z2 and I0 Y1 Z2 --> yes
// X0 Y1 Z2 and X0 --> yes
// X0 Y1 Z2 and Z0 I1 I2 --> no

// X0 Y1 and Y0 X1 --> *not* QWC, though they *do* commute.

// X0*Y0 = iZ0
// Y1*X1 = -iZ1
// (X0 Y1)*(Y0 X1) = -(-1) Z0 Z1
bool PauliOp::qubitwiseCommutation(std::vector<std::string>& pauliStrings) {
    PauliOp op_1;
    op_1.addTerm( pauliStrings );
    string charstring = op_1.getCharString();
    std::cout << "qubitWiseCommutation: " << charstring << std::endl;
    for(MapPString::iterator Pi = this->op_sum.begin(); Pi != this->op_sum.end(); ++Pi) {
    for(MapPString::const_iterator Pj = op_1.op_sum.begin(); Pj != op_1.op_sum.end(); ++Pj) {
        
    }
    }
    return true;
}