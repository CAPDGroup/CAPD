//============================================================================
//
//                              Program/Module
//                                   from
//                 C++ TOOLBOX FOR VERIFIED COMPUTING I
//                         Basic Numerical Problems
//
//      Copyright (c) 1995   Rolf Hammer, Matthias Hocks, Dietmar Ratz
//
// For details on theory, algorithms, and programs, see the book
//
//  R. Hammer, M. Hocks, U. Kulisch, D. Ratz:  C++ Toolbox for
//  Verified Computing I - Basic Numerical Problems. Springer-Verlag,
//  Heidelberg, New York, 1995.
//
//============================================================================
//----------------------------------------------------------------------------
// Example: Finding an enclosure for a root of a complex polynomial
//----------------------------------------------------------------------------
#include <cpzero.hpp>     // Roots (zeros) of a complex polynomial

#include "comparator.hpp"

using namespace cxsc;
using namespace std;

// ARGV-VERSION

int main(int argc, char *argv[]) {  
  int result = 1;
  
  if (argc > 3) {
    //data given from Parameters
    ifstream input(argv[1]); 
    ifstream compinput(argv[2]);  
    ifstream expected(argv[3]);
    
    //generated output to compare with compinput
    stringstream output;
  
    int       Err, n;
    complex   z;
    cinterval zz;

    do {
      output << "Enter the degree of the polynomial (>=0): ";
      input  >> n; output << endl;
    } while (n < 0);

    CPolynomial  p(n);
    CIPolynomial qq(n);

    output << "Enter the coefficients in increasing order:" << endl;
    input  >> p; output << endl;
    output << "Enter the starting approximation:" << endl;
    input  >> z; output << endl;

    CPolyZero(p,z,qq,zz,Err);

    if (!Err) {
      output << SetPrecision(15,7) << Scientific;   // Output format

      output << "Polynomial: "
	<< endl << p << endl;
      output << "Zero found in:"
	<< endl << zz << endl << endl;
      output << "The coefficients of the reduced polynomial are:"
	<< endl << qq << endl;
    }
    else
      output << CPolyZeroErrMsg(Err) << endl;
  
  
    ComparatorCTests cmp;
    cmp.debugmode(false);
    cmp.compare(output,compinput);    
    result = cmp.getResult(expected);
    //cout << "return code = " << result << endl;     
  }
  return result;

}
