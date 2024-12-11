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
// Example: Evaluation of a real polynomial with maximum accuracy
//----------------------------------------------------------------------------
#include <rpeval.hpp>     // Polynomial evaluation

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

  int      Err, No, n;
  real     t, y;
  interval yy;

  do {
    output << "Enter the degree of the polynomial (>=0): ";
    input  >> n; output << endl;
  } while (n < 0);

  RPolynomial p(n);

  output << SetPrecision(23,15) << Scientific;   // Output format

  output << "Enter the coefficients of p in increasing order:" << endl;
  input  >> p; output << endl;
  output << "Enter the argument t for evaluation:" << endl;
  input  >> t; output << endl;

  RPolyEval(p,t,y,yy,No,Err);

  if (!Err) {
    output << "Polynomial:" << endl << p <<endl;
    output << "Floating-point evaluation of p(t) using Horner's scheme:"
         << endl << " " << y << endl << endl;
    output << "Verified inclusion of p(t):"
         << endl << yy << endl << endl;
    output << "Number of iterations needed: " << No << endl;
  }
  else
    output << RPolyEvalErrMsg(Err) << endl;

    ComparatorCTests cmp;
    cmp.debugmode(false);
    cmp.compare(output,compinput);    
    result = cmp.getResult(expected);
    //cout << "return code = " << result << endl;     
  }
  return result;

}
