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
// Example: Linear systems of equations
//----------------------------------------------------------------------------
#include <matinv.hpp>     // Matrix inversion
#include <linsys.hpp>     // Linear system solver

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
    
  int  Err, n;
  real Cond;

  output << SetPrecision(23,15) << Scientific;   // Output format

  do {
    output << "Enter the dimension of the system: ";
    input  >> n;  output << endl;
  } while (n <= 0);

  rmatrix A(n,n), R(n,n);      // Dynamic allocation
  rvector b(n);
  ivector x(n);

  output << "Enter matrix A:" << endl;
  input  >> A;  output << endl;
  output << "Enter vector b:" << endl;
  input  >> b;  output << endl;

  LinSolve(A,b,x,Cond,Err);

  if (!Err) {
    // Compare the result to naive floating-point approximation
    MatInv(A,R,Err);
    output << "Naive floating-point approximation:" << endl
         << R*b << endl;
    output << "Verified solution found in:" << endl
         << x << endl;
  }
  else
    output << LinSolveErrMsg(Err) << endl;

  if (!Err || Err == 4)
    output << "Condition estimate: " << SetPrecision(9,1)
         << Cond << endl;

      ComparatorCTests cmp;
    cmp.debugmode(false);
    cmp.compare(output,compinput);    
    result = cmp.getResult(expected);
    //cout << "return code = " << result << endl;     
  }
  return result;
}
