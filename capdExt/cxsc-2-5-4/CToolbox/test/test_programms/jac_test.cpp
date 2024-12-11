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
// Example: Automatic differentiation for gradients
// This is an implementation of Newton's Method for computing a zero of
// a continuously differentiable multi-dimensional function.
//
//    given:      the function f(x) and the starting value x[0]
//    iteration:  x[n+1] := x[n] - InvJf(x[n])*f(x[n]) , n = 0,1,...
//
// where InvJf(x) denotes the inverse of the Jacobian matrix of f(x).
// All real operations are replaced by interval operations, all function and
// derivative evaluations are calculated by differentiation arithmetic.
//----------------------------------------------------------------------------
#include <matinv.hpp>       // Matrix inversion
#include <grad_ari.hpp>     // Gradient differentiation arithmetic

#include "comparator.hpp"

using namespace cxsc;
using namespace std;

const int
  fDim = 2,
  nmax = 100;

ivector  fx(fDim), x(fDim);
imatrix  Jfx(fDim,fDim);
rmatrix  InvJfx(fDim,fDim);
int      n, Error;

GTvector f ( const GTvector& x )
{
  GradType x1sqr(fDim);
  GTvector result(fDim);

  x1sqr = sqr(x[1]);
  result[1] = ((6.0*x1sqr - interval(252.0)/10.0)*x1sqr + 24.0)*x[1]
              - 6.0*x[2];
  result[2] = 12.0*x[2] - 6.0*x[1];
  return result;
}

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

  output << SetPrecision(23,15) << Scientific;   // Output format

  output << "Newton's method for finding a zero of Hansen's function:" << endl;
  output << "f1(x) = 6 x1^5 - 25.2 x1^3 + 24 x1 - 6 x2" << endl;
  output << "f2(x) = 12 x2 - 6 x1" << endl << endl;
  output << "Starting vector x: " << endl;
  input  >> x;
  output << "Iteration" << endl;

  fJEvalJ(f, x, fx, Jfx);

  n = 0;
  do {
    n++;
    output << "x:" << endl << x << "f(x):" << endl << fx << endl;
    MatInv(mid(Jfx), InvJfx, Error);
    if (!Error) {
      x = x - InvJfx * fx;
      fJEvalJ(f, x, fx, Jfx);
    }
  } while ( !( (in(0.0,fx[1]) && in(0.0,fx[2])) || (n >= nmax) || Error ) );
  output << endl;

  if (!Error) {
    output << "Zero:" << endl << x;
    output << "Function value:" << endl << fx << endl;
    output << "Expected zeros:" << endl;
    output << "-1.7475... or -1.0705... or 0.0 or 1.0705... or 1.7475..."
         << endl;
    output << "-0.8737...    -0.5352...    0.0    0.5352... or 0.8737..."
         << endl;
  }
  else
    output << MatInvErrMsg(Error) << endl;

    ComparatorCTests cmp;
    cmp.debugmode(false);
    cmp.compare(output,compinput);    
    result = cmp.getResult(expected);
    //cout << "return code = " << result << endl;     
  }
  return result;
}
