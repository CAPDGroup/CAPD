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
// Example: Differentiation arithmetic (one-dimensional)
// This is an implementation of Halley's method for approximating a zero of
// a twice continuously differentiable function.
//
//    given:      the function f(x) and the starting value x[0]
//    iteration:  x[k+1] := x[k] + a[k]/(1+b[k]/2)          |
//    with:       a[k]   := - f(x[k])/f'(x[k])              | k=0,1,2,...
//                b[k]   := a[k]*f''(x[k])/f'(x[k])         |
//
// All real operations are replaced by interval operations. All function and
// derivative evaluations are calculated by differentiation arithmetic.
// This approximate method can NOT guarantee to find and enclose a zero of
// the function 'f' because we have enclosed only roundoff errors, not
// truncation errors.
//----------------------------------------------------------------------------
#include"dot.hpp" // workaround for ioflags!!
#include <ddf_ari.hpp>     // Differentiation arithmetic


using namespace cxsc;
using namespace std;


const int MaxIter = 100;   // Maximum number of iterations

interval  x, a, b, fx, dfx, ddfx;
int       n;
real      start;

DerivType f ( const DerivType &xD )                // Sample function
  { return( exp(xD) * sin(4.0*xD) ); }

int main ( )
{
  start = 1.25;

  cout << SetPrecision(23,15) << Scientific;   // Output format

  cout << "Halley's method for the function  f(x) = exp(x) * sin(4*x)"
       << endl << endl;
  cout << "Starting value = " << start << endl;
  cout << "Iteration:" << endl;

  x = start;
  ddfEval(f,x,fx,dfx,ddfx);
  n = 0;
  do {
    n++;
    cout << "x              : " << x << endl;
    cout << "f(x)           : " << fx << endl;
    a = - fx / dfx;
    b = ddfx / dfx * a;
    x = x + a / (1.0 + 0.5*b);
    ddfEval(f,x,fx,dfx,ddfx);
  } while( !( in(0.0,fx) || (n >= MaxIter) ) );

  cout << endl
       << "Approximate zero : " << x << endl;
  cout << "Function value   : " << fx << endl << endl;
  cout << "Expected zero    : 1.570796326794896619202..." << endl;

  return 0;
}
