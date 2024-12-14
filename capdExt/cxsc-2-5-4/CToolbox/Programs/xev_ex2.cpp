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
// Example for the evaluation of an arithmetic expression using void function
// 'Eval()'. The second order difference quotient for a real function f(x) is
// defined by Df(x,h) := (f(x-h) - 2f(x) + f(x+h)) / h^2. This quotient is
// used to approximate the second derivative of
//
//   f(x) := 540 (x^4 - 23x^3 + 159x^2 - 2x + 45) / (x^3 + 18x^2 + 501x + 20)
//
// at x = 1. Since f is twice differentiable, Df(1,h) should tend to
// f''(1) = 36 if h tends to zero.
//----------------------------------------------------------------------------
#include <expreval.hpp>     // Expression evaluation


using namespace cxsc;
using namespace std;


Staggered f ( Staggered x )
{
  return( 540.0 * (Power(x,4) - 23.0*Power(x,3) + 159.0*Power(x,2) -
                   2.0*x + 45.0) /
                  (Power(x,3) + 18.0*Power(x,2) + 501.0*x + 20.0) );
}

Staggered Df ( StaggArray& v )
{
  Staggered x, h;

  x = v[1];
  h = v[2];
  return( (f(x-h) - 2.0*f(x) + f(x+h)) / (h*h) );
}

int main ( )
{
  real     Eps, Approx;
  int      StaggPrec, Err;
  rvector  Arg(2);
  interval Encl;

  cout << SetPrecision(23,15) << Scientific;   // Output format

  cout << "Evaluation of the second order difference quotient " << endl
       << "   Df(x,h) = (f(x-h) - 2f(x) + f(x+h)) / h^2" << endl
       << "for the function" << endl
       << "   f(x) = 540(x^4 - 23x^3 + 159x^2 - 2x + 45) / "
       << "(x^3 + 18x^2 + 501x + 20)" << endl
       << "Note: f''(1) = 36." << endl << endl;

  cout << "Enter the arguments:" << endl;
  cout << "   x = " ; cin >> Arg[1];
  cout << "   h = ";  cin >> Arg[2];
  cout << endl;

  cout << "Desired accuracy:" << endl;
  cout << "   Eps = " ; cin >> Eps;
  cout << endl;

  Eval(Df,Arg,Eps,Approx,Encl,StaggPrec,Err);

  if (!Err) {
    cout << "Floating-point evaluation:   " << Approx << endl;
    cout << "Interval enclosure:         " << Encl << endl;
    cout << "Defect corrections needed:    " << StaggPrec << endl;
  }
  else
    cout << EvalErrMsg(Err) << endl;
  return 0;
}
