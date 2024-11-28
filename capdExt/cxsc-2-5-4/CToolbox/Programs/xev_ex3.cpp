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
// Example for the evaluation of an arithmetical expression using function
// 'Eval()'. The exact result of the expression
//
//   f(x,y) = (1682xy^4 + 3x^3 + 29xy^2 - 2x^5 + 832) / 107751
//
// at x = 192119201 and y = 35675640 is 1783.
//----------------------------------------------------------------------------
#include <expreval.hpp>     // Expression evaluation


using namespace cxsc;
using namespace std;

Staggered f ( StaggArray& v )
{
  Staggered x, y;

  x = v[1];
  y = v[2];
  return( (1682.0*x*Power(y,4) + 3.0*Power(x,3) + 29.0*x*Power(y,2) -
           2.0*Power(x,5) + 832.0) / 107751.0 );
}

int main ( )
{
  real     Eps, Approx;
  int      StaggPrec, Err;
  rvector  Arg(2);
  interval Encl;

  cout << SetPrecision(23,15) << Scientific;   // Output format

  cout << "Evaluation of the expression:" << endl
       << "   f(x,y) = (1682xy^4 + 3x^3 + 29xy^2 - 2x^5 + 832) / 107751"
       << endl
       << "at x = 192119201 and y = 35675640. Note, the exact result is"
       << endl
       << "   f(x,y) = 1783." << endl << endl;

  cout << "Enter the arguments:" << endl;
  cout << "   x = " ; cin >> Arg[1];
  cout << "   y = ";  cin >> Arg[2];
  cout << endl;

  cout << "Desired accuracy:" << endl;
  cout << "   Eps = " ; cin >> Eps;
  cout << endl;

  Eval(f,Arg,Eps,Approx,Encl,StaggPrec,Err);

  if (!Err) {
    cout << "Floating-point evaluation:   " << Approx << endl;
    cout << "Interval enclosure:         " << Encl << endl;
    cout << "Defect corrections needed:    " << StaggPrec << endl;
  }
  else
    cout << EvalErrMsg(Err) << endl;
  return 0;
}
