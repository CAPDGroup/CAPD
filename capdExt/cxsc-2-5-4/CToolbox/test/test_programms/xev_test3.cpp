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

#include "comparator.hpp"

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

  real     Eps, Approx;
  int      StaggPrec, Err;
  rvector  Arg(2);
  interval Encl;

  output << SetPrecision(23,15) << Scientific;   // Output format

  output << "Evaluation of the expression:" << endl
       << "   f(x,y) = (1682xy^4 + 3x^3 + 29xy^2 - 2x^5 + 832) / 107751"
       << endl
       << "at x = 192119201 and y = 35675640. Note, the exact result is"
       << endl
       << "   f(x,y) = 1783." << endl << endl;

  output << "Enter the arguments:" << endl;
  output << "   x = " ; input >> Arg[1];
  output << "   y = ";  input >> Arg[2];
  output << endl;

  output << "Desired accuracy:" << endl;
  output << "   Eps = " ; input >> Eps;
  output << endl;

  Eval(f,Arg,Eps,Approx,Encl,StaggPrec,Err);

  if (!Err) {
    output << "Floating-point evaluation:   " << Approx << endl;
    output << "Interval enclosure:         " << Encl << endl;
    output << "Defect corrections needed:    " << StaggPrec << endl;
  }
  else
    output << EvalErrMsg(Err) << endl;
    ComparatorCTests cmp;
  
    cmp.debugmode(false);
    cmp.compare(output,compinput);    
    result = cmp.getResult(expected);
    //cout << "return code = " << result << endl;     
  }
  return result;

}
