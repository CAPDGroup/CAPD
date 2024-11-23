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
// Example for the evaluation of an arithmetic expression using function
// 'Eval()'. Evaluate the following expression:
//
//    f(x,y) := 1 / ( y^6 - 3*x*y^5 + 5*x^3*y^3 - 3*x^5*y - x^6 ).
//
// With x and y being successive Fibonacci numbers defined by
//
//    a_0 := 0,  a_1 := 1,  a_i := a_i-1 + a_i-2, i >= 2,
//
// we always get f(a_i,a_i+1) = (-1)^i. For instance a_66 = 27777890035288
// and a_67 := 44945570212853 are such two successive Fibonacci numbers.
//----------------------------------------------------------------------------
#include <expreval.hpp>     // Expression evaluation
#include <stacksz.hpp>      // To increase stack size for some
                            // special C++ compilers
#include "comparator.hpp"
                            
using namespace cxsc;
using namespace std;

Staggered f ( StaggArray& v )
{
  Staggered x, y, z;

  x = v[1];
  y = v[2];

  z = 1.0 / ( Power(y,6) - 3.0*x*Power(y,5) + 5.0*Power(x,3)*Power(y,3) -
              3.0*Power(x,5)*y - Power(x,6) );
  return(z);
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

  output << "Evaluation of 1 / ( y^6 - 3*x*y^5 + 5*x^3*y^3 - 3*x^5*y - x^6 )"
       << endl << endl;

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
