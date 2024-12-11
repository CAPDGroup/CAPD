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
// This program uses module 'nlinsys' to compute the solution of a system of
// nonlinear equations given by Hansen. A starting interval vector and a
// tolerance must be entered.
//----------------------------------------------------------------------------
#include <nlinsys.hpp>     // Nonlinear system solver
#include <stacksz.hpp>     // To increase stack size for some
                           // special C++ compilers
#include "comparator.hpp"

using namespace cxsc;
using namespace std;


const char* name = "Hansen\'s Function";

//----------------------------------------------------------------------------
// Definition of Hansen's function in 'fDim' variables.
//----------------------------------------------------------------------------
GTvector f ( const GTvector& x )
{
  int      i, n = x.Dim();
  GradType SqrSum(n);
  GTvector Result(n);

  SqrSum = 0.0;
  for (i = 1; i <= n; ++i)
    SqrSum = SqrSum + sqr(x[i]);
  for (i = 1; i <= n; ++i)
    Result[i] = 6.0*x[i]/10.0 - 2.0 + 49.0*x[i]*SqrSum/100.0;
  return Result;
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
    
  int  fDim;

  // Get the actual value of 'fDim'
  //-------------------------------
  output << "Desired dimension    : ";  input >> fDim;

  ivector   SearchInterval(fDim);
  real      Tolerance;
  imatrix   Solu;
  intvector Unique;
  int       NumberOfSolus, i, Error;

  output << endl
       << "Computing all solutions for " << name << " in " << fDim
       << " variables" << endl << endl;
  output << "Search interval      : ";  input >> SearchInterval;
  output << "Tolerance (relative) : ";  input >> Tolerance;
  output << endl;

  AllNLSS(f, SearchInterval, Tolerance, Solu, Unique, NumberOfSolus, Error);

  output << SetPrecision(23,15) << Scientific;   // Output format
  for (i = 1; i <= NumberOfSolus; i++) {
    output << Solu[i];
    if (Unique[i])
      output << "encloses a locally unique zero!";
    else
      output << "may contain a zero (not verified unique)!";
    output << endl << endl;
  }

  output << NumberOfSolus << " interval enclosure(s)" << endl;

  if (Error)
    output << endl << AllNLSSErrMsg(Error) << endl;
  else if ( (NumberOfSolus == 1) && Unique[1] )
    output << endl << "We have validated that there is a globally unique zero!"
         << endl;

    ComparatorCTests cmp;
    cmp.debugmode(false);
    cmp.compare(output,compinput);    
    result = cmp.getResult(expected);
    //cout << "return code = " << result << endl;     
  }
  return result;

}
