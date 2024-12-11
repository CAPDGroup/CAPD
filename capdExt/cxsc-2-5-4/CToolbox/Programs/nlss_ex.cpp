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


int main ( )
{
  int  fDim;

  // Get the actual value of 'fDim'
  //-------------------------------
  cout << "Desired dimension    : ";  cin >> fDim;

  ivector   SearchInterval(fDim);
  real      Tolerance;
  imatrix   Solu;
  intvector Unique;
  int       NumberOfSolus, i, Error;

  cout << endl
       << "Computing all solutions for " << name << " in " << fDim
       << " variables" << endl << endl;
  cout << "Search interval      : ";  cin >> SearchInterval;
  cout << "Tolerance (relative) : ";  cin >> Tolerance;
  cout << endl;

  AllNLSS(f, SearchInterval, Tolerance, Solu, Unique, NumberOfSolus, Error);

  cout << SetPrecision(23,15) << Scientific;   // Output format
  for (i = 1; i <= NumberOfSolus; i++) {
    cout << Solu[i];
    if (Unique[i])
      cout << "encloses a locally unique zero!";
    else
      cout << "may contain a zero (not verified unique)!";
    cout << endl << endl;
  }

  cout << NumberOfSolus << " interval enclosure(s)" << endl;

  if (Error)
    cout << endl << AllNLSSErrMsg(Error) << endl;
  else if ( (NumberOfSolus == 1) && Unique[1] )
    cout << endl << "We have validated that there is a globally unique zero!"
         << endl;

  return 0;
}
