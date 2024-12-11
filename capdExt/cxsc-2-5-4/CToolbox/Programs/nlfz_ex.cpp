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
// Example: Nonlinear equations in one variable
// This program uses module 'nlfzero' to compute the zeros of the function
//
//              exp(-3*x) - power(sin(x), 3)
//
// A starting interval and a tolerance must be entered.
//----------------------------------------------------------------------------
#include <nlfzero.hpp>     // Nonlinear equations
#include <stacksz.hpp>     // To increase stack size for some
                           // special C++ compilers


using namespace cxsc;
using namespace std;


DerivType f ( const DerivType &x )                    // Sample function
  { return exp(-3.0*x) - power(sin(x),3); }

//----------------------------------------------------------------------------
// Function for prompting and reading information to call the function
// 'AllZeros()'. This function must be called with the function 'f'
// and a string 'Name' containing a textual description of that function.
//----------------------------------------------------------------------------
void compute ( ddf_FctPtr f, const char *Name )
{
  interval  SearchInterval;
  real      Tolerance;
  ivector   Zero;
  intvector Unique;
  int       NumberOfZeros, i, Error;

  cout << SetPrecision(23,15) << Scientific;   // Output format

  cout << "Computing all zeros of the function  " << Name << endl;
  cout << "Search interval      : ";
  cin  >> SearchInterval;
  cout << "Tolerance (relative) : ";
  cin  >> Tolerance;
  cout << endl;

  AllZeros(f,SearchInterval,Tolerance,Zero,Unique,NumberOfZeros,Error);

  for ( i = 1; i <= NumberOfZeros; i++) {
    cout << Zero[i] << endl;
    if (Unique[i])
      cout << "encloses a locally unique zero!" << endl;
    else
      cout << "may contain a zero (not verified unique)!" << endl;
  }

  cout << endl << NumberOfZeros << " interval enclosure(s)" << endl;

  if (Error) cout << endl << AllZerosErrMsg(Error) << endl;
}

//--------------------
// Body of the program
//--------------------
int main ( )
 { compute(f, "EXP(-3x)-POWER(SIN(x),3)"); return 0;}
