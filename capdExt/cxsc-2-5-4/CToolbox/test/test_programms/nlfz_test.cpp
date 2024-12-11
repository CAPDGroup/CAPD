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
#include "comparator.hpp"

using namespace cxsc;
using namespace std;


DerivType f ( const DerivType &x )                    // Sample function
  { return exp(-3.0*x) - power(sin(x),3); }

//----------------------------------------------------------------------------
// Function for prompting and reading information to call the function
// 'AllZeros()'. This function must be called with the function 'f'
// and a string 'Name' containing a textual description of that function.
//----------------------------------------------------------------------------

stringstream output;

void compute ( ddf_FctPtr f, const char *Name, istream &input )
{
  interval  SearchInterval;
  real      Tolerance;
  ivector   Zero;
  intvector Unique;
  int       NumberOfZeros, i, Error;

  output << SetPrecision(23,15) << Scientific;   // Output format

  output << "Computing all zeros of the function  " << Name << endl;
  output << "Search interval      : ";
  input  >> SearchInterval;
  output << "Tolerance (relative) : ";
  input  >> Tolerance;
  output << endl;

  AllZeros(f,SearchInterval,Tolerance,Zero,Unique,NumberOfZeros,Error);

  for ( i = 1; i <= NumberOfZeros; i++) {
    output << Zero[i] << endl;
    if (Unique[i])
      output << "encloses a locally unique zero!" << endl;
    else
      output << "may contain a zero (not verified unique)!" << endl;
  }

  output << endl << NumberOfZeros << " interval enclosure(s)" << endl;

  if (Error) output << endl << AllZerosErrMsg(Error) << endl;
}

//--------------------
// Body of the program
//--------------------

// ARGV-VERSION

int main(int argc, char *argv[]) {  
  int result = 1;
  
  if (argc > 3) {
   //data given from Parameters
  ifstream input(argv[1]); 
  ifstream compinput(argv[2]);  
  ifstream expected(argv[3]);  
  
  compute(f, "EXP(-3x)-POWER(SIN(x),3)", input); 
  
  ComparatorCTests cmp;
  cmp.debugmode(true);
  cmp.compare(output,compinput);    
  result = cmp.getResult(expected);
  cout << "return code = " << result << endl;
  }
   
  return 0;  
}
