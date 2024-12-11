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
// Example: Global optimization (one-dimensional)
// This program uses module 'gop1' to compute the global optimizers of the
// functions
//
//    f(x) = (x + sin(x)) * exp(-sqr(x)).
//
// and
//              5
//    g(x) = - sum ( k * sin((k+1)*x + k) )
//             k=1
//
// A starting interval and a tolerance must be entered.
//----------------------------------------------------------------------------
#include <gop1.hpp>     // Global optimization

#include "comparator.hpp"

using namespace cxsc;
using namespace std;

stringstream output;

DerivType f ( const DerivType &x )
  { return (x + sin(x)) * exp(-sqr(x)); }

DerivType g ( const DerivType &x )
{
  DerivType  s;
  int        k;

  s = DerivConst(0.0);
  for (k = 1; k <= 5; k++)
    s = s + double(k) * sin ( double(k+1) * x + double(k) );
  return -s;
}

//---------------------------------------------------------------------------
// Function for printing and reading information to call the function
// 'AllGOp1'. It must be called with the function parameter 'f' and a
// string 'Name' containing a textual description of that function.
//---------------------------------------------------------------------------
void compute ( ddf_FctPtr f, char *Name, istream &input )
{
  interval   SearchInterval, Minimum;
  real       Tolerance;
  ivector    Opti;
  intvector  Unique;
  int        NumberOfOptis, i, Error;

  output << SetPrecision(23,15) << Scientific;   // Output format

  output << "Computing all global minimizers of the function  " << Name << endl;
  output << "Search interval      : "; input >> SearchInterval;
  output << "Tolerance (relative) : "; input >> Tolerance;
  output << endl;

  AllGOp1(f,SearchInterval,Tolerance,
            Opti,Unique,NumberOfOptis,Minimum,Error);

  for (i = 1; i <= NumberOfOptis; i++) {
    output << Opti[i] << endl;
    if (Unique[i])
      output << "encloses a locally unique candidate for a global minimizer!";
    else
      output << "may contain a local or global minimizer!";
    output << endl;
  }

  if (NumberOfOptis != 0)
    output << endl << Minimum << endl
         << "encloses the global minimum value!" << endl;
  output << endl << NumberOfOptis << " interval enclosure(s)" << endl;

  if (Error)
    output << endl << AllGOp1ErrMsg(Error) << endl;
  else if ( (NumberOfOptis == 1) && Unique[1] )
    output << endl << "We have validated that there is "
                    "a unique global optimizer!" << endl;
}
// ARGV-VERSION
int main(int argc, char *argv[]) { 
  int result = 1;
  if (argc > 3) {
    ifstream input(argv[1]); 
    ifstream compinput(argv[2]);   
    ifstream expected(argv[3]);
    
    compute(f,(char*)"(x + SIN(x))*EXP(-x^2)",input);
    output << endl << endl;
    compute(g,(char*)"-SUM(k*SIN((k+1)*x+k),k,1,5)",input);
  
    ComparatorCTests cmp;
    cmp.debugmode(false);
    cmp.compare(output,compinput);   
    result = cmp.getResult(expected);
     //cout << "return code = " << result << endl; 
  }
  return result;
}
