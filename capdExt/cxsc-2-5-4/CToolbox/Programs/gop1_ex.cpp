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


using namespace cxsc;
using namespace std;


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
void compute ( ddf_FctPtr f, char *Name )
{
  interval   SearchInterval, Minimum;
  real       Tolerance;
  ivector    Opti;
  intvector  Unique;
  int        NumberOfOptis, i, Error;

  cout << SetPrecision(23,15) << Scientific;   // Output format

  cout << "Computing all global minimizers of the function  " << Name << endl;
  cout << "Search interval      : "; cin >> SearchInterval;
  cout << "Tolerance (relative) : "; cin >> Tolerance;
  cout << endl;

  AllGOp1(f,SearchInterval,Tolerance,
            Opti,Unique,NumberOfOptis,Minimum,Error);

  for (i = 1; i <= NumberOfOptis; i++) {
    cout << Opti[i] << endl;
    if (Unique[i])
      cout << "encloses a locally unique candidate for a global minimizer!";
    else
      cout << "may contain a local or global minimizer!";
    cout << endl;
  }

  if (NumberOfOptis != 0)
    cout << endl << Minimum << endl
         << "encloses the global minimum value!" << endl;
  cout << endl << NumberOfOptis << " interval enclosure(s)" << endl;

  if (Error)
    cout << endl << AllGOp1ErrMsg(Error) << endl;
  else if ( (NumberOfOptis == 1) && Unique[1] )
    cout << endl << "We have validated that there is "
                    "a unique global optimizer!" << endl;
}

int main ( )
{
  compute(f,(char*)"(x + SIN(x))*EXP(-x^2)");
  cout << endl << endl;
  compute(g,(char*)"-SUM(k*SIN((k+1)*x+k),k,1,5)");
  return 0;
}
