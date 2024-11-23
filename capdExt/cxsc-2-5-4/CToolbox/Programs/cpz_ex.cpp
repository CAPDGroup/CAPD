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
// Example: Finding an enclosure for a root of a complex polynomial
//----------------------------------------------------------------------------
#include <cpzero.hpp>     // Roots (zeros) of a complex polynomial


using namespace cxsc;
using namespace std;


int main ( )
{
  int       Err, n;
  complex   z;
  cinterval zz;

  do {
    cout << "Enter the degree of the polynomial (>=0): ";
    cin  >> n; cout << endl;
  } while (n < 0);

  CPolynomial  p(n);
  CIPolynomial qq(n);

  cout << "Enter the coefficients in increasing order:" << endl;
  cin  >> p; cout << endl;
  cout << "Enter the starting approximation:" << endl;
  cin  >> z; cout << endl;

  CPolyZero(p,z,qq,zz,Err);

  if (!Err) {
    cout << SetPrecision(15,7) << Scientific;   // Output format

    cout << "Polynomial: "
         << endl << p << endl;
    cout << "Zero found in:"
         << endl << zz << endl << endl;
    cout << "The coefficients of the reduced polynomial are:"
         << endl << qq << endl;
  }
  else
    cout << CPolyZeroErrMsg(Err) << endl;

  return 0;
}
