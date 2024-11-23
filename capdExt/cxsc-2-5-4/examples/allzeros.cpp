// Compute all zeros of the function 
//
//              (x-1)*(exp(-3*x) - power(sin(x), 3))
//
// This example needs the CToolbox sources additionally!

#include "nlfzero.hpp"     // Nonlinear equations module
#include "stacksz.hpp"     // To increase stack size for some
                           // special C++ compiler
using namespace cxsc;
using namespace std;


DerivType f ( const DerivType& x )                    // Sample function
{ 
   return (x-1)*( exp(-3*x) - power(sin(x),3) ); 
}

// The class DerivType allows the computation of f, f', and f'' using
// automatic differentiation; see "C++ Toolbox for Verified Computing"

int main()
{
  interval  SearchInterval;
  real      Tolerance;
  ivector   Zero; 
  intvector Unique;
  int       NumberOfZeros, i, Error;

  cout << SetPrecision(23,15) << Scientific;   // Output format

  cout << "Search interval     : ";
  cin  >> SearchInterval;
  cout << "Tolerance (relative): ";
  cin  >> Tolerance;
  cout << endl;

  // Call the function 'AllZeros()' from the C++ Toolbox
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
  return 0;
}

/* --------------------------- Output ------------------------------
Search interval      : [-1,15]
Tolerance (relative) : 1e-10

[ 5.885327439818601E-001, 5.885327439818619E-001]
encloses a locally unique zero!
[ 9.999999999999998E-001, 1.000000000000001E+000]
encloses a locally unique zero!
[ 3.096363932404308E+000, 3.096363932416931E+000]
encloses a locally unique zero!
[ 6.285049273371415E+000, 6.285049273396501E+000]
encloses a locally unique zero!
[ 9.424697254738511E+000, 9.424697254738533E+000]
encloses a locally unique zero!
[ 1.256637410119757E+001, 1.256637410231546E+001]
encloses a locally unique zero!

6 interval enclosure(s)
------------------------------------------------------------------*/
