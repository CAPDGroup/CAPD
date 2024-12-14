// Interval Newton method using ordinary interval arithmetic
// Verified computation of a zero of the function f() 

#include <iostream>
#include "interval.hpp"          // Include interval arithmetic package
#include "imath.hpp"             // Include interval standard functions
using namespace std;
using namespace cxsc;

interval f(const real x)
{                                // Function f
  interval y(x);                 // y is a thin interval initialized by x
  return sqrt(y) + (y+1)*cos(y); // Interval arithmetic is used
}

interval deriv(const interval& x)
{                                // Derivative function f'
   return 1/(2*sqrt(x)) + cos(x) - (x+1)*sin(x);
}

bool criter(const interval& x)    // Computing: f(a)*f(b) < 0   and
{                                 //            not 0 in f'([x])?
   return Sup( f(Inf(x))*f(Sup(x)) ) < 0.0  &&  !(0.0 <= deriv(x)); 
}                                 // '<=' means 'element of'
 
int main(void)
{
   interval x, xOld;
   cout << SetPrecision(20,15); // Number of mantissa digits in I/O
   x= interval(2,3);
   cout << "Starting interval is [2,3]" << endl;
   if (criter(x))
   {  // There is exactly one zero of f in the interval x
      do {
         xOld = x;
         cout << "Actual enclosure is " << x << endl;
         x = (mid(x)-f(mid(x))/deriv(x)) & x; // Iteration formula
      } while (x != xOld);   
      cout << "Final enclosure of the zero: " << x << endl;
   } 
   else
      cout << "Criterion not satisfied!" << endl;
   return 0;
}

/* Output:

Starting interval is [2,3]
Actual enclosure is [   2.000000000000000,   3.000000000000000]
Actual enclosure is [   2.000000000000000,   2.218137182953809]
Actual enclosure is [   2.051401462380920,   2.064726329907714]
Actual enclosure is [   2.059037791936965,   2.059053011233253]
Actual enclosure is [   2.059045253413831,   2.059045253416460]
Actual enclosure is [   2.059045253415142,   2.059045253415145]
Final enclosure of the zero: [ 2.059045253415142, 
                               2.059045253415145]

*/
