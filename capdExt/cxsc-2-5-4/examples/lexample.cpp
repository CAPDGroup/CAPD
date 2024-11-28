#include "l_interval.hpp"  // interval staggered arithmetic
#include <iostream>
using namespace cxsc;
using namespace std;

int main() 
{
  l_interval a, b;         // Multiple-precision intervals 
  a = 1.0;
  b = 3.0;
  stagprec = 2;            // global integer variable      
  cout << SetDotPrecision(16*stagprec, 16*stagprec-3) << RndNext;
  // I/O for variables of type l_interval is done using
  // the long accumulator (i.e. a dotprecision variable)   

  cout << "a/b = " << a/b << endl;  
  return(0); 
}

/* --------------------------- Output ------------------------------
a/b = [ 0.33333333333333333333333333333, 0.33333333333333333333333333334]
------------------------------------------------------------------*/
