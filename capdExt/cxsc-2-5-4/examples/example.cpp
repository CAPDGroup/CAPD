#include "interval.hpp"  // predefined interval arithmetic
#include <iostream>
using namespace cxsc;
using namespace std;

int main()
{
  interval a, b;            // Standard intervals     
  a = 1.0;                  // a   = [1.0,1.0]       
  b = 3.0;                  // b   = [3.0,3.0]        
  cout << "a/b = " << a/b << endl;
  return 0;
}

/* --------------------------- Output ------------------------------
a/b = [  0.333333,  0.333334]
------------------------------------------------------------------*/
