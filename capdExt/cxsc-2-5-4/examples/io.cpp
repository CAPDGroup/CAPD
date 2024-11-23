#include <iostream>
#include "interval.hpp"
using namespace cxsc;
using namespace std;

int main()
{
  real a, b;

  cout << "Please enter real a: ";
  cout << RndDown;             // set rounding mode 
  cin  >> a;                   // read a rounded downwards                  
  cout << SetPrecision(7,4);   // set field width and number of digits 
                               // for output
  cout << a << endl;           // print a rounded downwards to 4 digits  
  cout << RndUp;              
  cout << a << endl;           // print a rounded upwards to 4 digits 
    
  "0.3" >> b;                  // convert the string "0.3" to a floating 
                               // point value b using rounding down    

  cout << SetPrecision(18,15); // from now on print 15 digits

  cout << b << endl;           // print b rounded upwards to 15 digits 
  cout << RndDown;
  cout << b << endl;           // print b rounded downwards to 15 digits 
  interval x;
  "[1.5, 2]" >> x;             // string to interval conversion 
  cout << x << endl;           // print interval x using default setting 

  cout << SaveOpt;             // push I/O parameters to internal stack 
  cout << SetPrecision(10,7);  // modifies output field width and 
                               // number of digits to be print 
  cout << x << endl;           // print x in the modified format 
  cout << RestoreOpt;          // pop parameters from internal stack   
  cout << x << endl;           // again, print x using the former format 
  return 0;
}

/* --------------------------- Output ------------------------------
Please enter real a: 0.3
 0.2999
 0.3000
 0.300000000000001
 0.300000000000000
[ 1.500000000000000, 2.000000000000000]
[ 1.5000000, 2.0000000]
[ 1.500000000000000, 2.000000000000000]
------------------------------------------------------------------*/
