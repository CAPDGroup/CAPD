// Runge-Kutta Method

#include <iostream>
#include "rvector.hpp"             // Include dynamic arrays (real vectors)
using namespace std;
using namespace cxsc;

rvector F(const real x, const rvector& Y) // Function definition
{  
   rvector Z(3);                   // Constructor call
                               
   Z[1] =  Y[2]  * Y[3];           // F is independent of x
   Z[2] = -Y[1]  * Y[3];
   Z[3] = -0.522 * Y[1] * Y[2];
   return Z;
}

void Init(real& x, real& h, rvector& Y)
{                                  // Initialization
   x    = 0.0;   h    = 0.1;
   Y[1] = 0.0;   Y[2] = 1.0;  Y[3] = 1.0;
}

int main(void)
{
   real x, h;                      // Declarations and dynamic
   rvector Y(3), K1, K2, K3, K4;   // memory allocation
                                   // Automatic resize of Ki below 
   Init(x, h, Y);                                          
   for (int i=1; i<=3; i++) {             // 3 Runge-Kutta steps
      K1 = h * F(x, Y);                   // with array result
      K2 = h * F(x + h / 2, Y + K1 / 2);
      K3 = h * F(x + h / 2, Y + K2 / 2);
      K4 = h * F(x + h, Y + K3);
      Y  = Y + (K1 + 2 * K2 + 2 * K3 + K4) / 6;
      x += h;
      cout << SetPrecision(10,6) << Dec;  // I/O modification
      cout << "Step: " << i << ",  "
           << "x =   " << x << endl;
      cout << "Y =   " << endl << Y << endl;
   }
   return 0;
}

/* --------------------------- Output ------------------------------
Step: 1,  x =     0.100000
Y =   
  0.099747
  0.995013
  0.997400

Step: 2,  x =     0.200000
Y =   
  0.197993
  0.980203
  0.989716

Step: 3,  x =     0.300000
Y =   
  0.293320
  0.956014
  0.977286
------------------------------------------------------------------*/
