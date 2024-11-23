// Trace of a (complex) matrix product  
// Let C denote the matrix product A*B. 
// Then the diagonal entries of C are added to get the trace of C. 

#include <iostream>
#include "cmatrix.hpp"       // Use the complex matrix package
using namespace std;
using namespace cxsc;
 
int main()
{
   int n;
   cout << "Please enter the matrix dimension n: ";  cin  >> n;
   cmatrix A(n,n), B(n,n);   // Dynamic allocation of A, B
   cdotprecision accu;       // Allows exact computation of dotproducts
   cout << "Please enter the matrix A:" << endl;  cin  >> A;
   cout << "Please enter the matrix B:" << endl;  cin  >> B;
   accu = 0.0;               // Clear accumulator
   for (int i=1; i<=n; i++) accumulate(accu, A[i], B[Col(i)]);
   // A[i] and B[Col(i)] are subarrays of type rvector.

   // The exact result stored in the complex dotprecision variable accu 
   // is rounded to the nearest complex floating point number:
   complex result = rnd(accu);  
   cout << SetPrecision(12,6) << RndNext << Dec;
   cout << "Trace of product matrix: " << result << endl;
   return 0;
}

/* --------------------------- Output ------------------------------
Please enter the matrix dimension n: 3
Please enter the matrix A:
(1,0) (2,0) (3,0)
(4,0) (5,0) (6,0)
(7,0) (8,0) (9,0)
Please enter the matrix B:
(10,0) (11,0) (12,0)
(13,0) (14,0) (15,0)
(16,0) (17,0) (18,0)
Trace of product matrix: (  666.000000,    0.000000) 
------------------------------------------------------------------*/
