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
// Example: Linear systems of equations
//----------------------------------------------------------------------------
#include <matinv.hpp>     // Matrix inversion
#include <linsys.hpp>     // Linear system solver


using namespace cxsc;
using namespace std;


int main ( )
{
  int  Err, n;
  real Cond;

  cout << SetPrecision(23,15) << Scientific;   // Output format

  do {
    cout << "Enter the dimension of the system: ";
    cin  >> n;  cout << endl;
  } while (n <= 0);

  rmatrix A(n,n), R(n,n);      // Dynamic allocation
  rvector b(n);
  ivector x(n);

  cout << "Enter matrix A:" << endl;
  cin  >> A;  cout << endl;
  cout << "Enter vector b:" << endl;
  cin  >> b;  cout << endl;

  LinSolve(A,b,x,Cond,Err);

  if (!Err) {
    // Compare the result to naive floating-point approximation
    MatInv(A,R,Err);
    cout << "Naive floating-point approximation:" << endl
         << R*b << endl;
    cout << "Verified solution found in:" << endl
         << x << endl;
  }
  else
    cout << LinSolveErrMsg(Err) << endl;

  if (!Err || Err == 4)
    cout << "Condition estimate: " << SetPrecision(9,1)
         << Cond << endl;

  return 0;
}
