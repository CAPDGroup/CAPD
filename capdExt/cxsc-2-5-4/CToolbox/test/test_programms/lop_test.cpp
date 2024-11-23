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
// Example for the use of module 'lop' to solve the linear programming problem
// (A,b,c) with:
//                                                  ( 50 )
//       (   1   0 1 0 0 )       (   50 )           (  9 )
//   A = (   0   1 0 1 0 ),  b = (  200 ), and  c = (  0 )
//       ( 100  18 0 0 1 )       ( 5000 )           (  0 )
//                                                  (  0 )
//----------------------------------------------------------------------------
#include <rev_simp.hpp>     // Revised Simplex algorithm
#include <lop.hpp>          // Linear optimization

#include "comparator.hpp"

using namespace cxsc;
using namespace std;


int main(int argc, char *argv[]) {  
  int result = 1;
  
  if (argc > 3) {
    //data given from Parameters
    ifstream input(argv[1]); 
    ifstream compinput(argv[2]);  
    ifstream expected(argv[3]);
    
    //generated output to compare with compinput
    stringstream output;
  int m, n;

  output << "Enter dimensions of linear "                // Read dimension of
       << "optimization problem (m and n): " << endl;  // optimization problem
  input >> m >> n;

  rmatrix    A(m,n);
  rvector    c(n), b(m), x_start;
  intvector  B_start;
  intmatrix  V;
  imatrix    X;
  real       z_start;
  interval   z;
  int        Error, i, NoOfSolutions;

  // Read optimization problem P = (A,b,c)
  //--------------------------------------
  output << endl << "Enter objective function (c[1]..c[n]): " << endl;
  input  >> c;
  output << endl << "Enter matrix of constraints (A[1,1]..A[m,n]): " << endl;
  input  >> A;
  output << endl << "Enter right-hand side (b[1]..b[m]): " << endl;
  input  >> b;
  output << endl;

  // Call revised simplex function
  //------------------------------
  RevSimplex(A, b, c, x_start, B_start, z_start, Error);

  if (Error)
    { output << RevSimplexErrMsg(Error) << endl; return -1; }

  // Display results of calculation
  //-------------------------------
  output << SetPrecision(23,15) << Scientific;  // output format

  output << "Results of approximation:" << endl
       << "Optimal value :  " << z_start << endl << endl
       << "Initial basic index set : " << B_start << endl
       << "Initial optimal solution : " << endl << x_start << endl;

  // Call verifying linear optimization solver
  //------------------------------------------
  LinOpt(A, b, c, B_start, z, V, X, NoOfSolutions, Error);

  if (Error)
    output << LinOptErrMsg(Error) << endl;
  else {
    // Display results of calculation
    //-------------------------------
    output << "Results of verification:" << endl;
    output << "Optimal value interval :  " << z << endl << endl
         << "List of optimal basic index sets : " << endl;
    for (i = 1; i <= NoOfSolutions; i++)
      output << "B" << i << " = " << V[i];
    output << endl
         << "List of optimal basic feasible solutions : " << endl;
    for (i = 1; i <= NoOfSolutions; i++)
      output << "x of B" << i << " = " << endl << X[i];
    output << endl;
  }
    ComparatorCTests cmp;
    cmp.debugmode(false);
    cmp.compare(output,compinput,"lop_test");    
    result = cmp.getResult(expected);
    //cout << "return code = " << result << endl;     
  }
  return result;

}
