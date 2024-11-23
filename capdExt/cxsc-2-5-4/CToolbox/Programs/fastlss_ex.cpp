/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.5.4)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2014 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

//----------------------------------------------------------------------------
// Example: Linear systems of equations
//----------------------------------------------------------------------------
#include <fastlss.hpp>     // Linear system solver
#include <matinv_aprx.hpp>

using namespace cxsc;
using namespace std;


int main ( )
{
  int  Err, n, type;

  cout << SetPrecision(23,15) << Scientific;   // Output format

  do {
    cout << "Please select data type:" << endl;
    cout << "(1) real" << endl;
    cout << "(2) interval" << endl; 
    cout << "(3) complex" << endl;
    cout << "(4) cinterval" << endl;
    cin >> type;
  } while(type < 1 || type > 4);

  do {
    cout << "Enter the dimension of the system: ";
    cin  >> n;  cout << endl;
  } while (n <= 0);

  
  if(type == 1) {

    rmatrix A(n,n), R(n,n);      // Dynamic allocation
    rmatrix b(n,1);
    imatrix x(n,1);

    cout << "Enter matrix A:" << endl;
    cin  >> A;  cout << endl;
    cout << "Enter vector b:" << endl;
    cin  >> b;  cout << endl;

    /* Call to solver
    *
    * Mandatory parameters:
    *  -System matrix (rmatrix)
    *  -Matrix of right hand sides (rmatrix)
    *  -Solution matrix x (imatrix)
    *  -Error variable (int)
    *
    * Optional paramters:
    *  -Dot product precision (int)
    *  -Status message output (bool)
    *  -Solver algorithm (int, 
    *                    possible values: LSS_ONLY_PART_ONE (fast, works with condition numbers of up to about 1e15)
    *                                     LSS_ONLY_PART_TWO (slow, works with badly conditioned systems (up to about 1e30),
    *                                     LSS_BOTH_PARTS (default, try algorithm one, if this fails, try algorithm two) )
    *
    * Examples: 
    * Default call, solving the system with dot product precision k=2, no status messages and using both algorithms
    * lss(A,b,x,Err);
    *
    * Call solving the system using precision k=0 (maximum precision), activating status messages,using only algorithm one
    * lss(A,b,x,Err,0,true,LSS_ONLY_PART_ONE);
    */
    lss(A,b,x,Err);

    if (!Err) {
      // Compare the result to naive floating-point approximation
      MatInvAprx(A,R,Err);
      cout << "Naive floating-point approximation:" << endl
           << R*b << endl;
      cout << "Verified solution found in:" << endl << x << endl;
    }

  } else if(type == 2) {

    imatrix A(n,n);      // Dynamic allocation
    imatrix b(n,1);
    imatrix x(n,1);

    cout << "Enter matrix A:" << endl;
    cin  >> A;  cout << endl;
    cout << "Enter vector b:" << endl;
    cin  >> b;  cout << endl;

    //Call to solver, see above for details
    ilss(A,b,x,Err);

    if(!Err)
      cout << "Verified solution found in:" << endl << x << endl;
  
  } else if(type == 3) {

    cmatrix A(n,n), R(n,n);      // Dynamic allocation
    cmatrix b(n,1);
    cimatrix x(n,1);

    cout << "Enter matrix A:" << endl;
    cin  >> A;  cout << endl;
    cout << "Enter vector b:" << endl;
    cin  >> b;  cout << endl;

    //Call to solver, see above for details
    clss(A,b,x,Err);

    if (!Err) {
      // Compare the result to naive floating-point approximation
      MatInvAprx(A,R,Err);
      cout << "Naive floating-point approximation:" << endl
           << R*b << endl;
      cout << "Verified solution found in:" << endl << x << endl;
    }

  } else if(type == 4) {
  
    cimatrix A(n,n);      // Dynamic allocation
    cimatrix b(n,1);
    cimatrix x(n,1);

    cout << "Enter matrix A:" << endl;
    cin  >> A;  cout << endl;
    cout << "Enter vector b:" << endl;
    cin  >> b;  cout << endl;

    //Call to solver, see above for details
    cilss(A,b,x,Err);

    if(!Err)
      cout << "Verified solution found in:" << endl << x << endl;

  }

  
  if(Err) cout << LinSolveErrMsg(Err) << endl;

  return 0;
}
