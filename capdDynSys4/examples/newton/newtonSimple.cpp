
/////////////////////////////////////////////////////////////////////////////
/// @file newtonSimple.cpp
///  
///  This file contains examples how to use Newton and Krawczyk classes
///
///  This file differs from newtontst.cpp only in using facade classes 
///  instead of general templates
/// 
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iostream>

#include "capd/capdlib.h"
#include "capd/newton/Newton.h"
#include "capd/newton/Krawczyk.h"
using namespace capd;
using namespace std;

 
/* ******************************************************************
 *  Example 1: Non-rigorous (ordinary) Newton and Krawczyk Method
 *             Each iteration gives new estimate (hopefully better)
 *             for zero of given map (R^n -> R^n)
*/


// For Henon map we iterate nonrigorous Newton Operator.
//  Iteration is not rigorous because we pass non rigorous map (DMap)
void FloatHenonMap()
{
	// Creating Hennon map and setting its parameters
   cout << "\n HENON map:  H(x,y) = (1 - a*x*x + y, b*x)   where a=1.4,  b=0.3 ";
   DMap henon("par:a,b;var:x,y;fun:1-a*x*x+y,b*x;");
   henon.setParameter("a",1.4);
   henon.setParameter("b",0.3);

   // Setting initial conditions
   DVector x(2);
   x[0]=-1300; x[1] = 1000;
   cout << "\n\n Non rigorous NEWTON OPERATOR "
             << "\n Starting point x0 = " << x;

   // We make 5 iterations of Newton method for zero finding
   for(int i=1; i<=5; ++i)
   {
      x = newton::NewtonOperator(x, x, henon);
      cout << "\n x "<< i << " = " << x;
   }
}

/* ******************************************************************
 *  Example 2: Rigorous Newton and Krawczyk Method
 *             Each iteration gives new set (hopefully smaller)
 *             for zero of given map (R^n -> R^n)
*/



// For Henon map we iterate Krawczyk and Newton Interval Operators.
// In this example Newton Method requires much smaller initial set then Krawczyk to work.
void HenonMap()
{
   cout << "\n\n HENON map:  H(x,y) = (1 - a*x*x + y, b*x)   where a=1.4,  b=0.3 ";

   // We define henon map and set its parameters
   IMap Henon("par:a,b;var:x,y;fun:1-a*x*x+y,b*x;");
   Henon.setParameter("a", interval(1.4));
   Henon.setParameter("b", interval(0.3));

   IVector x0(2), K(2), N(2);

   // We define an initial set for iterations of taking Krawczyk Operator  
   K[0] = interval(-100,110); K[1]=interval(-100,110);

   cout << "\n\n Interval KRAWCZYK OPERATOR "
             << "\n for set X  = " << K;

   for(int i=1; i<=5; ++i)
   {
      x0=midVector(K);
      K = newton::KrawczykOperator(x0, K, Henon);
      cout << "\n iteration "<< i << " = " << K;
   }

   // We define an initial set for iterations of taking Newton Operator  
   N[0]=interval(-.1,.1); N[1] = interval(0.9,1.1);

   cout << "\n\n Interval NEWTON OPERATOR "
        << "\n for set X = " << N;

   for(int i=1; i<=5; ++i)
   {
      x0=midVector(N);
      N = newton::NewtonOperator(x0, N, Henon);
      cout << "\n iteration "<< i << " = " << N;
   }
}


/* ******************************************************************************
 * Example 3: Rigorous proof of the existence of zero of 
 *            a second iteration of the Henon map
 *            
 *            In this example we also define special class that can be passed as 
 *            a parameter to Newton or Krawczyk operator (instaed of Map class).
 *            It is useful e.g. when part of the map definition is 
 *            an integration of ODE or map is a composition of several other maps.
 */
 
/*  Class MapIterator computes n-th iterate
 *            of given map.
 */ 
class MapIterator
{
 public:
 // These internal types need to be defined
  typedef IVector  VectorType;
  typedef IMatrix  MatrixType;
  typedef IMap MapType;

  MapIterator(const MapType &f, int iter = 1): map(f)
  {
     numberOfIteration = iter;
     dim = 2;
  }

  /* *** Needed Operators and Functions **************/
   
  // Value of a map - used in Newton for computation in a point
  VectorType operator()(VectorType x)
  {
     for(int i=0; i<numberOfIteration; ++i)
        x = map(x);
     return x;
  }

  // Value and derivative of a map - used in Krawczyk for computation in a point
  VectorType operator()(VectorType x, MatrixType &dF)
  {
     dF = MatrixType::Identity(dim);
     for(int i=0; i<numberOfIteration; ++i)
     {
        dF = map[x]* dF;
        x = map(x);
     }
     return x;
  }

  // Derivative of a map - used in both Krawczyk and Newton for computation on a set X
  MatrixType operator[](VectorType x)
  {
     MatrixType dF = MatrixType::Identity(dim);
     for(int i=0; i<numberOfIteration; ++i)
     {
        dF = map[x]* dF;
        x = map(x);
     }
     return dF;
  }
  
  int dimension() const
  {
     return dim;
  }
  /* *************************** */ 
private:  
// members needed only by this particular class
  int dim;
  int numberOfIteration;
  MapType map;
};


// Rigorous proof of the existence of zero for second iteration of Henon map
// Instead of IMap class we use defined above class MapIterator.
void HenonProof()
{
    std::cout << "\n\n SECOND ITERATION OF HENON IMap:  \n  H(x,y) = (1 - a*x*x + y, b*x)   where a=1.4,  b=0.3 ";

    IMap Henon("par:a,b;var:x,y;fun:1-a*x*x+y,b*x;");
    Henon.setParameter("a",DInterval(1.4));
    Henon.setParameter("b",DInterval(0.3));

    MapIterator IMap(Henon, 2);

    newton::Krawczyk<MapIterator> henon(IMap);

    cout << "\n\n KRAWCZYK PROOF \n\n ";

    DVector x(2);               // Good quess for zero for second iteration of Henon map
    x[0] = -10./3; x[1]=131./9; // It is also a center of the set.
    double size = 1.e-5;        // Size of the set in which we will search for a zero

  try{
    newton::KrawczykResult code = henon.proof(x, size);
    cout << resultToText(code) ;
    cout << "\n\n afer " << henon.numberOfIterations << " iteration of Krawczyk method"
              << "\n set X " << henon.X
              << "\n Krawczyk operator K " << henon.K << std::endl;
  }
  catch(exception& e)
  {
     cout << e.what();
  }

}

// -------------------------------------------------------------------------

int main(int, char *[])
{
   try
   {
      FloatHenonMap();
      HenonMap();
      HenonProof();
   }catch(exception& e)
   {
      cout << "\nException caught:\n" << e.what() << endl;
      return 1;
   }
   return 0;
}

