/////////////////////////////////////////////////////////////////////////////
//
/// @file LorenzPeriodicOrbit.cpp
///
/// @author Daniel Wilczak
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

/*
 * This example reproduces a computer assisted proof of the existence od 116 periodic orbits for the Lorenz system.
 * The original proof has been given in the paper
 *  Z. Galias, W. Tucker,
 *  Validated study of the existence of short cycles for chaotic systems using symbolic dynamics and interval tools.
 *  Int. J. Bifurcation and Chaos, 21(2):551-563, 2011.
 *
 * This example shows also how to speed up computations by means of usage of static memory allocation.
 * This proof executes within less than 6 seconds on a laptop-type computer.
 *
 * Allocating objects on storage (malloc or new) can be significant overhead.
 * If the dimension is fixed and known at compile time one can use static allocation on the stack.
 * The CAPD library provides a special header file
 *    capd/fdcapdlib.h
 * that can be used instead of generic
 *    capd/capdlib.h
 * Before including this file you MUST define two macros
 * 1. CAPD_USER_NAMESPACE
 *    - is the namespace in which most important types (like vectors, matrices, ODE solvers, Poincare maps)
 *      will be defined for you.
 * 2. CAPD_DEFAULT_DIMENSION
 *    - is a nonnegative number that stands for the dimension.
 *
 * NOTE: if you need few different dimensions in the same program you can undef these macros and define them again.
 *       For example:
 *
 * #define CAPD_USER_NAMESPACE capd3
 * #define CAPD_DEFAULT_DIMENSION 3
 * #include "capd/fdcapdlib.h
 * #undef CAPD_USER_NAMESPACE
 * #undef CAPD_DEFAULT_DIMENSION
 *
 * #define CAPD_USER_NAMESPACE capd5
 * #define CAPD_DEFAULT_DIMENSION 5
 * #include "capd/fdcapdlib.h
 * #undef CAPD_USER_NAMESPACE
 * #undef CAPD_DEFAULT_DIMENSION
 *
 * Now capd3::IVector and capd5::IVector are 3 and 5 dimensional vectors, respectively.
 * The same convention apply to matrices, ODE solvers, etc.
 *
 * In this example we will use 3-dimensional vectors when integration of the Lorenz system
 * and variable size vectors and matrices when computing Interval Newton Operator
 * (this is not time-critical operation in this example).
 */

#include <iostream>
#include <sstream>

// This file defines all types in default namespace "capd".
// Objects are of arbitrary dimension.
#include "capd/capdlib.h"

// Here we define 3-dimensional objects
#define CAPD_DEFAULT_DIMENSION 3
#define CAPD_USER_NAMESPACE capd3
#include "capd/fdcapdlib.h"
#undef CAPD_DEFAULT_DIMENSION
#undef CAPD_USER_NAMESPACE

// data for periodic orbits
#include "LorenzPeriodicOrbits.dat"

using namespace capd;
using namespace std;

/*
  * The following function computes the Interval Newton Operator for the map
  * 
  * F: (x0,x1,x2,...x_{n-1}) -> (x0-P(x_{n-1}),x1-P(x0), .... , x_{n-1} - P(x_{n-2}))
  * 
  * where P is the Poincare map for the Lorenz system. It verifies the existence 
  * and uniqueness of a periodic point in a given set X.
  *
  * Parameters are:
  * @param pm - an instance of PoincareMap for the Lorenz system
  * @param X[] - array of intervals that contains our candidate for periodic orbits.
  * @param period - period of the point with respect to (half) Poincare map
  * 
  * @returns - true if the orbit has been validated 
 */

bool provePeriodicOrbit(capd3::IPoincareMap& pm, interval X[], int dim)
{
  IVector imCenter(dim), Y(dim,X);
  IMatrix D(dim,dim);
  IVector center = midVector(Y);

  for(int i=0;i<dim;i+=2){
    int prev = (dim+i-2)%dim;
		
    // First we compute image at center, i.e center_i - P(center_{i-1}).
    // Note that the points on Poincare section have the form (x,y,z=27).
    capd3::C0HORect2Set s1({center[prev],center[prev+1],interval(27.)});
    capd3::IVector x = pm(s1);
    imCenter[i] = center[i] - x[0];
    imCenter[i+1] = center[i+1] - x[1];
		
    // Then we compute derivative at the set X
    capd3::C1Rect2Set s2({X[prev],X[prev+1],interval(27.)});
    // here we compute monodromy matrix, not derivative of Poincare map
    capd3::IMatrix DP(3,3);
    x = pm(s2,DP);
    // This member function recomputes monodromy matrix into derivatives of Poincare map
    DP = pm.computeDP(x,DP);
    // The matrix DP is 3x3. Take the 2x2 slice that corresponds to (x,y) variables.
    D[i][i] = 1;
    D[i + 1][i + 1] = 1;
    D[i][prev] -= DP[0][0];
    D[i][prev + 1] -= DP[0][1];
    D[i + 1][prev] -= DP[1][0];
    D[i + 1][prev + 1] -= DP[1][1];
  }

  // compute interval Newton operator
  IVector N = center - matrixAlgorithms::krawczykInverse(D) * imCenter;
  // check the inclusion
  return subsetInterior(N, Y);
}

// ----------------------------------- MAIN ----------------------------------------

int main()
{
  int order = 11;
  double tolerance = 1e-7;
  
  try{
    // Here we define vector field, ODE solver, Poincare section and Poicnare map
    // The Poincare section is z=27, i.e. index 2 coordinate of 3 is equal to 27
    capd3::IMap vectorField("par:q;var:x,y,z;fun:10*(y-x),x*(28-z)-y,x*y-q*z;");
    capd3::IOdeSolver solver(vectorField,order);
    capd3::ICoordinateSection section(3,2,27.);
    capd3::IPoincareMap pm(solver,section);
    vectorField.setParameter("q", interval(8.) / interval(3.));
    solver.setAbsoluteTolerance(tolerance);
    solver.setRelativeTolerance(tolerance);

    const int numberOfPeriodicOrbits = 116;
    double xMin, xMax, yMin,yMax, temp;
    int period, validated=0;

    // Array of 40 intervals for storing initial conditions of periodic orbits
    // Maximal period is 20, each point on section has two coordinates (x,y).
    interval X[40];

    for(int i = 0;i<numberOfPeriodicOrbits;++i){
      istringstream in(lorenzPOData[i]);
      in >> period;
      period*=2;
      for(int j=0;j<period;++j){
        in >> xMin >> xMax >> yMin >> yMax >> temp >> temp;
        X[2*j] = interval(xMin,xMax);
        X[2*j+1] = interval(yMin,yMax);
      }
      validated += provePeriodicOrbit(pm,X,2*period);
      cout << "already validated: " << validated << " out of " << (i+1) << endl;
    }
  } catch(exception& e) {
    cout << "\n\nException caught: "<< e.what() << endl;
  }
  return 0;
} // END
