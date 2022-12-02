/**
 *  @file diffInclExample.cpp
 *
 *  Created on: 2009-11-13
 *      Author: kapela
 */

#include <iostream>
#include "capd/capdlib.h"
using namespace capd;

/**
 *  In this example we integrate differential inclusion based on rossler equation
 *     x' in -(y+z) + [e]
 *     y' in x+b*y + [e]
 *     z' in b+z*(x-a) + [e]
 *  where
 *    a=5.7,  b=0.2,
 *    [e] = [-eps, eps]
 *  We use two different method based on component wise estimates and on logarithmic norm
 *  and compare results.
 */
void RosslerExample() {

  // f is an unperturbed vector field
  IMap f("par:a,b;var:x,y,z;fun:-(y+z),x+b*y,b+z*(x-a);");
  f.setParameter("a", DInterval(57.) / DInterval(10.));                 // note that 5.7 is not a representable number
  f.setParameter("b", DInterval(2.) / DInterval(10.));

  // we define a perturbation
  double eps = 1.0e-4;      //  maximal forcing (perturbation)
  IMap perturb("par:e;var:x,y,z;fun:e,e,e;");
  perturb.setParameter("e", DInterval(-eps, eps));

  // We set right hand side of differential inclusion to f + perturb
  IMultiMap rhs(f, perturb);

  double timeStep = 1. / 512.; //  time step for integration
  int order = 10; //  order of Taylor method

  // We set up two differential inclusions (they differ in the way they handle perturbations)
  CWDiffInclSolver cwDiffInclSolver(rhs, order, IMaxNorm());
  LNDiffInclSolver lnDiffInclSolver(rhs, order, IEuclLNorm());
  cwDiffInclSolver.setStep(timeStep);
  lnDiffInclSolver.setStep(timeStep);

  // Starting point for computations
  IVector x1(3);
  x1[0] = 0.0;    x1[1] = -10.3;      x1[2] = 0.03;

  // We prepare sets that know how to propagate themselves with differential inclusions
  InclRect2Set lnSet(x1),
               cwSet(x1);

  // We do the numberOfSets steps
  int numberOfSteps = 10;
  for(int i = 0; i < numberOfSteps; ++i) {
    lnSet.move(lnDiffInclSolver);
    cwSet.move(cwDiffInclSolver);
  }

  // We compute interval vector that covers given set.
  IVector lnResult = IVector(lnSet),
          cwResult = IVector(cwSet);
  std::cout.precision(16);
  std::cout << "\n\n Method based on logarithmic norms : \n " << lnResult
      << "\n  diam = " << maxDiam(lnResult) << "\n";
  std::cout << "\n\n Method based on component wise estimates : \n " << cwResult
      << "\n  diam = " << maxDiam(cwResult) << "\n";

}

int main() {
  RosslerExample();
  return 0;
} // END
