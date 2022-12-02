/**
 *  @file diffInclExample.cpp
 *
 *  Created on: 2009-11-13
 *      Author: kapela
 */

// it provides IVector, IMatrix, IMap
#include "capd/capdlib.h"

#include "capd/diffIncl/DiffInclusionLN.hpp"
#include "capd/diffIncl/DiffInclusionCW.hpp"
#include "capd/diffIncl/InclRect2Set.hpp"
#include "capd/diffIncl/MultiMap.h"
#include "capd/vectalg/Norm.hpp"
#include "capd/poincare/PoincareMap.hpp"
#include <iostream>

using namespace capd;

typedef capd::diffIncl::InclRect2Set<IMatrix> InclRect2Set;
typedef capd::diffIncl::MultiMap<IMap> IMultiMap;
typedef capd::diffIncl::DiffInclusionCW<IMultiMap> DiffInclusionCW;
typedef capd::diffIncl::DiffInclusionLN<IMultiMap> DiffInclusionLN;
typedef capd::poincare::PoincareMap<DiffInclusionCW> DiffPoincare;

/**
 *  In this test we integrate differential inclusion based on harmonic oscillator
 *     x' in  y + [e]
 *     y' in -x + [e]
 *   
 *  where
 *    [e] = [-eps, eps]
 * 
 *  We compute sharp enclosures for selection 
 *  x' = y
 *  y' = -x + eps sin(t) 
 *   
 */
void harmonicOscillatorTest() {

  // f is an unperturbed vector field
  IMap f("var:x,y;fun:y,-x;");
  
  // we define a perturbation
  double eps = 1.0e-10;      //  maximal forcing (perturbation)
  IMap perturb("par:e;var:x,y;fun:e,e;");
  perturb.setParameter("e", DInterval(-eps, eps));

  // We set right hand side of differential inclusion to f + perturb
  IMultiMap rhs(f, perturb);

  // We do the numberOfSets steps
  int numberOfSteps = 1000;
  double timeStep =  M_PI * 2/numberOfSteps; //  time step for integration
  int order = 10; //  order of Taylor method

  // We set up two differential inclusions (they differ in the way they handle perturbations)
  DiffInclusionCW diffInclCW(rhs, order, IMaxNorm());
  DiffInclusionLN diffInclLN(rhs, order, IEuclLNorm());
 
  diffInclCW.setStep(timeStep);
  diffInclLN.setStep(timeStep);

  IMap selection("par:e;var:x,y;time:t;fun:y,-x+e*sin(t);");
  selection.setParameter("e",DInterval(eps));
  
  ITaylor selectionSolver(selection, order);
  selectionSolver.setStep(timeStep);
  
  // Starting point for computations
  IVector x1(2);
  x1[0] = 1.0;    x1[1] = 0.0;   

  // We prepare sets that know how to propagate themselves with differential inclusions
  InclRect2Set setLN(x1),
               setCW(x1);
  
  C0Rect2Set selSet(x1);

  
  for(int i = 0; i < numberOfSteps; ++i) {
    setLN.move(diffInclLN);
    setCW.move(diffInclCW);
    selSet.move(selectionSolver);
  //}

  // We compute interval vector that covers given set.
    IVector resultLN = IVector(setLN),
            resultCW = IVector(setCW);
    IVector selResult = IVector(selSet);
//    std::cerr
//              <<  "\n selection  : " << selResult
//              <<  "\n diffInclCW : " << resultCW
//              <<  "\n diffInclLN : " << resultLN << "\n---------\n";
    if(!(subset(selResult, resultLN)) or !(subset(selResult, resultCW) )){
      
      std::cerr<< "Possible error in DiffIncl: "
              <<  "\n selection  : " << selResult
              <<  "\n diffInclCW : " << resultCW
              <<  "\n diffInclLN : " << resultLN << "\n";
      IVector intersect = intersection(intersection(selResult, resultLN), resultCW);  
    
    }
  
  
  }
}

int main() {
  std::cerr.precision(16);
  std::cout.precision(16);
  harmonicOscillatorTest();
  return 0;
} // END
