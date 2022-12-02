/////////////////////////////////////////////////////////////////////////////
//
/// @file cndemo1.cpp
///
/// @author kapela  @date 2010-02-09
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group 
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#include <cmath>
#include <iostream>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;

const double PI = 3.14159265358979323846264338327;

int main(){
  std::cout << "\n Harmonic oscilator \n----------------------\n";

  const int degree = 3;
  // We define a vector field
  IMap f("var:x,y;fun:-y,x;",degree);

  int order = 10;  // order of numerical method
  int numberOfSteps = 128;
  double timeStep = PI/numberOfSteps;

  // We define ICnSolver - numerical method for ODE integration
  ICnOdeSolver solver(f, order);
  solver.setStep(timeStep); // fix time step and turn off step control

  // We set an initial condition
  IVector v(2);
  double radius = 0.1;
  v[0] = interval(1-radius, 1+radius);
  v[1] = interval(-radius, radius);

  CnMultiMatrixRect2Set set(v, degree);

  // We move set with given dynamical system (we make numberOfSteps with given timeStep)
  for(int i=0; i<numberOfSteps; ++i){
    set.move(solver);
  }

  // We extract information from the CnRect2Set.
  std::cout << "\n initial condition : " << v
            << "\n time of integration : " << numberOfSteps*timeStep
            << "\n value : \n  " << IVector(set)
            << "\n C^1 derivative :\n" << IMatrix(set)
            << "\n all derivatives" << set.currentSet().toString();
  return 0;
}

