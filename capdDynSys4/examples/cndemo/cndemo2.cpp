/////////////////////////////////////////////////////////////////////////////
//
/// @file cndemo2.cpp
///
/// @author kapela  @date 2010-02-09
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#include <iostream>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;

int main(){
  const int degree = 3;
  std::cout << "\n Michelson equation \n----------------------\n";
  // vector field
  IMap f("par:c;var:x,y,z;fun:y,z,c^2-y-0.5*x*x;",degree);
  f.setParameter("c",interval(1.0));

  // ICnSolver - numerical method for ODE integration
  int order = 20;                 // order of numerical method
  ICnOdeSolver solver(f, order);

  // initial condition
  IVector v(3);
  double sizeOfSet = 0.005;
  v[0]=0.0;    v[1]= interval(1.563-sizeOfSet, 1.563+sizeOfSet);  v[2]=0.0;
  CnRect2Set set(v, degree);

  ICnTimeMap timeMap(solver);
  interval time = 2;
  timeMap(time, set);

  cout << set.currentSet().toString();

  return 0;
}

