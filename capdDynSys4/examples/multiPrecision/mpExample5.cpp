///////////////////////////////////////////////////////////////////////
//
/// @file mpExample5.cpp
/// Example how to use multiple precision library
///
/// @author kapela  @date 2021-09-29
///
// /////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include <iostream>
#include "capd/capdlib.h"
#include "capd/mpcapdlib.h"
#include "capd/dynsys/MpStepControl.h"
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/poincare/TimeMap.hpp"

using namespace std;
using namespace capd;

typedef  capd::dynsys::OdeSolver<capd::MpIMap, capd::dynsys::MpLastTermsStepControl> MyMpIOdeSolver;
//typedef  capd::dynsys::OdeSolver<capd::MpIMap, capd::dynsys::ILastTermsStepControl> MyMpIOdeSolver;

int main(){
  MpFloat::setDefaultPrecision(1000);
  /// We set MpFloat to round numbers up in each operation

  cout.precision(60);
  try{
    // This is vector field for the Rossler system
    MpIMap vectorField("par:a,b;var:x,y,z;fun:-(y+z),x+b*y,b+z*(x-a);");

    // set chaotic parameter values
    vectorField.setParameter("a",MpInterval(57.)/10.); // a=5.7
    vectorField.setParameter("b",MpInterval(2.)/10.); // b=0.2


    // the solver, is uses high order enclosure method to verify the existence
    // of the solution. The order is set to 20.
    MyMpIOdeSolver solver(vectorField,800);

    MpFloat tol("1.e-400");
    //double tol = 1.e-20;
    solver.setAbsoluteTolerance(tol);
    solver.setRelativeTolerance(tol);

    capd::poincare::TimeMap<MyMpIOdeSolver> timeMap(solver);
    // this is a good approximation of a periodic point
    MpIVector c{ MpInterval{0.}, MpInterval{-8.3809417428298}, MpInterval{0.029590060630665}};
    // take some box around c
    c += MpInterval(-1,1)*1e-100;

    // define a doubleton representation of the interval vector c
    MpC0Rect2Set s(c);

    // we integrate the set s over the time T
    MpInterval T(1.);

    cout << "\ninitial set: " << c;
    cout << "\ndiam(initial set): " << diam(c) << endl;

    MpIVector result = timeMap(T,s);

    cout << "\n\nafter time=" << T << " the image is: " << result;
    cout << "\ndiam(image): " << diam(result) << endl << endl;
  }catch(exception& e)
  {
    cout << "\n\nException caught!\n" << e.what() << endl << endl;
  }

  return 0;
}
