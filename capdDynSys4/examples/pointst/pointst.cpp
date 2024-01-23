/// @addtogroup poincare
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file pointst.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include <stdexcept>
#include <cmath>

#include "capd/krak/krak.h"
#include "capd/capdlib.h"
using namespace capd;

double minx = -2.2;
double maxx = 1.6;
double miny = -2.6;
double maxy = 2.;

Frame fr, txt;

// -----------------------------------------------------------------

void initGraphics() {
  openGW(900, 700, WHITE, BLACK);
  rootFrame.Clear();
  txt = Frame(0, 0, 595, 130);
  fr = Frame(5, 135, 555, 585, WHITE, BLACK);
  fr.setWorldCoord(minx, miny, maxx, maxy);
  fr.line(0.0, miny, 0.0, maxy, BLACK);
  fr.line(minx, 0.0, maxx, 0.0, BLACK);
}

// -----------------------------------------------------------------

int main(int, char *[]) {
  initGraphics();
  try {
    double grid = 30;
    int order = 20;
    IMap iVectorField = "var:x,y,z;fun:y,z,1-y-x^2/2;";
    ICoordinateSection iSection(3, 2); // z=0
    IOdeSolver iSolver(iVectorField, order);
    IPoincareMap iPM(iSolver, iSection);

    DMap dVectorField = "var:x,y,z;fun:y,z,1-y-x^2/2;";
    // z=0 just to demonstrate another possibility of defining a section
    DNonlinearSection dSection = "var:x,y,z;fun:z;";
    DOdeSolver dSolver(dVectorField, order);
    DPoincareMap dPM(dSolver, dSection);

    IVector iv(3);
    iv[0] = iv[2] = interval(0.0);
    iv[1] = interval(15, 16) / interval(10);

    fr << At(11, 55) << "y";
    fr << At(0, 33) << "y'";
    txt << "Test of class PoincareMap.\n";
    txt << "ODE: Michelson system, order:" << order << ", Poincare section: y\"=0\n\n";
    txt << "P - Poincare return map,    The set s=" << iv << "\n";
    txt << "We compute P(s) and P^2(s)\n";

    int c = RED;
    txt.SetFgColor(c);
    txt << "using set arithmetic and class PoincareMap\n";
    fr << At(7, 20) << "P^2(s)";
    fr << At(26, 26) << "P(s)";

    interval part = interval(iv[1].rightBound() - iv[1].leftBound()) / interval(grid);

    for (int i = 0; i < grid; i++) {
      IVector w(interval(0.), iv[1].leftBound() + interval(i, i + 1) * part, interval(0.));
      C0HOTripletonSet set(w);
      for (int j = 0; j < 2; ++j) {
        w = iPM(set);
        fr.boxFill(w[0].leftBound(), w[1].leftBound(), w[0].rightBound(), w[1].rightBound(), c);
        fr.box(w[0].leftBound(), w[1].leftBound(), w[0].rightBound(), w[1].rightBound());
      }
    }

    c = GREEN;
    txt.SetFgColor(c);
    txt.SetBgColor(BLACK);
    txt << "using vector arithmetic and class PoincareMap\n";
    grid *= 30;

    for (int i = 0; i <= grid; i++) {
      DVector v(0, 1.5, 0.);
      v[1] += 0.1 * i / grid;
      for (int j = 0; j < 2; ++j) {
        v = dPM(v);
        fr.dot(v[0], v[1], c);
      }
    }
    waitBt();
  } catch (std::exception& e) {
    rootFrame << "\n" << e.what();
    waitBt();
  }

  closeGW();
  return 0;
}

/// @}
