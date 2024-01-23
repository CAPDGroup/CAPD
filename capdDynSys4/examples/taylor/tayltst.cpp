/// @addtogroup tayltst
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file tayltst.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#include <cmath>
#include <stdexcept>
#include <fstream>

#include "capd/krak/krak.h"
#include "capd/capdlib.h"
using namespace capd;

double minx=-2.6;
double maxx=2.6;
double miny=-2.6;
double maxy=2.6;

Frame fr, txt;
double step=0.25;
int order = 20;

void axes(Frame &fr)
{
  fr.line(0.0,miny,0.0,maxy,BLACK);
  fr.line(minx,0.0,maxx,0.0,BLACK);
}

// -------------------------------------------------------------------------

void initGraph(double step, int order)
{
  openGW(1000,800,WHITE,BLACK);
  fr = Frame(5,115,545,655,WHITE,BLACK);
  txt = Frame(0,660,860,780,WHITE,BLACK);
  fr.setWorldCoord(minx,miny,maxx,maxy);
  rootFrame.Clear();
  axes(fr);
  rootFrame << "Comparison of various epresentations and algorithms for IVP using class IMap and IOdeSolver.\n";
  rootFrame << "Evaluation of the set and its projection onto plane y\"=0.\n\n";
  rootFrame << "ODE: Michelson system\n";
  rootFrame << "Step control turned off. Constant time step: " << step << ",  order: " << order << "\n";
}

// -------------------------------------------------------------------------
template<class SetT>
void computeRigorousTrajectory(SetT &set, IOdeSolver &solver, int color, frstring description)
{
  int length=static_cast<int>(8./step);

  txt.SetFgColor(color);
  txt << "interval test: " << description;

  int i;
  for(i=0;i<length;i++)
  {
    IVector w = IVector(set);
    fr.boxFill(w[0].leftBound(),w[1].leftBound(),w[0].rightBound(),w[1].rightBound(),color);
    fr.box(w[0].leftBound(),w[1].leftBound(),w[0].rightBound(),w[1].rightBound());
    set.move(solver);
  }
  IVector w(set);
  txt.SetFgColor(BLACK);
  txt << ", diam:" << diam(w) << "\n"; 
}

// -------------------------------------------------------------------------

void computeNonrigorousTrajectory(const DVector &center, double sizeOfSet, int grid)
{
  DMap f = "var:x,y,z;fun:y,z,1-y-x^2/2;";
  DOdeSolver solver(f,order);
  // disable step control and set fixed time step
  solver.setStep(step);
  int length= static_cast<int>(8./step);

  txt.SetFgColor(BLUE);
  txt << "vector test: nonrigorous simulation in yellow";
  int i,j;
  DVector left = center + DVector(0.,-sizeOfSet,0.);
  double step = 2*sizeOfSet/(double)grid;
  for(j=0;j<=grid;j++)
  {
    DVector P = left + DVector(0.,j*step,0.);
    for(i=0;i<length;i++)
    {
      fr.dot(P[0],P[1],YELLOW);
      P = solver(P);
    }
  }
}

// -------------------------------------------------------------------------

int main(int, char *[])
{
  initGraph(step,order);

  try
  {
    DVector z(0.,1.563,0.);
    IVector v(0.,1.563,0.);
    double sizeOfSet = 0.0075;

    DInterval iv(-sizeOfSet,sizeOfSet);
    IVector r(0.,iv,0.);

    C0Rect2Set doubleton(v,r);
    C0HORect2Set ho_doubleton(v,r);
    C0TripletonSet tripleton(v,r);
    C0HOTripletonSet ho_tripleton(v,r);

    IMap f="par:c;var:x,y,z;fun:y,z,c^2-y-x^2/2;";
    f.setParameter("c", DInterval(1.0));
    IOdeSolver solver(f,order);
    solver.setStep(step);

    computeRigorousTrajectory(doubleton,solver,RED,"C0Rect2Set");
    computeRigorousTrajectory(tripleton,solver,GREEN,"C0TripletonSet");
    computeRigorousTrajectory(ho_doubleton,solver,MAGENTA,"C0HORect2Set");
    computeRigorousTrajectory(ho_tripleton,solver,VIOLET,"C0HOTripletonSet");
    computeNonrigorousTrajectory(z,sizeOfSet,40);
    waitBt();

  }catch(std::exception& e)
  {
    std::ofstream fileStream;
    fileStream.open("report");
    fileStream << e.what();
    fileStream.close();
    rootFrame << "\n\nException caught! See 'report' file for details.";
    waitBt();
  }
  closeGW();
  return 0;
}


/// @}
