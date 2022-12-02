/// @addtogroup nonAutTest
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file nonAutTest.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2009 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <string>
#include <cstdlib>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;

// -------------------------------------------------------------------------

bool checkVectorIntersection(const IVector& u1, const IVector& u2, interval time)
{
  interval temp;
  if( ! (intersection(u1[0],u2[0],temp) && intersection(u1[1],u2[1],temp) && intersection(u1[2],time,temp)) )
  {
    cout << "FATAL ERROR: trajectories have empty intersection!\n";
    cout << u1 << endl;
    cout << u2 << ", " << time << endl;
    return false;
  }
  return true;
}

bool checkMatrixIntersection(const IMatrix& m1, const IMatrix& m2)
{
  interval temp;
  if( !(
         intersection(m1(1,1),m2(1,1),temp) &&
         intersection(m1(1,2),m2(1,2),temp) &&
         intersection(m1(2,1),m2(2,1),temp) &&
         intersection(m1(2,2),m2(2,2),temp)
       )
    )
  {
    cout << "FATAL ERROR: the matrices have empty intersection!\n";
    cout << m1 << endl;
    cout << m2 << endl;
    return false;
  }
  return true;
}

template<class C2CoeffType>
bool checkHessianIntersection(const C2CoeffType& c1, const C2CoeffType& c2)
{
  interval temp;
  for(int i=0;i<2;++i)
    for(int j=0;j<2;++j)
      for(int k=j;k<2;++k)
      {
        if(!intersection(c1(i,j,k),c2(i,j,k),temp))
        {
          cout << "FATAL ERROR: empty intersection of a coefficient in hessian!\n";
          cout << "c1(i,j,k) = " << c1(i,j,k) << endl;
          cout << "c2(i,j,k) = " << c2(i,j,k) << endl;
          return false;
        }
      }
  return true;
}

template<class CnCoeffType>
bool checkThirdOrderIntersection(const CnCoeffType& c1, const CnCoeffType& c2)
{
  interval temp;
  for(int i=0;i<2;++i)
    for(int j=0;j<2;++j)
      for(int k=j;k<2;++k)
        for(int s=k;s<2;++s)
        {
          if(!intersection(c1(i,j,k,s),c2(i,j,k,s),temp))
          {
            cout << "FATAL ERROR: empty intersection of a coefficient in hessian!\n";
            cout << "c1(i,j,k,s) = " << c1(i,j,k,s) << endl;
            cout << "c2(i,j,k,s) = " << c2(i,j,k,s) << endl;
            return false;
          }
        }
  return true;
}

// -------------------------------------------------------------------------

void showC0Information(const IVector& u1, const IVector& u2, int numberOfSteps, interval step, interval time)
{
  cout << "During " << numberOfSteps << " steps with step=" << step << " trajectories have nonempty intersection.\nImages at the end are:\n";
  cout << "Solver:      " << u1 << endl;
  cout << "NonAutSolver:" << u2 << ", " << time << endl;
}

void showC1Information(const IMatrix& m1, const IMatrix& m2)
{
  cout << "Monodromy matrix Solver:      " << m1 << endl;
  cout << "Monodromy matrix NonAutSolver:" << m2 << endl;
}

template<class C2CoeffType>
void printHessian(const C2CoeffType& c)
{
  for(int i=0;i<2;++i)
  {
    for(int j=0;j<2;++j)
      for(int k=j;k<2;++k)
        cout << c(i,j,k) << " ";
    cout << endl;
  }
}

template<class C2CoeffType>
void showC2Information(const C2CoeffType& c1, const C2CoeffType& c2)
{
  cout << "\nHessian Solver:\n";
  printHessian(c1);
  cout << "\nHessian NonAutSolver:\n";
  printHessian(c2);
}

template<class CnCoeffType>
void printC3Coeff(const CnCoeffType& c)
{
  for(int i=0;i<2;++i)
  {
    for(int j=0;j<2;++j)
      for(int k=j;k<2;++k)
        for(int s=k;s<2;++s)
          cout << c(i,j,k,s) << " ";
    cout << endl;
  }
}

template<class CnCoeffType>
void showC3Information(const CnCoeffType& c1, const CnCoeffType& c2)
{
  cout << "\nThird order derivatives Solver:\n";
  printC3Coeff(c1);
  cout << "\nThird order derivatives NonAutSolver:\n";
  printC3Coeff(c2);
}

// -------------------------- C^0 tests -------------------------------------

bool checkIntersection(C0Rect2Set&s, C0Rect2Set& nas)
{
  return checkVectorIntersection(IVector(s),IVector(nas),nas.getCurrentTime());
}

void showFinalInformation(C0Rect2Set&s, C0Rect2Set& nas, int numberOfSteps, interval step)
{
  showC0Information(IVector(s),IVector(nas),numberOfSteps,step,nas.getCurrentTime());
  cout << endl;
}

// -------------------------- C^1 tests -------------------------------------

bool checkIntersection(C1Rect2Set&s, C1Rect2Set& nas)
{
   return checkVectorIntersection(IVector(s),IVector(nas),nas.getCurrentTime())
         && checkMatrixIntersection(IMatrix(s),IMatrix(nas));
}

void showFinalInformation(C1Rect2Set&s, C1Rect2Set& nas, int numberOfSteps, interval step)
{
  showC0Information(IVector(s),IVector(nas),numberOfSteps,step,nas.getCurrentTime());
  showC1Information(IMatrix(s),IMatrix(nas));
  cout << endl;
}

// -------------------------- C^2 tests -------------------------------------

bool checkIntersection(C2Rect2Set&s, C2Rect2Set& nas)
{
   return checkVectorIntersection(IVector(s),IVector(nas),nas.getCurrentTime())
         && checkMatrixIntersection(IMatrix(s),IMatrix(nas))
         && checkHessianIntersection(IHessian(s),IHessian(nas));
}

void showFinalInformation(C2Rect2Set&s, C2Rect2Set& nas, int numberOfSteps, interval step)
{
  showC0Information(IVector(s),IVector(nas),numberOfSteps,step,nas.getCurrentTime());
  showC1Information(IMatrix(s),IMatrix(nas));
  showC2Information(IHessian(s),IHessian(nas));
  cout << endl;
}

// -------------------------- C^3 tests -------------------------------------

bool checkIntersection(CnRect2Set&s, CnRect2Set& nas)
{
   return checkVectorIntersection(IVector(s),IVector(nas),nas.getCurrentTime())
         && checkMatrixIntersection(IMatrix(s),IMatrix(nas))
         && checkHessianIntersection(s.currentSet(),nas.currentSet())
         && checkThirdOrderIntersection(s.currentSet(),nas.currentSet());
}

void showFinalInformation(CnRect2Set&s, CnRect2Set& nas, int numberOfSteps, interval step)
{
  showC0Information(IVector(s),IVector(nas),numberOfSteps,step,nas.getCurrentTime());
  showC1Information(IMatrix(s),IMatrix(nas));
  showC2Information(s.currentSet(),nas.currentSet());
  showC3Information(s.currentSet(),nas.currentSet());
  cout << endl;
}

// -------------------------------------------------------------------------

template<class Solver,class Set1, class Set2>
void Test(Solver& T, Solver& naT, Set1& s, Set2& nas, int order, int numberOfSteps, string description)
{
  naT.setOrder(order);
  T.setOrder(order);

  cout << "\n------------------------------------------------------------------------------------------------------------\n\n";
  cout << "Computing " << description << " trajectory using Solver and NonAutSolver for the forced pendulum.\nInitial point: " << IVector(s) << endl;
  for(int i=0;i<numberOfSteps;++i)
  {
    s.move(T);
    nas.move(naT);
    if( ! checkIntersection(s,nas) )
      exit(1);
  }
  showFinalInformation(s,nas,numberOfSteps,T.getStep());
}

// -------------------------------------------------------------------------

int main(int, char *[])
{
  try
  {
    cout.precision(12);
    int order=6;
    interval step(0.125);

    string nonAutPendulumFormula = "time:t;par:beta;var:x,dx;fun:dx,cos(t)-beta*x-sin(x);";
    string pendulumFormula = "par:beta;var:x,dx,t;fun:dx,cos(t)-beta*x-sin(x),1;";
    IMap naVF(nonAutPendulumFormula,3);
    IMap VF(pendulumFormula,3);
    naVF.setParameter("beta",0.1);
    VF.setParameter("beta",0.1);

    ICnOdeSolver T(VF,order);
    ICnOdeSolver naT(naVF,order);
    T.setStep(step);
    naT.setStep(step);

    interval coeff[] = {interval(1),interval(2),interval(0)};
    IVector v(3,coeff);
    IVector nav(2,coeff);

// C0 Test
    C0Rect2Set s(v);
    C0Rect2Set nas(nav,0.);
    Test(T,naT,s,nas,10,100,"C0");

// C1Test
    C1Rect2Set c1s(v);
    C1Rect2Set c1nas(nav,0.);
    Test(T,naT,c1s,c1nas,10,100,"C1");

// C2Test
    C2Rect2Set c2s(v);
    C2Rect2Set c2nas(nav,0.);
    Test(T,naT,c2s,c2nas,10,100,"C2");

// C3Test
    CnRect2Set c3s(v,3);
    CnRect2Set c3nas(nav,3,0.);
    Test(T,naT,c3s,c3nas,10,100,"C3");

  }catch(exception& e)
  {
    cout << "Exception caught: " << e.what() << endl;
    return -1;
  }
  return 0;
}


/// @}
