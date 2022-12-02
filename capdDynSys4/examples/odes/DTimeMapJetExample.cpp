#include <iostream>
#include "capd/capdlib.h"

using namespace capd;

int main(){
  // maximal degree of truncated Taylor series
  const int degree = 4;
  // define an instance of class DMap that describes the vector field
  DMap pendulum("time:t;par:omega;var:x,dx;fun:dx,sin(omega*t)-sin(x);",degree);
  pendulum.setParameter("omega",1.);

  int order = 20;
  DCnOdeSolver solver(pendulum,order);
  DCnTimeMap timeMap(solver);

  double initTime = 4;
  double finalTime = 8;

  // specify initial condition
  DJet initJet(2,degree);
  // set initial conditions
  // main ODE
  initJet(0) = 1; initJet(1) = 2;
  // init matrix
  initJet(0,0) = 1; initJet(0,1) = 1; initJet(1,0) = 3; initJet(1,1) = 4;
  // init hessian
  initJet(0,0,0) = 1; initJet(0,0,1) = 3; initJet(0,1,1) = 3; initJet(1,0,0) = 0; initJet(1,0,1) = 2; initJet(1,1,1) = 1;
  // higher order coefficients remain zero

  DJet J(2,degree);

  // integrate
  timeMap(finalTime,initJet,J,initTime);

  // print result
  std::cout << J.toString() << std::endl;
  return 0;
}

/* output:

value :
   {0,0}  : {-1.32974,0.0567577}
Taylor coefficients of order 1 :
   {1,0}  : {4.38935,-7.13066}
   {0,1}  : {5.76801,-9.14253}
Taylor coefficients of order 2 :
   {2,0}  : {33.0763,-27.229}
   {1,1}  : {86.681,-74.3629}
   {0,2}  : {54.8328,-48.3479}
Taylor coefficients of order 3 :
   {3,0}  : {181.037,-106.984}
   {2,1}  : {735.217,-438.7}
   {1,2}  : {975.287,-580.696}
   {0,3}  : {423.543,-248.531}
Taylor coefficients of order 4 :
   {4,0}  : {1035.45,-207.342}
   {3,1}  : {5608.26,-1220.75}
   {2,2}  : {11214.2,-2581.6}
   {1,3}  : {9815.9,-2339.67}
   {0,4}  : {3175.13,-771.844}
*/
