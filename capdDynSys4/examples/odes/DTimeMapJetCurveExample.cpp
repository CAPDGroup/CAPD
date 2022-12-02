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
  double finalTime = 10;

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
  // define functional object
  DCnTimeMap::SolutionCurve solution(initTime);

  // integrate
  timeMap(finalTime,initJet,J,solution);

  // print result
  std::cout << "solution(7)=" << solution.jet(7).toString() << std::endl;
  std::cout << "solution(10)="<< J.toString() << std::endl;
  return 0;
}

/* output:

solution(7)=
value :
   {0,0}  : {-0.451828,-1.67523}
Taylor coefficients of order 1 :
   {1,0}  : {10.1969,-3.06072}
   {0,1}  : {13.1894,-3.86086}
Taylor coefficients of order 2 :
   {2,0}  : {42.6352,17.2114}
   {1,1}  : {114.667,42.4899}
   {0,2}  : {73.2755,26.8991}
Taylor coefficients of order 3 :
   {3,0}  : {145.438,187.641}
   {2,1}  : {605.182,743.839}
   {1,2}  : {812.539,969.994}
   {0,3}  : {353.441,416.38}
Taylor coefficients of order 4 :
   {4,0}  : {389.64,1142.97}
   {3,1}  : {2220.14,6172.14}
   {2,2}  : {4601.83,12302}
   {1,3}  : {4123.34,10734.3}
   {0,4}  : {1351.71,3462}

solution(10)=
value :
   {0,0}  : {1.20484,1.17041}
Taylor coefficients of order 1 :
   {1,0}  : {-7.80566,-2.6189}
   {0,1}  : {-10.001,-3.48359}
Taylor coefficients of order 2 :
   {2,0}  : {-29.1662,-14.8471}
   {1,1}  : {-79.9856,-38.9769}
   {0,2}  : {-52.1961,-24.2494}
Taylor coefficients of order 3 :
   {3,0}  : {-142.07,-97.0567}
   {2,1}  : {-575.954,-389.605}
   {1,2}  : {-757.59,-513.922}
   {0,3}  : {-323.487,-223.165}
Taylor coefficients of order 4 :
   {4,0}  : {-554.512,-714.868}
   {3,1}  : {-3064.1,-3836.76}
   {2,2}  : {-6216.23,-7630.76}
   {1,3}  : {-5495.57,-6668.73}
   {0,4}  : {-1790.96,-2161.71}
*/
