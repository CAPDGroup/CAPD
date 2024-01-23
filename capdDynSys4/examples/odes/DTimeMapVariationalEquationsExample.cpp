#include <cstdio>
#include "capd/capdlib.h"

using namespace capd;

int main(){
  // define an instance of class DMap that describes the vector field
  DMap pendulum("time:t;par:omega;var:x,dx;fun:dx,sin(omega*t)-sin(x);");
  pendulum.setParameter("omega",1.);

  int order = 20;
  DOdeSolver solver(pendulum,order);
  DTimeMap timeMap(solver);

  // specify initial condition
  DVector u(2);
  u[0] = 1.;
  u[1] = 2.;
  double initTime = 4;
  double finalTime = 8;
  double data[] = {1,1,3,4}; // {{1,1},{3,4}}
  DMatrix initMatrix(2,2,data);
  DMatrix D(2,2);

  // integrate
  u = timeMap(finalTime,u,initMatrix,D,initTime);

  // print result
  printf("t=%6f, x=%6f, y=%6f, du1/dx=%6f, du1/dy=%f6, du2/dx=%6f, du2/dy=%6f\n",finalTime,u[0],u[1],D[0][0],D[0][1],D[1][0],D[1][1]);

  return 0;
}

// output: t=8.000000, x=-1.329739, y=0.056758, du1/dx=4.389353, du1/dy=5.7680156, du2/dx=-7.130663, du2/dy=-9.142525
