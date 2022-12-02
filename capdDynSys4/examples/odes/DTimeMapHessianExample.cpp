#include <cstdio>
#include "capd/capdlib.h"

using namespace capd;

int main(){
  // define an instance of class DMap that describes the vector field
  DMap pendulum("time:t;par:omega;var:x,dx;fun:dx,sin(omega*t)-sin(x);");
  pendulum.setParameter("omega",1.);

  int order = 20;
  DC2OdeSolver solver(pendulum,order);
  DC2TimeMap timeMap(solver);

  // specify initial condition
  DVector u(2);
  u[0] = 1.;
  u[1] = 2.;
  double initTime = 4;
  double finalTime = 8;
  double initDerData[] = {1,1,3,4}; // {{1,1},{3,4}}
  double initHessData[] = {1,3,3,0,2,1}; // {dxdx,dxdy,dydy} = {{1,3,3},{0,2,1}}
  DMatrix initMatrix(2,2,initDerData);
  DHessian initHessian(2,2,initHessData);
  DMatrix D(2,2);
  DHessian H(2,2);

  // integrate
  u = timeMap(finalTime,u,initMatrix,initHessian,D,H,initTime);

  // print result
  printf("t=%6f, x=%6f, y=%6f\n",finalTime,u[0],u[1]);

  // print derivative
  printf("du1/dx=%6f, du1/dy=%f6, du2/dx=%6f, du2/dy=%6f\n",D[0][0],D[0][1],D[1][0],D[1][1]);

  // print normalized hessian
  printf(
      "d2u1/dxdx=%6f, d2u1/dxdy=%f6, d2u1/dydy=%f6, d2u1/dxdx=%6f, d2u1/dxdy=%f6, d2u1/dydy=%f6\n",
      H(0,0,0),H(0,0,1),H(0,1,1),H(1,0,0),H(1,0,1),H(1,1,1)
      );
  return 0;
}

/* output:

t=8.000000, x=-1.329739, y=0.056758
du1/dx=4.389353, du1/dy=5.7680156, du2/dx=-7.130663, du2/dy=-9.142525
d2u1/dxdx=33.076326, d2u1/dxdy=86.6810126, d2u1/dydy=54.8328166, d2u1/dxdx=-27.228978, d2u1/dxdy=-74.3629156, d2u1/dydy=-48.3478676

*/
