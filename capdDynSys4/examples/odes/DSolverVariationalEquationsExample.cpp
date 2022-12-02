#include <cstdio>
#include "capd/capdlib.h"

using namespace capd;

int main(){
  // Define an instance of class DMap that describes the vector field
  DMap pendulum("time:t;par:omega;var:x,y;fun:y,sin(omega*t)-sin(x);");
  pendulum.setParameter("omega",1.);

  // define an instance of ODE solver
  int order = 20;
  DOdeSolver solver(pendulum,order);

  // We set a fixed time step - this disables automatic step control.
  // If one really wants to use fixed time step, it is recommended to put a floating point number with many tailing zeros in the mantissa.
  solver.setStep(0.125);

  // specify initial condition
  DVector u(2);
  u[0] = 1.;
  u[1] = 2.;
  double t = 4;
  DMatrix D = DMatrix::Identity(2);

  // use one step Taylor method to integrate
  do{
    u = solver(t,u,D,D);
    printf("t=%6f, x=%6f, y=%6f, du1/dx=%6f, du1/dy=%f6, du2/dx=%6f, du2/dy=%6f\n",t,u[0],u[1],D[0][0],D[0][1],D[1][0],D[1][1]);
  }while(t<5);

  return 0;
}

/* output:
t=4.125000, x=1.236999, y=1.788258, du1/dx=0.996334, du1/dy=0.1248596, du2/dx=-0.054193, du2/dy=0.996888
t=4.250000, x=1.446319, y=1.558578, du1/dx=0.987553, du1/dy=0.2491456, du2/dx=-0.082043, du2/dy=0.991905
t=4.375000, x=1.626223, y=1.318748, du1/dx=0.976821, du1/dy=0.3730036, du2/dx=-0.085989, du2/dy=0.990894
t=4.500000, x=1.775834, y=1.074735, du1/dx=0.966890, du1/dy=0.4972216, du2/dx=-0.069940, du2/dy=0.998278
t=4.625000, x=1.894909, y=0.830721, du1/dx=0.959990, du1/dy=0.6230426, du2/dx=-0.038190, du2/dy=1.016892
t=4.750000, x=1.983628, y=0.589383, du1/dx=0.957828, du1/dy=0.7519716, du2/dx=0.005228, du2/dy=1.048134
t=4.875000, x=2.042431, y=0.352264, du1/dx=0.961636, du1/dy=0.8856116, du2/dx=0.056802, du2/dy=1.092206
t=5.000000, x=2.071902, y=0.120138, du1/dx=0.972242, du1/dy=1.0255256, du2/dx=0.113523, du2/dy=1.148294
}
*/
