#include <cstdio>
#include "capd/capdlib.h"

using namespace capd;

int main(){
  // define an instance of class DMap that describes the vector field
  DMap pendulum("time:t;par:omega;var:x,dx;fun:dx,sin(omega*t)-sin(x);");
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

  // use one step Taylor method to integrate
  do{
    u = solver(t,u);
    printf("t=%6f, x=%6f, dx=%6f\n",t,u[0],u[1]);
  }while(t<5);

  return 0;
}

/* output:
t=4.125000, x=1.236999, dx=1.788258
t=4.250000, x=1.446319, dx=1.558578
t=4.375000, x=1.626223, dx=1.318748
t=4.500000, x=1.775834, dx=1.074735
t=4.625000, x=1.894909, dx=0.830721
t=4.750000, x=1.983628, dx=0.589383
t=4.875000, x=2.042431, dx=0.352264
t=5.000000, x=2.071902, dx=0.120138
*/
