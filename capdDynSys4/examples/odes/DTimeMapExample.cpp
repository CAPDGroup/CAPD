#include <cstdio>
#include "capd/capdlib.h"

using namespace capd;

int main(){
  // define an instance of class DMap that describes the vector field
  DMap pendulum("time:t;par:omega;var:x,dx;fun:dx,sin(omega*t)-sin(x);");
  pendulum.setParameter("omega",1.);

  int order = 30;
  DOdeSolver solver(pendulum,order);
  DTimeMap timeMap(solver);

  // specify initial condition
  DVector u(2);
  u[0] = 1.;
  u[1] = 2.;
  double initTime = 4;
  double finalTime = 8;

  // use Taylor method with step control to integrate
  timeMap.stopAfterStep(true);

  do{
    u = timeMap(finalTime,u,initTime);
    printf("t=%6f, x=%6f, dx=%6f\n",initTime,u[0],u[1]);
  }while(!timeMap.completed());

  return 0;
}

/* output:
t=4.343750, x=1.584065, dx=1.379298
t=4.718750, x=1.964273, dx=0.649367
t=5.171875, x=2.065724, dx=-0.190317
t=5.703125, x=1.723517, dx=-1.077336
t=6.125000, x=1.143133, dx=-1.641492
t=6.500000, x=0.469677, dx=-1.898192
t=6.875000, x=-0.234645, dx=-1.793385
t=7.234375, x=-0.809953, dx=-1.363382
t=7.609375, x=-1.202215, dx=-0.706059
t=8.000000, x=-1.329739, dx=0.056758
*/
