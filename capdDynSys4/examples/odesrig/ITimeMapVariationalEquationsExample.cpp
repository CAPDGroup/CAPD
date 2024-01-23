#include <iostream>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;

int main(){
  cout.precision(17);

  // define an instance of class IMap that describes the vector field
  IMap pendulum("time:t;par:omega;var:x,dx;fun:dx,sin(omega*t)-sin(x);");
  pendulum.setParameter("omega",1.);

  int order = 20;
  IOdeSolver solver(pendulum,order);
  ITimeMap timeMap(solver);

  // specify initial condition
  IVector u(2);
  u[0] = 1.;
  u[1] = 2.;
  interval initTime = 4;
  C1Rect2Set set(u,initTime);

  // integrate
  interval finalTime = 8;
  u = timeMap(finalTime,set);

  // print result
  cout << "u=" << u << endl;
  cout << "monodromyMatrix=" << (IMatrix)set << endl;
  return 0;
}

/*output:
u={[-1.3297388770241241, -1.3297388770240786],[0.056757688397800897, 0.05675768839784965]}
monodromyMatrix={
{[0.25336614973537314, 0.25336614973554084],[1.3786621618302277, 1.3786621618305588]},
{[-1.0950745933673036, -1.0950745933671608],[-2.0118627006388627, -2.011862700638559]}
}
*/
