#include <iostream>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;

int main(){
  // define an instance of class IMap that describes the vector field
  IMap pendulum("time:t;par:omega;var:x,dx;fun:dx,sin(omega*t)-sin(x);");
  pendulum.setParameter("omega",1.);

  int order=20;
  IOdeSolver solver(pendulum,order);
  ITimeMap timeMap(solver);

  // specify initial condition
  IVector u(2);
  u[0] = 1.;
  u[1] = 2.;
  interval initTime = 4;
  C1Rect2Set set(u,initTime);

  // define functional object
  ITimeMap::SolutionCurve solution(initTime);

  // and integrate
  interval finalTime = 10;
  timeMap(finalTime,set,solution);

  // then you can ask for any intermediate point on the trajectory
  cout.precision(17);
  cout << "domain = [" << solution.getLeftDomain() << "," << solution.getRightDomain() << "]\n";
  cout << "solution(4) =" << solution(4) << endl;
  cout << "solution(5,5.125) =" << solution(interval(5,5.125)) << endl;
  cout << "solution(10)=" << solution(10) << endl;

  cout << "monodromyMatrix(4)=" << solution.derivative(4) << endl;
  cout << "monodromyMatrix(5,5.125)=" << solution.derivative(interval(5,5.125)) << endl;
  cout << "monodromyMatrix(10)=" << solution.derivative(10) << endl;
  return 0;
}

/* output:
domain = [4,10]
solution(4) ={[1, 1],[2, 2]}
solution(5,5.125) ={[2.0659018710648156, 2.0826100235006182],[-0.10794601246272083, 0.12123454798834073]}
solution(10)={[1.2048370333482177, 1.2048370333482994],[1.1704148266177106, 1.1704148266177841]}
monodromyMatrix(4)={
{[1, 1],[-0, 0]},
{[-0, 0],[1, 1]}
}
monodromyMatrix(5,5.125)={
{[0.97223163741647423, 0.99012219196297169],[1.0255143836882481, 1.1731157012509985]},
{[0.11326407278286743, 0.17284240065817016],[1.1480297669989168, 1.2147644133774944]}
}
monodromyMatrix(10)={
{[-1.2196146888865109, -1.2196146888862074],[-2.1953473294485577, -2.1953473294479529]},
{[-0.024860798371007034, -0.024860798370730811],[-0.86468135954864922, -0.86468135954812031]}
}
*/
