#include <iostream>
#include "capd/capdlib.h"
using namespace std;
using namespace capd;

int main(){
  cout.precision(17);

  // define vector field
  IMap pendulum("time:t;par:omega;var:x,dx;fun:dx,sin(omega*t)-sin(x);");
  pendulum.setParameter("omega",1.);

  // define solver and TimeMap
  IOdeSolver solver(pendulum,20);
  ITimeMap timeMap(solver);

  IVector u(2);
  u[0] = 1.;
  u[1] = 2.;
  interval initTime = 4;

  C0HORect2Set set(u,initTime);
  timeMap.stopAfterStep(true);
  interval finalTime = 5.;
  do{
    timeMap(finalTime,set);
    cout << "currentTime: " << set.getCurrentTime() << endl;
    cout << "x(t) = " << (IVector)set << endl;
    cout << "enclosure over last time step: " << set.getLastEnclosure() << endl << endl;
  }while(!timeMap.completed());

  return 0;
}
/* output:
currentTime: [4.140625, 4.140625]
x(t) = {[1.2647227942438677, 1.2647227942438684],[1.7603531149933569, 1.7603531149933582]}
enclosure over last time step: {[0.99999999999999989, 1.2647227942438684],[1.7603531149933569, 2.0000000000000004]}

currentTime: [4.3046875, 4.3046875]
x(t) = {[1.5287156355947062, 1.528715635594708],[1.4545060500284652, 1.4545060500284688]}
enclosure over last time step: {[1.2647227942438675, 1.528715635594708],[1.4545060500284654, 1.7603531149933584]}

currentTime: [4.46875, 4.46875]
x(t) = {[1.7412935160292706, 1.741293516029274],[1.1358729403345047, 1.1358729403345105]}
enclosure over last time step: {[1.528715635594706, 1.741293516029274],[1.1358729403345049, 1.454506050028469]}

currentTime: [4.6484375, 4.6484375]
x(t) = {[1.9138452882401127, 1.913845288240118],[0.78520666257286698, 0.78520666257287419]}
enclosure over last time step: {[1.7412935160292704, 1.913845288240118],[0.78520666257286686, 1.1358729403345107]}

currentTime: [4.8359375, 4.8359375]
x(t) = {[2.0272354552591754, 2.0272354552591838],[0.42585209942088714, 0.42585209942089625]}
enclosure over last time step: {[1.9138452882401125, 2.0272354552591838],[0.42585209942088709, 0.7852066625728743]}

currentTime: [5, 5]
x(t) = {[2.0719023157504193, 2.0719023157504308],[0.12013839737212174, 0.12013839737213276]}
enclosure over last time step: {[2.0272354552591749, 2.0719023157504308],[0.12013839737212177, 0.4258520994208963]}

*/
