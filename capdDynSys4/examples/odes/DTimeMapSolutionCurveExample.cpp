#include <iostream>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;

int main(){
  cout.precision(21);

  // define an instance of class DMap that describes the vector field
  LDMap pendulum("time:t;par:omega;var:x,dx;fun:dx,sin(omega*t)-sin(x);");
  pendulum.setParameter("omega",1.);

  int order=20;
  LDOdeSolver solver(pendulum,order);
  LDTimeMap timeMap(solver);

  // specify initial condition
  LDVector u(2);
  u[0] = 1.;
  u[1] = 2.;
  long double initTime = 4;
  long double finalTime = 10;

  // define functional object
  LDTimeMap::SolutionCurve solution(initTime);

  // and integrate
  timeMap(finalTime,u,solution);

  // then you can ask for any intermediate point on the trajectory
  cout << "domain = [" << solution.getLeftDomain() << "," << solution.getRightDomain() << "]\n";
  cout << "solution(4) =" << solution(4) << endl;
  cout << "solution(7) =" << solution(7) << endl;
  cout << "solution(10)=" << solution(10) << endl;

  return 0;
}

/* output:
domain = [4,10]
solution(4) ={1,2}
solution(7) ={-0.451828248033488868443,-1.67522823683216929268}
solution(10)={1.20483703334825757551,1.17041482661774720001}
*/
