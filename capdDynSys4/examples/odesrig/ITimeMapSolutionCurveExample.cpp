#include <iostream>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;

int main(){
  cout.precision(17);
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
  double initTime = 4;
  double finalTime = 10;
  C0Rect2Set set(u,initTime);

  // define functional object
  ITimeMap::SolutionCurve solution(initTime);

  // and integrate
  timeMap(finalTime,set,solution);

  // then you can ask for any intermediate point on the trajectory
  cout << "domain = [" << solution.getLeftDomain() << "," << solution.getRightDomain() << "]\n";
  cout << "solution(4) =" << solution(4) << endl;
  cout << "solution(5.,5.125) = " << solution(interval(5,5.125)) << endl;
  cout << "solution(7,8) =" << solution(interval(7,8)) << endl;
  cout << "solution(10)=" << solution(10) << endl;

  return 0;
}

/* output:
domain = [4,10]
solution(4) ={[1, 1],[2, 2]}
solution(5.,5.125) = {[2.061890364489789, 2.0922917471924896],[-0.10693311935888665, 0.12198962338481042]}
solution(7,8) ={[-1.3710853587244543, -0.44825760563602185],[-1.6756032292973486, 0.056820499411927523]}
solution(10)={[1.2048370333482143, 1.2048370333483036],[1.1704148266177143, 1.170414826617781]}
*/
