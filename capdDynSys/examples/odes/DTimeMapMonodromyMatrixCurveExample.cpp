#include <iostream>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;

int main(){
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
  long double data[] = {1,1,3,4};
  LDMatrix initMatrix(2,2,data); // {{1,1},{3,4}}
  LDMatrix monodromyMatrix(2,2);

  // define functional object
  LDTimeMap::SolutionCurve solution(initTime);

  // and integrate
  timeMap(finalTime,u,initMatrix,monodromyMatrix,solution);

  // then you can ask for any intermediate point on the trajectory
  cout.precision(21);
  cout << "domain = [" << solution.getLeftDomain() << "," << solution.getRightDomain() << "]\n";
  cout << "solution(4) =" << solution(4) << endl;
  cout << "solution(7) =" << solution(7) << endl;
  cout << "solution(10)=" << solution(10) << endl;

  cout << "monodromyMatrix(4)=" << solution.derivative(4) << endl;
  cout << "monodromyMatrix(7)=" << solution.derivative(7) << endl;
  cout << "monodromyMatrix(10)=" << solution.derivative(10) << endl;
  cout << monodromyMatrix << endl;
  return 0;
}

/* output:
domain = [4,10]
solution(4) ={1,2}
solution(7) ={-0.451828248033488868443,-1.67522823683216929268}
solution(10)={1.20483703334825757551,1.17041482661774720001}
monodromyMatrix(4)={
{1,1},
{3,4}
}
monodromyMatrix(7)={
{10.1969455920725050385,13.1893571903452595699},
{-3.06072164458348302597,-3.8608572219154975219}
}
monodromyMatrix(10)={
{-7.80565667723112838715,-10.0010040066793847611},
{-2.61890487701602830734,-3.48358623656441465785}
}
{
{-7.80565667723112838715,-10.0010040066793847611},
{-2.61890487701602830734,-3.48358623656441465785}
}
*/
