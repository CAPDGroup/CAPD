#include <iostream>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;

void print(const LDHessian& h, long double t){
  cout << "hessian(" << t << ")={";
  for(LDHessian::size_type i=0;i<h.imageDimension();++i)
    for(LDHessian::size_type j=0;j<h.dimension();++j)
      for(LDHessian::size_type c=j;c<h.dimension();++c)
        cout << " " << h(i,j,c);

  cout << " }\n";
}

int main(){
  // define an instance of class DMap that describes the vector field
  LDMap pendulum("time:t;par:omega;var:x,dx;fun:dx,sin(omega*t)-sin(x);");
  pendulum.setParameter("omega",1.);

  int order=20;
  LDC2OdeSolver solver(pendulum,order);
  LDC2TimeMap timeMap(solver);

  // specify initial condition
  LDVector u(2);
  u[0] = 1.;
  u[1] = 2.;
  long double initTime = 4;
  long double finalTime = 10;
  long double initDerData[] = {1,1,3,4};
  long double initHessData[] = {1,3,3,0,2,1}; // {dxdx,dxdy,dydy} = {{1,3,3},{0,2,1}}
  LDMatrix initMatrix(2,2,initDerData); // {{1,1},{3,4}}
  LDHessian initHessian(2,2,initHessData);

  LDMatrix D(2,2);
  LDHessian H(2);

  // define functional object
  LDC2TimeMap::SolutionCurve solution(initTime);

  // and integrate
  timeMap(finalTime,u,initMatrix,initHessian,D,H,solution);

  // then you can ask for any intermediate point on the trajectory
  cout.precision(21);
  cout << "domain = [" << solution.getLeftDomain() << "," << solution.getRightDomain() << "]\n";
  cout << "solution(4) =" << solution(4) << endl;
  cout << "solution(7) =" << solution(7) << endl;
  cout << "solution(10)=" << solution(10) << endl;

  cout << "monodromyMatrix(4)=" << solution.derivative(4) << endl;
  cout << "monodromyMatrix(7)=" << solution.derivative(7) << endl;
  cout << "monodromyMatrix(10)=" << solution.derivative(10) << endl;
  cout << D << endl;

  print(solution.hessian(4),4);
  print(solution.hessian(7),7);
  print(solution.hessian(10),10);
  print(H,10);
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
hessian(4)={ 1 3 3 0 2 1 }
hessian(7)={ 42.6352224242584004026 114.667151261346096119 73.2754844459227481224 17.2114064486085741219 42.4899350405156879182 26.8990700123135739949 }
hessian(10)={ -29.1661902479786202484 -79.9855736889665852526 -52.1961182105122759693 -14.8470974997925887675 -38.9769249314583959565 -24.2494060090330983866 }
hessian(10)={ -29.1661902479786202484 -79.9855736889665852526 -52.1961182105122759693 -14.8470974997925887675 -38.9769249314583959565 -24.2494060090330983866 }
*/
