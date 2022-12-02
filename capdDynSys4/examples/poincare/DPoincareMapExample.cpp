#include <iostream>
#include "capd/capdlib.h"

using namespace std;
using namespace capd;

int main(){
  cout.precision(17);

  // Define vector field for the Lorenz system
  DMap lorenz("par:s,r,q;var:x,y,z;fun:s*(y-x),x*(r-z)-y,x*y-q*z;");
  lorenz.setParameter("s",10.);
  lorenz.setParameter("r",28.);
  lorenz.setParameter("q",8./3.);

  // Define solver
  DOdeSolver solver(lorenz,20);

  // Define Poincare section
  DNonlinearSection section("par:r;var:x,y,z;fun:z-r+1;");
  section.setParameter("r",28);

  // Define Poincare Map
  DPoincareMap pm(solver,section,poincare::PlusMinus);

  double returnTime=0;

  // Initial condition close to a periodic orbit
  double data[] = {-2.1473674741633753,2.0780481172873415,27};
  DVector x(3,data);

  cout << "     x=" << x << endl;
  x = pm(x,returnTime);
  cout << "  P(x)=" << x << ", halfReturnTime=" << returnTime<< endl;
  x = pm(x,returnTime);
  cout << "P^2(x)=" << x << ", fullReturnTime=" << returnTime<< endl;
  return 0;
}

/* output:
     x={-2.1473674741633753,2.0780481172873415,27}
  P(x)={2.1473675875613338,-2.0780482713004851,27}, halfReturnTime=0.77932610853932338
P^2(x)={-2.1473677282179087,2.078048081555901,27}, fullReturnTime=1.5586522078971192
*/

