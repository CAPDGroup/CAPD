#include <iostream>

#include "capd/capdlib.h"

using namespace capd;
using namespace std;

// ----------------------------------- MAIN ----------------------------------------

int main()
{
  cout.precision(17);
  try{
    DMap vectorField("par:a,b;var:x,y,z;fun:-(y+z),x+b*y,b+z*(x-a);");
    vectorField.setParameter("a",5.7);
    vectorField.setParameter("b",0.2);

    DOdeSolver solver(vectorField,20);

    // MinusPlus means that the section is crossed with 'x' changing sign from minus to plus
    DCoordinateSection section(3,0); //  the section is x=0, 0-index from 3
    DPoincareMap pm(solver,section,poincare::MinusPlus);

    // this point is very close to fixed point for the Poincare map
    double data[] = {0.,-8.3809417428298 , 0.029590060630665};
    DVector u(3, data);
    double returnTime = 0.;
    DMatrix monodromyMatrix(3,3);

    // compute Poincare map
    DVector P = pm(u,monodromyMatrix,returnTime);

    // recompute monodromy matrix to derivative of Poincare map
    DMatrix DP = pm.computeDP(P,monodromyMatrix,returnTime);

    // print results
    cout << "   u=" << u << endl;
    cout << "P(u)=" << P << endl;
    cout << "return time = " << returnTime << endl;
    cout << "monodromyMatrix=" << monodromyMatrix << endl;
    cout << "DP(u)=" << DP << endl;
   }catch(exception& e)
  {
    cout << "\n\nException caught: "<< e.what() << endl;
  }
  return 0;
}

/* output:
   u={0,-8.3809417428297994,0.029590060630665001}
P(u)={-5.5511151231257827e-17,-8.3809417428300641,0.029590060630667017}
return time = 5.8810884555538996
monodromyMatrix={
{0.50685970275960168,-2.4490247655232884,0.4263784329575836},
{-0.59178625525395556,-1.9133051621539718,1.8817251176642307},
{0.0016796765683430427,-0.010279868615685349,0.0024919275429165304}
}
DP(u)={
{0,0,0},
{-0.49005513908815546,-2.404845565855263,1.9673029484873146},
{-0.00022220565883626208,-0.0010904289144987797,0.00089203400377833864}
}
*/
