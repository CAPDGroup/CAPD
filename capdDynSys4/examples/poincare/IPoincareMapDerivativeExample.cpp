#include <iostream>

#include "capd/capdlib.h"

using namespace capd;
using namespace std;

// ----------------------------------- MAIN ----------------------------------------

int main()
{
  cout.precision(17);
  try{
    IMap vectorField("par:a,b;var:x,y,z;fun:-(y+z),x+b*y,b+z*(x-a);");
    vectorField.setParameter("a",interval(57)/interval(10));
    vectorField.setParameter("b",interval(2)/interval(10));

    IOdeSolver solver(vectorField,20);

    // MinusPlus means that the section is crossed with 'x' changing sign from minus to plus
    ICoordinateSection section(3,0); // the section is x=0, i.e. 0-index coordinate of 3
    IPoincareMap pm(solver,section,poincare::MinusPlus);

    // this point is very close to fixed point for the Poincare map
    interval data[] = {0.,-8.3809417428298 , 0.029590060630665};
    IVector u(3, data);
    C1Rect2Set set(u);
    interval returnTime;
    IMatrix monodromyMatrix(3,3);

    // compute Poincare map
    IVector P = pm(set,monodromyMatrix,returnTime);

    // recompute monodromy matrix to derivative of Poincare map
    IMatrix DP = pm.computeDP(P,monodromyMatrix,returnTime);

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
   u={[0, 0],[-8.3809417428297994, -8.3809417428297994],[0.029590060630665001, 0.029590060630665001]}
P(u)={[-7.3856090901169207e-13, 7.3856090901169561e-13],[-8.3809417428304709, -8.380941742829652],[0.029590060630664217, 0.029590060630669816]}
return time = [5.8810884555538525, 5.8810884555539493]
monodromyMatrix={
{[0.50685970275883696, 0.50685970276038284],[-2.4490247655247228, -2.4490247655218584],[0.42637843295482369, 0.42637843295635924]},
{[-0.59178625525463191, -0.59178625525322748],[-1.9133051621553399, -1.9133051621525994],[1.881725117656958, 1.881725117658358]},
{[0.0016796765683386335, 0.0016796765683475363],[-0.010279868615699331, -0.010279868615671371],[0.00249192754290111, 0.002491927542910677]}
}
DP(u)={
{[-1.596056620201125e-12, 1.596056620201125e-12],[-3.1072922013208881e-12, 3.1068481121110381e-12],[-1.57773794029481e-12, 1.5777934514460412e-12]},
{[-0.49005513908904019, -0.49005513908721571],[-2.4048455658571841, -2.4048455658533383],[1.9673029484794418, 1.9673029484812423]},
{[-0.00022220565884599518, -0.00022220565882650795],[-0.0010904289145296838, -0.00109042891446786],[0.00089203400376549973, 0.00089203400378485296]}
}
*/
