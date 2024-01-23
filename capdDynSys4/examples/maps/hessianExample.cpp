/**
 * \file hessianExample.cpp
 * This files demonstrates how to compute second order derivatives for maps.
 */
#include <iostream>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;
using capd::autodiff::Node;

// ####################################################

/*
 * This is a map we will evaluate and differentiate (dimIn = 4, dimOut = 2, noParams = 5)
 * @param in is an array of independent variables
 * @param out is an array of dependent variables, i.e. out = f(in)
 * @param params - parameters of the map.
 */
void _f(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int noParams)
{
  out[0] = in[0]+in[1];
  out[1] = -(in[2]+in[3]);

  for(int i=0;i<noParams;++i){
    Node temp = params[i]*sqrt(sqr(out[0] - in[0]) + sqr(out[1] - in[3]));
    out[1] = params[i]*sqrt(sqr(out[0] + in[1]) + sqr(out[1] - in[2]));
    out[0] = temp;
  }
}

// ####################################################

int main(){
  int dimIn=4, dimOut=2, noParam=5;

  // this is the maximal order of derivative we request
  // default value is set to 1 if the argument is skipped
  int maxDerivativeOrder = 2;
  IMap f(_f,dimIn,dimOut,noParam,maxDerivativeOrder);

  // set parameter value that encloses 1/9, 2/9, etc
  for(int i=0;i<noParam;++i)
    f.setParameter(i, interval(i+1)/interval(9));

  interval v[] = {2,3,4,5};
  IVector x(dimIn,v);

  // declare an object for storing derivative and hessian
  IMatrix Df(dimOut,dimIn);
  IHessian Hf(dimOut,dimIn);

  // simultaneous computation of value, derivative and normalized hessian
  // NOTE! Hf contains second order Taylor coefficients of f at x, i.e. normalized derivatives.
  IVector y = f(x,Df,Hf);

  // print value and derivative of f at x
  cout.precision(17);
  cout << "y=" << y << endl;
  cout << "Df=" << Df << endl;

  // print normalized second order derivatives
  for(int fi=0;fi<dimOut;++fi)
    for(int dx1=0;dx1<dimIn;++dx1)
      for(int dx2=dx1;dx2<dimIn;++dx2)
        cout << "Hf(" << fi << "," << dx1 << "," << dx2 << ")=" << Hf(fi,dx1,dx2) << endl;
}

// ####################################################

/* output:
y={[1.5666313012267541, 1.5666313012267572],[2.7160497654235312, 2.7160497654235356]}
Df={
{[0.060484484806883045, 0.060484484806884051],[-0.16642557749344458, -0.16642557749344236],[-0.086398550591301623, -0.086398550591300416],[0.45810665329170269, 0.45810665329170663]},
{[0.029096263766104483, 0.02909626376610501],[0.40574188883571083, 0.40574188883571372],[0.13126754262517129, 0.13126754262517268],[0.18311228017669873, 0.18311228017670086]}
}
Hf(0,0,0)=[0.074118944476562462, 0.074118944476563184]
Hf(0,0,1)=[0.021315780883087665, 0.021315780883088681]
Hf(0,0,2)=[0.01769262837118854, 0.017692628371189165]
Hf(0,0,3)=[-0.08623872680805543, -0.086238726808053057]
Hf(0,1,1)=[-0.0053748330978258384, -0.0053748330978252252]
Hf(0,1,2)=[0.019554811781480465, 0.019554811781482002]
Hf(0,1,3)=[-0.017720362061031517, -0.017720362061027721]
Hf(0,2,2)=[-0.0091193387205848984, -0.0091193387205845741]
Hf(0,2,3)=[-0.004218996464429756, -0.0042189964644276561]
Hf(0,3,3)=[0.024251452565690394, 0.024251452565692031]
Hf(1,0,0)=[0.020245143945059545, 0.020245143945059847]
Hf(1,0,1)=[0.0087801679885119114, 0.0087801679885125602]
Hf(1,0,2)=[-0.0021131490236710707, -0.0021131490236707558]
Hf(1,0,3)=[-0.019773696730218816, -0.019773696730217935]
Hf(1,1,1)=[0.015870027950691565, 0.015870027950692377]
Hf(1,1,2)=[-0.04225008508394984, -0.042250085083948029]
Hf(1,1,3)=[0.01124396733092266, 0.011243967330925161]
Hf(1,2,2)=[0.028717542333811172, 0.028717542333811609]
Hf(1,2,3)=[-0.019752757074261142, -0.019752757074259841]
Hf(1,3,3)=[0.0084826519764703637, 0.0084826519764710315]
*/
