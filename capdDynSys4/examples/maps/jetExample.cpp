/**
 * \file jetExample.cpp
 * This files demonstrates how to compute third order derivatives for maps.
 */
#include <iostream>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;
using capd::autodiff::Node;

// ####################################################

/*
 * This is a map we will evaluate and differentiate (dimIn = 3, dimOut = 2, noParams = 5)
 * @param in is an array of independent variables
 * @param out is an array of dependent variables, i.e. out = f(in)
 * @param params - parameters of the map.
 */
void _f(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int noParams)
{
  out[0] = in[0]+in[1];
  out[1] = -(in[1]+in[2]);

  for(int i=0;i<noParams;++i){
    Node temp = params[i]*sqrt(sqr(out[0] - in[1]) + sqr(out[1] - in[2]));
    out[1] = params[i]*sqrt(sqr(out[0] + in[2]) + sqr(out[1] - in[0]));
    out[0] = temp;
  }
}

// ####################################################

int main(){
  int dimIn=3, dimOut=2, noParam=5;

  // this is the maximal order of derivative we request
  // default value is set to 1 if the argument is skipped
  int maxDerivativeOrder = 3;
  DMap f(_f,dimIn,dimOut,noParam,maxDerivativeOrder);

  // set parameter values
  for(int i=0;i<noParam;++i)
    f.setParameter(i, double(i+1)/9.);

  double v[] = {2,3,4};
  DVector x(dimIn,v);

  // declare an object for storing jets
  DJet jet(dimOut,dimIn,maxDerivativeOrder);
  f(x,jet);

  // print the result
  // NOTE: indices of variables must form a nondecreasing sequence!
  cout << "df1_dx0_dx0_dx2 = " << jet(1,0,0,2) << endl;


  cout << jet.toString();
}

// ####################################################

/* output:
df1_dx0_dx0_dx2 = -0.0066946

value :
   {0,0,0}  : {1.31271,2.96826}
Taylor coefficients of order 1 :
   {1,0,0}  : {-0.00665595,-0.041866}
   {0,1,0}  : {0.279001,0.116691}
   {0,0,1}  : {0.122254,0.67548}
Taylor coefficients of order 2 :
   {2,0,0}  : {-0.00986967,0.0438603}
   {1,1,0}  : {0.00505905,-0.00594214}
   {1,0,1}  : {0.00607538,-0.0394037}
   {0,2,0}  : {0.0322789,0.0145537}
   {0,1,1}  : {-0.0509478,-0.0188595}
   {0,0,2}  : {0.0175866,0.0169232}
Taylor coefficients of order 3 :
   {3,0,0}  : {0.00112803,-0.00261786}
   {2,1,0}  : {0.00424835,-0.000458252}
   {2,0,1}  : {-0.00241089,-0.0066946}
   {1,2,0}  : {-0.000898502,-0.000433624}
   {1,1,1}  : {-0.00416536,0.00259422}
   {1,0,2}  : {0.00200803,0.00729993}
   {0,3,0}  : {-0.0088823,-0.00260142}
   {0,2,1}  : {0.0123647,0.00243158}
   {0,1,2}  : {-0.0018637,-0.000114807}
   {0,0,3}  : {-0.0013343,-0.00259822}
*/
