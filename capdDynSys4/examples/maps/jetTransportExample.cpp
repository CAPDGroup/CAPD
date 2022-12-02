/**
 * \file jetTransportExample.cpp
 * This files demonstrates how to compute composition of jets (or jet of composition of maps at given point)
 */
#include <iostream>
#include "capd/capdlib.h"

using namespace capd;
using namespace std;
using capd::autodiff::Node;

// ####################################################

/*
 * This is a map we will evaluate and differentiate (dimIn = 2, dimOut = 2, noParams = 0)
 * @param in is an array of independent variables
 * @param out is an array of dependent variables, i.e. out = f(in)
 * @param params - parameters of the map.
 */
void _f(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node /*params*/[], int/* noParams*/)
{
  Node d = sqr(in[0])+sqr(in[1]);
  out[0] = (sqr(in[0])+in[1])/d;
  out[1] = (sqr(in[1])+in[0])/d;
}

// ####################################################

int main(){
  int dimIn=2, dimOut=2, noParam=0;

  // this is the maximal order of derivative we request
  int maxDerivativeOrder = 3;
  DMap f(_f,dimIn,dimOut,noParam,maxDerivativeOrder);

  double v[] = {2,3};
  DVector x(dimIn,v);

  // declare an object for storing jets
  DJet fx(dimOut,dimIn,maxDerivativeOrder);

  // and compute jet of f(x);
  f(x,fx);
  cout << "f(x):\n" << fx.toString() << endl;

  // now we compute jet of f(f(x))
  DJet ffx = f(fx);
  cout << "f(f(x)):\n" << ffx.toString();
  return 0;
}

// ####################################################

/* output:
f(x):

value :
   {0,0}  : {0.538462,0.846154}
Taylor coefficients of order 1 :
   {1,0}  : {0.142012,-0.183432}
   {0,1}  : {-0.171598,0.0710059}
Taylor coefficients of order 2 :
   {2,0}  : {-0.00819299,-0.00864816}
   {1,1}  : {-0.0127447,0.0628129}
   {0,2}  : {0.0377788,-0.0209376}
Taylor coefficients of order 3 :
   {3,0}  : {-0.00840307,0.0167711}
   {2,1}  : {0.0209026,-0.0207976}
   {1,2}  : {-0.0166661,-0.00843808}
   {0,3}  : {-0.00423655,0.00420153}

f(f(x)):

value :
   {0,0}  : {1.12941,1.24706}
Taylor coefficients of order 1 :
   {1,0}  : {0.146505,0.0278201}
   {0,1}  : {-0.0405536,0.0289965}
Taylor coefficients of order 2 :
   {2,0}  : {-0.0084657,-0.0225223}
   {1,1}  : {-0.0165577,0.0513845}
   {0,2}  : {0.00167963,-0.0400586}
Taylor coefficients of order 3 :
   {3,0}  : {-0.023039,-0.00636561}
   {2,1}  : {0.0387984,0.0137109}
   {1,2}  : {-0.0165008,-0.0181941}
   {0,3}  : {0.00349705,0.0127827}
*/
