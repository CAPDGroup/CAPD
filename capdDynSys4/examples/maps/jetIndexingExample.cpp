/**
 * \file jetIndexingExample.cpp
 * This files demonstrates how to access higher order derivatives stored in jets.
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
  int maxDerivativeOrder = 4;
  DMap f(_f,dimIn,dimOut,noParam,maxDerivativeOrder);

  // set parameter values
  for(int i=0;i<noParam;++i)
    f.setParameter(i, double(i+1)/9.);

  double v[] = {2,3,4};
  DVector x(dimIn,v);

  // declare an object for storing jets
  DJet jet(dimOut,dimIn,maxDerivativeOrder);
  // compute partial derivatives
  f(x,jet);

  // print the result
  int data[] = {0,2,2};
  Multiindex m(3,data);
  cout << "multiindex: " << m << endl;

  cout << "vector of Taylor coefficients D^{m}f = " << jet(m) << endl;
  cout << "selected Taylor coefficient D^{m}f_1 = " << jet(1,m) << endl << endl;

  Multipointer mp(m);
  cout << "conversion to multipointer: " << mp << endl;
  cout << "vector of Taylor coefficients D^(mp}f = " << jet(m) << endl;
  cout << "selected Taylor coefficient D^(mp}f_1 = " << jet(1,m) << endl << endl;

  // print all Taylor coefficients of degree 4
  Multipointer p = jet.first(4);
  do{
    cout << p << ": " << jet(p) << endl;
  }while(jet.hasNext(p));
  return 0;
}

// ####################################################

/* output:
multiindex: {0,2,2}
vector of Taylor coefficients D^{m}f = {-0.00197489,-0.000377739}
selected Taylor coefficient D^{m}f_1 = -0.000377739

conversion to multipointer: {1,1,2,2}
vector of Taylor coefficients D^(mp}f = {-0.00197489,-0.000377739}
selected Taylor coefficient D^(mp}f_1 = -0.000377739

{0,0,0,0}: {6.64042e-05,6.03996e-05}
{0,0,0,1}: {-0.000418493,0.000195404}
{0,0,0,2}: {-0.000382953,0.00104158}
{0,0,1,1}: {-0.00100131,-0.00011813}
{0,0,1,2}: {5.53662e-06,0.000113215}
{0,0,2,2}: {0.00088786,0.000850009}
{0,1,1,1}: {-8.04876e-05,0.000219013}
{0,1,1,2}: {0.00163166,-0.000157837}
{0,1,2,2}: {-0.000185176,-0.000586786}
{0,2,2,2}: {-0.000584331,-0.00135329}
{1,1,1,1}: {0.00194541,0.00045357}
{1,1,1,2}: {-0.00135485,-0.000169508}
{1,1,2,2}: {-0.00197489,-0.000377739}
{1,2,2,2}: {0.00132892,0.000305802}
{2,2,2,2}: {-9.34493e-06,0.000436602}
*/
