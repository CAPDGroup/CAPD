//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file mittagLeffler12Test.cpp
///
/// @author Jonathan Jaquette
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

//#define BOOST_TEST_MODULE mittagLeffler12Test
#include "compare.h"
BOOST_AUTO_TEST_SUITE(mittagLeffler12Suite)

// f(x,y) = mittagLeffler12(x+y) = E_{1,2}(x+y) = (exp(x+y)-1)/(x+y)
//
// Since the argument is u = x+y (linear), all mixed partials reduce to:
//   D^n_x D^m_y f = E_{1,2}^{(n+m)}(x+y)
//   jet[n,m] = E_{1,2}^{(n+m)}(x0+y0) / (n! * m!)
//
// Closed form for derivatives of E_{1,2}(u) = (exp(u)-1)/u:
//   k=0: (exp(u)-1)/u
//   k=1: (exp(u)*(u-1)+1)/u^2
//   k=2: (exp(u)*(u^2-2u+2)-2)/u^3
//   k=3: (exp(u)*(u^3-3u^2+6u-6)+6)/u^4
//   k=4: (exp(u)*(u^4-4u^3+12u^2-24u+24)-24)/u^5
//   k=5: (exp(u)*(u^5-5u^4+20u^3-60u^2+120u-120)+120)/u^6
//
// These formulas are valid for u != 0.  Test points are chosen so that u != 0.

std::vector<double> computeML12Der(MapType::VectorType & u){
  capd::rounding::DoubleRounding::roundNearest();
  double x  = u[0].leftBound();
  double y  = u[1].leftBound();
  double u0 = x + y;
  double E  = exp(u0);
  double u2 = u0*u0, u3 = u2*u0, u4 = u3*u0, u5 = u4*u0, u6 = u5*u0;

  double d0 = (E - 1) / u0;
  double d1 = (E*(u0 - 1) + 1) / u2;
  double d2 = (E*(u2 - 2*u0 + 2) - 2) / u3;
  double d3 = (E*(u3 - 3*u2 + 6*u0 - 6) + 6) / u4;
  double d4 = (E*(u4 - 4*u3 + 12*u2 - 24*u0 + 24) - 24) / u5;
  double d5 = (E*(u5 - 5*u4 + 20*u3 - 60*u2 + 120*u0 - 120) + 120) / u6;

  // Ordering matches CAPD JetType: row m, column n gives jet[m-n, n] for n=0..m
  double r[] = {
    d0,                                      // [0,0]
    d1,       d1,                            // [1,0] [0,1]
    d2/2.,    d2,     d2/2.,                 // [2,0] [1,1] [0,2]
    d3/6.,    d3/2.,  d3/2.,  d3/6.,         // [3,0] [2,1] [1,2] [0,3]
    d4/24.,   d4/6.,  d4/4.,  d4/6., d4/24.,// [4,0] [3,1] [2,2] [1,3] [0,4]
    d5/120.,  d5/24., d5/12., d5/12., d5/24., d5/120. // [5,0]..[0,5]
  };
  return std::vector<double>(r, r + sizeof(r)/sizeof(double));
}

BOOST_AUTO_TEST_CASE(xmittagLeffler12)
{
  std::string txt = "var:x,y;fun:mittagLeffler12(x+y);",
              msg = "Function \"" + txt + "\"  (x,y) = ";
  MapType f(txt, 5);
  VectorType x(2);
  JetType df(1, 2, 5);

  x[0] = 0.5;  x[1] = 0.5;   // u = 1.0
  std::vector<double> expected = computeML12Der(x);
  f(x, df);
  compareResults(expected, df, msg + "(0.5, 0.5)");

  x[0] = 1.5;  x[1] = 0.5;   // u = 2.0
  expected = computeML12Der(x);
  f(x, df);
  compareResults(expected, df, msg + "(1.5, 0.5)");

  x[0] = 0.5;  x[1] = -1.5;  // u = -1.0
  expected = computeML12Der(x);
  f(x, df);
  compareResults(expected, df, msg + "(0.5, -1.5)");

  // Verify the alias ml12 produces identical results
  MapType g("var:x,y;fun:ml12(x+y);", 5);
  x[0] = 0.5;  x[1] = 0.5;
  expected = computeML12Der(x);
  g(x, df);
  compareResults(expected, df, "Function \"ml12(x+y)\"  (alias)  (x,y) = (0.5, 0.5)");
}

using capd::autodiff::Node;

void _ml12(Node /*t*/, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node /*params*/[], int /*noParams*/)
{
  out[0] = mittagLeffler12(in[0] + in[1]);
}

BOOST_AUTO_TEST_CASE(xmittagLeffler12node)
{
  std::string msg = "Function \"mittagLeffler12(x+y)\" [node]  (x,y) = ";
  MapType f(_ml12, 2, 1, 0, 5);
  VectorType x(2);
  JetType df(1, 2, 5);

  x[0] = 0.5;  x[1] = 0.5;
  std::vector<double> expected = computeML12Der(x);
  f(x, df);
  compareResults(expected, df, msg + "(0.5, 0.5)");

  x[0] = 1.5;  x[1] = 0.5;
  expected = computeML12Der(x);
  f(x, df);
  compareResults(expected, df, msg + "(1.5, 0.5)");

  x[0] = 0.5;  x[1] = -1.5;
  expected = computeML12Der(x);
  f(x, df);
  compareResults(expected, df, msg + "(0.5, -1.5)");
}

BOOST_AUTO_TEST_SUITE_END()
