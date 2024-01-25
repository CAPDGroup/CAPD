//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file FadBadReversibleSolverTest.cpp
///
/// @author Daniel Wilczak
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD group
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#include "doTest.h"
#include "capd/dynsys/FadMap.h"
#include "capd/dynsys/FadOdeSolver.hpp"
#include "capd/poincare/TimeMap.hpp"

/**
 * The Michelson system is reversible with respect to involution (x,y,z,t) -> (-x,y,-z,-t)
 * We check rigorous C^0-C^1 solvers for preserving this property.
 */

class MichelsonFadMap : public capd::dynsys::FadMap<interval,0>
{
public:
  typedef interval Scalar;
  typedef capd::dynsys::FadMap<interval,0>::MatrixType MatrixType;
  typedef capd::dynsys::FadMap<interval,0>::VectorType VectorType;
  typedef capd::dynsys::FadFunction<VectorType> FunctionType;

  /// This operator computes vector field of the Michelson system.
  template <typename AVector>
  AVector operator()(const AVector& in) const
  {
    AVector out(3);
    out[0] = in[1];
    out[1] = in[2];
    out[2] = 1. - in[1] - 0.5*sqr(in[0]);
    return out;
  }

  /// This operator computes vector field of the Lorenz system.
  /// The system is autonomous, but the operator must have two arguments in order to fit requested signature.
  template <typename TimeT, typename AVector>
  AVector operator()(const TimeT& /*t*/, const AVector& in) const
  {
    return (*this)(in);
  }

  MatrixType derivative(const Scalar& t, const VectorType& u) const
  {
    MatrixType der(3,3);
    computeDerivative(*this,t,u,der);
    return der;
  }

  /// computes simultaneously value and derivative of the map for a given vector and time
  VectorType operator()(Scalar t, const VectorType& x, MatrixType& o_der) const
  {
    return computeDerivative(*this,t,x,o_der);
  }

  MatrixType operator[](const VectorType& u) const
  {
    MatrixType der(3,3);
    computeDerivative(*this,TypeTraits<Scalar>::zero(),u,der);
    return der;
  }

  unsigned dimension() const {return 3u;}
};

typedef capd::dynsys::FadOdeSolver<MichelsonFadMap> MichelsonSolver;
typedef capd::poincare::TimeMap<MichelsonSolver> TimeMap;
// ##################################################################


BOOST_AUTO_TEST_SUITE(FadSolverSuite)

BOOST_AUTO_TEST_CASE(C0Test)
{
  MichelsonFadMap f;
  MichelsonSolver solver(f,10);
  TimeMap tm(solver);
  IVector x(0.,1.563,0.);

  doC0test<C0Rect2Set>(tm,x);
  doC0test<C0Pped2Set>(tm,x);
  doC0test<C0TripletonSet>(tm,x);
//  doC0test<C0HORect2Set>(tm,x);
//  doC0test<C0HOTripletonSet>(tm,x);
}

BOOST_AUTO_TEST_CASE(C1Test)
{
  MichelsonFadMap f;
  MichelsonSolver solver(f,10);
  TimeMap tm(solver);
  IVector x(0.,1.563,0.);

  doC1test<C1Rect2Set>(tm,x);
  doC1test<C1Pped2Set>(tm,x);
//  doC1test<C1HORect2Set>(tm,x);
//  doC1test<C1HOPped2Set>(tm,x);
}

BOOST_AUTO_TEST_SUITE_END()
