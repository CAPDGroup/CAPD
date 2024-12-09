/// @addtogroup vectalg
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file InvariantManifold.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_INVSET_INVARIANTMANIFOLD_H_
#define _CAPD_INVSET_INVARIANTMANIFOLD_H_

#include <stdexcept>
#include <vector>
#include "capd/capdlib.h"
#include "capd/matrixAlgorithms/capd2alglib.h"

namespace capd{
namespace invset{

/**
 * Given a map 'f' it computes the derivative of 'period' iteration
 * of 'f' at 'u'.
*/
template<class MapT, class V, class M>
V computeDerivative(MapT& f, V u, M& A, int period){
  A.setToIdentity();
  M temp(u.dimension(),u.dimension());
  for(int i=0;i<period;++i)
  {
    u = f(u,temp);
    A = temp*A;
  }
  return u;
}

template<class M>
std::pair<typename M::RowVectorType,typename M::RowVectorType> computeCoordSystem(const M& A, M& rVec)
{
  int dim = A.numberOfRows();
  typename M::RowVectorType rV(dim), iV(dim);  // real and imaginary parts of eigenvalues
  M iVec(dim,dim);                             // real and imaginary parts of eigenvectors
  typedef typename M::ScalarType Scalar;

  capd::alglib::computeEigenvaluesAndEigenvectors(A,rV,iV,rVec,iVec);

  // sort eigenvalues
  for(int j=0;j<dim-1;++j){
    for(int k=j+1;k<dim;++k){
      if(fabs(toDouble(rV[j])) < fabs(toDouble(rV[k]))){
        Scalar t = rV[j];
        rV[j] = rV[k];
        rV[k] = t;
        t = iV[j];
        iV[j] = iV[k];
        iV[k] = t;
        for(int s=0;s<dim;++s){
          t = rVec[s][j];
          rVec[s][j] = rVec[s][k];
          rVec[s][k] = t;
          t = iVec[s][j];
          iVec[s][j] = iVec[s][k];
          iVec[s][k] = t;
        }
      }
    }
  }

  // if eigenvalue is complex then take real and imag parts as eigenvectors
  for(int i=0;i<dim;++i){
    if(iV[i]!=0){
      rVec.column(i) = iVec.column(i);
      ++i;
    }
  }
  return std::make_pair<rV,iV>;
}

template<class MapT, class V, class M>
std::pair<V,V> computeCoordSystem(MapT& f, V u, int period, M& rVec){
  M A(u.dimension(),u.dimension(),false);
  computeDerivative(f,u,A,period);
  return computeCoordSystem(A,rVec);
}

/**
 * This is a generic algorithm for computing (nonrigorous) parameterization
 * of one-dimensional invariant manifold at a fixed point x.
 */
template<class Map, class V, class Jet>
void oneDimInvariantManifold(Map& f, typename Map::VectorType& x, Jet& jet, int period){
  const int d = x.dimension();
  typedef typename Jet::MatrixType M;
  Jet J(d,jet.degree()), temp(d,jet.degree());
  J() = x;
  J.setMatrix(M::Identity(d));
  int i;

  for(i=0;i<period;++i) J = f(J);
  M A = (M)J;
  M C(d,d,false);
  std::pair<V,V> computeCoordSystem(A,C);
  jet.clear();
  // set unstable direction
  for(i=0;i<d;++i)
    jet(i,0) = C[i][0];
  for(i=0;i<jet.degree();++i){
    substitutionPowerSeries(J,jet,temp);
    jet = temp;
  }
}

}}
#endif // _CAPD_INVSET_INVARIANTMANIFOLD_H_

/// @}
