/// @addtogroup diffAlgebra
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Hessian.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_HESSIAN_HPP_
#define _CAPD_DIFFALGEBRA_HESSIAN_HPP_

#include "capd/vectalg/Matrix.hpp"
#include "capd/diffAlgebra/Hessian.h"
#include "capd/vectalg/algebraicOperations.hpp"

namespace capd{
namespace diffAlgebra{


template<typename Scalar,__size_type M, __size_type D>
inline void setDimensionInArray(Hessian<Scalar,M,D>* , __size_type, __size_type, __size_type){
}

template<typename Scalar>
inline void setDimensionInArray(Hessian<Scalar,0,0>* t, __size_type N, __size_type M, __size_type D){
  for(__size_type i=0;i<N;++i) t[i].resize(M,D);
}

template<typename Scalar,__size_type M, __size_type D>
Hessian<Scalar,M,D>* Hessian<Scalar,M,D>::makeArray(size_type N, size_type _M, size_type _D)
{
  Hessian* result = new Hessian[N];
  setDimensionInArray(result,N,_M,_D);
  return result;
}

template<typename ScalarType, __size_type D>
Hessian<ScalarType,D,D> operator*(const capd::vectalg::Matrix<ScalarType,D,D>& m, const Hessian<ScalarType,D,D>& c2)
{
  typedef typename Hessian<ScalarType,D,D>::size_type size_type;
  const size_type dim = m.numberOfRows();
  if(dim==c2.dimension() and dim == c2.imageDimension())
  {
    Hessian<ScalarType,D,D> result(dim);
    for(size_type j=0;j<dim;++j)
      for(size_type c=j;c<dim;++c)
        for(size_type i=0;i<dim;++i)
          for(size_type r=0;r<dim;++r)
            result(i,j,c) += m[i][r]*c2(r,j,c);
    return result;
  }
  throw std::runtime_error("operator* - incompatible dimensions of matrix and Hessian");
}

template<typename ScalarType, __size_type D>
Hessian<ScalarType,D,D> operator*(const Hessian<ScalarType,D,D>& H, const capd::vectalg::Matrix<ScalarType,D,D>& M)
{
  typedef typename Hessian<ScalarType,D,D>::size_type size_type;
  const size_type dim = M.numberOfRows();
  if(dim==H.dimension() and dim==H.imageDimension())
  {
    Hessian<ScalarType,D,D> result(dim);
    size_type i, j, c, ss, r;
    for (j = 0; j < dim; ++j)
    {
      for (ss = 0; ss < dim; ++ss)
        for (r = ss; r < dim; ++r) {
          ScalarType temp = (r==ss) ? sqr(M(r+1,j+1)) : M(r+1,j+1) * M(ss+1,j+1);
          for (i = 0; i < dim; ++i)
            result(i, j, j) += H(i, ss, r) * temp;
        }
      for(c=j+1;c<dim;++c)
        for (ss = 0; ss < dim; ++ss)
          for (r = ss; r < dim; ++r) {
            ScalarType temp = (r==ss) ?
                    (2.*M(r+1,j+1))*M(r+1,c+1)
                  : M(r+1,j+1) * M(ss+1,c+1) + M(ss+1,j+1) * M(r+1,c+1);
            for (i = 0; i < dim; ++i)
              result(i, j, c) += H(i, ss, r) * temp;
          }
    }
    return result;
  }
  throw std::runtime_error("operator* - incompatible dimensions of Hessian and matrix");
}

// ---------------------------------------------------

template<typename ScalarType, __size_type M, __size_type D>
Hessian<ScalarType,M,D>& Hessian<ScalarType,M,D>::operator+=(const Hessian& s)
{
  if(m_domainDimension==s.m_domainDimension and m_imageDimension==s.m_imageDimension)
  {
    capd::vectalg::addObjects(*this,s,*this);
    return *this;
  }
  throw std::runtime_error("incompatible dimensions in Hessian::operator+=");
}

// ---------------------------------------------------

template<typename ScalarType, __size_type M, __size_type D>
Hessian<ScalarType,M,D>& Hessian<ScalarType,M,D>::operator*=(const ScalarType& s)
{
  return capd::vectalg::multiplyAssignObjectScalar(*this,s);
}

template<typename ScalarType, __size_type D>
Hessian<ScalarType,D,D> operator*(ScalarType c, const Hessian<ScalarType,D,D>& H){
  Hessian<ScalarType,D,D> result(H.dimension());
  capd::vectalg::multiplyObjectScalar(H,c,result);
  return result;
}

}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_C2COEFF_HPP_

/// @}
