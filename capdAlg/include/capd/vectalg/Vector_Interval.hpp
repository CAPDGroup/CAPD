/////////////////////////////////////////////////////////////////////////////
/// @file Vector_Interval.hpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_VECTALG_IVECTOR_HPP_
#define _CAPD_VECTALG_IVECTOR_HPP_

#include <stdexcept>
#include "capd/vectalg/iobject.hpp"
#include "capd/vectalg/EmptyIntersectionException.h"
#include "Vector.h"

namespace capd{
namespace vectalg{

// ---------------------------------------------------------------- //
// general algorithms for vectors and matrices

template<typename IVector>
inline
IVector leftVector(const IVector& u)
{
  return leftObject<IVector>(u);
}

template<typename IVector>
inline
IVector rightVector(const IVector& u)
{
  return rightObject<IVector>(u);
}

template<typename IVector>
typename IVector::ScalarType maxDiam(const IVector&);

template<typename IVector>
IVector intervalBall(const IVector&, const typename IVector::ScalarType &r);

template<typename IVector>
typename IVector::ScalarType solveAffineInclusion(const IVector&,const IVector&,const IVector&);

template<typename IVector>
typename IVector::ScalarType solveAffineInclusion(const IVector&,const IVector&,const IVector&,int&);

template<typename IVector>
IVector intervalHull(const IVector&, const IVector&);

template<typename IVector>
IVector midVector(const IVector&);

template<typename IVector>
IVector diam(const IVector&);

template<typename IVector>
IVector intersection(const IVector&, const IVector&);



//---------------  definitions for class intervalVector only -------------------//

template<typename IVectorType>
IVectorType diam(const IVectorType& v)
{
  IVectorType result(v.dimension(),true);
  diameter(v,result); // defined in iobject.hpp
  return result;
}

/******************************************************************************/
// intersection of two interval vectors
// halts if the intersection is empty
/******************************************************************************/

template<typename IVectorType>
IVectorType intersection(const IVectorType &v1, const IVectorType &v2)
{
  IVectorType result(v1.dimension(),true);
  if(!(intersection(v1,v2,result))) // defined in iobject.hpp
      throw EmptyIntersectionException("Intersection of two interval vectors is empty");
  return result;
}

template<typename IVectorType>
IVectorType midVector(const IVectorType& v)
{
  IVectorType result(v.dimension(),true);
  mid(v,result); // defined in iobject.hpp
  return result;
}

template<typename IVectorType>
IVectorType intervalHull(const IVectorType& v1, const IVectorType& v2)
{
  IVectorType result(v1.dimension(),true);
  intervalHull(v1,v2,result);
  return result;
}

template<typename IVectorType>
IVectorType intervalBall(const IVectorType &iv, const typename IVectorType::ScalarType &r)
{
  IVectorType result(iv.dimension(),true);
  typename IVectorType::iterator b=result.begin(), e=result.end();
  typename IVectorType::const_iterator i=iv.begin();

  while(b!=e)
  {
    *b = ball(*i,r);
    ++b;
    ++i;
  }
  return result;
}

template<typename IVectorType>
typename IVectorType::ScalarType solveAffineInclusion(
      const IVectorType& a,
      const IVectorType& p,
      const IVectorType& c
   )
{
  typedef typename IVectorType::ScalarType ScalarType;

  if(a.dimension()!=p.dimension() || a.dimension()!=c.dimension())
    throw std::runtime_error("Incompatible vectors in function solveAffineInclusion");

  ScalarType result = solveAffineInclusion(a[0],p[0],c[0]);
  for(typename IVectorType::size_type i=1; i<a.dimension(); ++i)
  {
    ScalarType iv=solveAffineInclusion(a[i],p[i],c[i]);
    result = capd::min( result.leftBound(), iv.leftBound() );
  }
  return result;
}

template<typename IVectorType>
typename IVectorType::ScalarType solveAffineInclusion(
      const IVectorType& a,
      const IVectorType& p,
      const IVectorType& c,
      int& dir
   )
{
  if(a.dimension()!=p.dimension() || a.dimension()!=c.dimension())
    throw std::runtime_error("Incompatible vectors in function solveAffineInclusion");

  typedef typename IVectorType::ScalarType ScalarType;
  typedef typename IVectorType::size_type size_type;
  ScalarType result = solveAffineInclusion(a[0],p[0],c[0]);
  dir=0;

  for(size_type i=1; i<a.dimension(); ++i)
  {
    ScalarType iv = solveAffineInclusion(a[i],p[i],c[i]);
    if( iv.leftBound() < result.leftBound() )
    {
      result = iv.leftBound();
      dir = i;
    }
  }
  return result;
}

}} // namespace capd::vectalg

#endif // _CAPD_VECTALG_IVECTOR_HPP_
