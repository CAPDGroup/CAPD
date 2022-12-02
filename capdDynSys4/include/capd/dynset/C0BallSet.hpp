
/////////////////////////////////////////////////////////////////////////////
/// @file C0BallSet.hpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0BALLSET_HPP_
#define _CAPD_DYNSET_C0BALLSET_HPP_

#include "capd/dynset/C0BallSet.h"
#include "capd/vectalg/Vector_Interval.hpp"

namespace capd{
namespace dynset{

template<typename MatrixType>
C0BallSet<MatrixType>::C0BallSet(const VectorType& a_x, const ScalarType& a_r, NormType& a_n, ScalarType t)
  : C0Set<MatrixType>(intervalBall(a_x,a_r),VectorType(a_x.dimension()),t),
    x(a_x), r(a_r), n(a_n.clone())
{}

template<typename MatrixType>
C0BallSet<MatrixType>::C0BallSet(const VectorType& a_x, NormType& a_n, ScalarType t)
  : C0Set<MatrixType>(a_x,VectorType(a_x.dimension()),t),
    x(a_x), r(0.), n(a_n.clone())
{}

template<typename MatrixType>
C0BallSet<MatrixType>::C0BallSet(const C0BallSet& a_set)
  : C0Set<MatrixType>(a_set),
    x(a_set.x), r(a_set.r), n(a_set.n->clone())
{}

template<typename MatrixType>
const C0BallSet<MatrixType>& C0BallSet<MatrixType>::operator=(const C0BallSet& a_set)
{
  C0Set<MatrixType>::operator=(a_set);
  x = a_set.x;
  r = a_set.r;
  delete n;
  n = a_set.n->clone();
  return *this;
}

template<typename MatrixType>
C0BallSet<MatrixType>::~C0BallSet()
{
  delete n;
}

template<typename MatrixType>
void C0BallSet<MatrixType>::move(capd::dynsys::DynSys<MatrixType> &dynsys)
{
  this->move(dynsys,*this);
}

template<typename MatrixType>
void C0BallSet<MatrixType>::move(capd::dynsys::DynSys<MatrixType> &dynsys, C0BallSet &result) const
{
  VectorType y(x.dimension()),s(x.dimension());

/*   
  Radius is multipled by the Lipschitz constant of the generator of DynSys
  on the IntervalType enclosure of the ball.
  Note that the proper computation of Lipschitz
  for flows is hidden in the implementation of Lipschitz for flows (class odenum)
*/
  result.r = r*dynsys.Lipschitz(this->getCurrentTime(),intervalBall(x,r),*n);

/* x is replaced by the center of its image under the generator */
  y = dynsys.Phi(this->getCurrentTime(),x) + dynsys.Remainder(this->getCurrentTime(),x,result.m_lastEnclosure);
  split(y,s);
  result.x = y;
/* and the norm of the arising deviation s is added to the radius */
  result.r += (*n)(s);
  result.r = result.r.rightBound();
  result.m_currentSet = intervalBall(result.x,result.r);
  result.setCurrentTime(this->getCurrentTime()+dynsys.getStep());
}

template<typename MatrixType>
std::string C0BallSet<MatrixType>::show(void) const
{
  std::ostringstream descriptor;
  descriptor << "C0BallSet(" << n->show() << "): x=";
  descriptor << x << " r=";
  descriptor << r << " ";
  return descriptor.str();
}

template<typename MatrixType>
typename C0BallSet<MatrixType>::VectorType C0BallSet<MatrixType>::affineTransformation(
    const MatrixType& A_M,
    const VectorType& A_C
  ) const
{
  return A_M*(intervalBall(x,r)-A_C);
}


}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C0BALLSET_HPP_
