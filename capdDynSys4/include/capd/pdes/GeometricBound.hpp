/// @addtogroup pdes
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file GeometricBound.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008-2022 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#include "capd/pdes/GeometricBound.h"
#include "capd/vectalg/Vector.hpp"

#ifndef _CAPD_PDES_GeometricBound_HPP_
#define _CAPD_PDES_GeometricBound_HPP_

namespace capd {
namespace pdes {

// ---------------------------- constructors ----------------------------------
template<class ScalarT>
GeometricBound<ScalarT>::GeometricBound(size_type dim)
  : m_x(dim), m_C(0.), m_decay(2.)
{
  this->check();
}

template<class ScalarT>
GeometricBound<ScalarT>::GeometricBound(size_type dim, capd::interval C, capd::interval decay)
  : m_x(dim), m_C(C.rightBound()), m_decay(decay.leftBound())
{
  this->check();
}

template<class ScalarT>
GeometricBound<ScalarT>::GeometricBound(size_type dim, capd::interval C, capd::interval q, const ScalarType* coeff)
  : m_x(dim,coeff), m_C(C.rightBound()), m_decay(q.leftBound())
{
  this->check();
}

template<class ScalarT>
GeometricBound<ScalarT>::GeometricBound(capd::interval C, capd::interval decay, const FiniteVectorType& u)
  : m_x(u), m_C(C.right()), m_decay(decay.left())
{
  this->check();
}

template<class ScalarT>
void GeometricBound<ScalarT>::check(){
  if( !(m_decay>1.) or !(m_C>=0) )
  {
    throw std::runtime_error("GeometricBound::check - object in an inconsistent state (q<=1) or C<0");
  }
}

// --------------------------- operator+ --------------------------------------

template<class ScalarT>
GeometricBound<ScalarT> operator+(const GeometricBound<ScalarT>& x, const GeometricBound<ScalarT>& y) {
  return GeometricBound<ScalarT>(
      x.getConstant()+y.getConstant(),
      capd::min(x.getGeometricDecay().leftBound(),y.getGeometricDecay().leftBound()),
      x.getExplicitCoefficients() + y.getExplicitCoefficients()
  );
}

template<class ScalarT>
GeometricBound<ScalarT>& GeometricBound<ScalarT>::operator+=(const GeometricBound& y) {
  this->m_C += y.getConstant();
  this->m_decay = capd::min(this->getGeometricDecay().leftBound(),y.getGeometricDecay().leftBound());
  this->m_x += y.getExplicitCoefficients();
  return *this;
}

// --------------------------- operator- --------------------------------------
template<class ScalarT>
GeometricBound<ScalarT> operator-(const GeometricBound<ScalarT>& x, const GeometricBound<ScalarT>& y) {
  return GeometricBound<ScalarT>(
      x.getConstant()+y.getConstant(),
      capd::min(x.getGeometricDecay().leftBound(),y.getGeometricDecay().leftBound()),
      x.getExplicitCoefficients() - y.getExplicitCoefficients()
  );
}

template<class ScalarT>
GeometricBound<ScalarT>& GeometricBound<ScalarT>::operator-=(const GeometricBound& y) {
  this->m_C += y.getConstant();
  this->m_decay = capd::min(this->getGeometricDecay().leftBound(),y.getGeometricDecay().leftBound());
  this->m_x -= y.getExplicitCoefficients();
  return *this;
}

// ----------------------------------------------------------------------------

template<class ScalarT>
GeometricBound<ScalarT> operator*(const ScalarT& s, const GeometricBound<ScalarT>& x) {
  return GeometricBound<ScalarT>(
      capd::abs(x.getConstant()*s).rightBound(),
      x.getGeometricDecay(),
      s*x.getExplicitCoefficients()
  );
}

template<class ScalarT>
GeometricBound<ScalarT> operator*(const GeometricBound<ScalarT>& x, const ScalarT& s) {
  return s*x;
}

template<class ScalarT>
GeometricBound<ScalarT>& GeometricBound<ScalarT>::operator*=(const ScalarType& s) {
  this->m_C = capd::abs(this->m_C*s).rightBound();
  this->m_x *= s;
  return *this;
}

// ----------------------------------------------------------------------------
/**
 * This operator realizes the following Matrix by vector multiplication
 * (A*mainVariables,Id*tail)
 */
template<class ScalarT>
GeometricBound<ScalarT> operator*(const capd::vectalg::Matrix<ScalarT,0,0>& A, const GeometricBound<ScalarT>& v) {
  if (A.numberOfRows() > v.dimension() || A.numberOfColumns() > v.dimension())
    throw std::runtime_error("Cannot multiply matrix by GeometricBound - incorrect dimensions");

  GeometricBound<ScalarT> result = v;
  result.projection(A.numberOfRows()) = A*v.projection(A.numberOfColumns());
  return result;
}

// ----------------------------------------------------------------------------

template<class ScalarT>
std::ostream& operator << (std::ostream& out, const GeometricBound<ScalarT>& x)
{
  out << "{"
      << x.getExplicitCoefficients() << ","
      << x.getConstant() << ","
      << x.getGeometricDecay() << "}";
  return out;
}

// ----------------------------------------------------------------------------

template<class ScalarT>
void split(const GeometricBound<ScalarT>& X, GeometricBound<ScalarT>& x, GeometricBound<ScalarT>& dx){
  x.setConstant(0.);
  x.setGeometricDecay(X.getGeometricDecay());
  dx.setConstant(X.getConstant());
  dx.setGeometricDecay(X.getGeometricDecay());
  split(X.getExplicitCoefficients(),x.getExplicitCoefficients(),dx.getExplicitCoefficients());
}

// ----------------------------------------------------------------------------

template<class ScalarT>
GeometricBound<ScalarT> intersection(const GeometricBound<ScalarT>& x, const GeometricBound<ScalarT>& y){
  double qX = x.getGeometricDecay().leftBound();
  double qY = y.getGeometricDecay().leftBound();
  double cX = x.getConstant().rightBound();
  double cY = y.getConstant().rightBound();

  // if the same decay, take minimum of constants cX and cY
  if(qX==qY)
    return GeometricBound<ScalarT>(
        capd::min(cX,cY),
        qX,
        intersection(x.getExplicitCoefficients(),y.getExplicitCoefficients())
    );
  // otherwise, choose the tail with faster decay as the intersection
    return GeometricBound<ScalarT>(
        capd::max(cX,cY),
        capd::max(qX,qY),
        intersection(x.getExplicitCoefficients(),y.getExplicitCoefficients())
    );
}

} // end of namespace pdes
} // end of namespace capd

#endif // _CAPD_PDES_GeometricBound_HPP_

/// @}


