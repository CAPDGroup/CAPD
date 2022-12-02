/// @addtogroup pdes
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file GeometricBound.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008-2016 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#include "capd/pdes/GeometricBound.h"
#include "capd/vectalg/Vector.hpp"
namespace capd {
namespace pdes {

// ---------------------------- constructors ----------------------------------

GeometricBound::GeometricBound(size_type dim)
  : m_x(dim), m_C(0.), m_decay(2.)
{
  this->check();
}

GeometricBound::GeometricBound(size_type dim, ScalarType C, ScalarType decay)
  : m_x(dim), m_C(C.rightBound()), m_decay(decay.leftBound())
{
  this->check();
}

GeometricBound::GeometricBound(size_type dim, ScalarType C, ScalarType q, const ScalarType* coeff)
  : m_x(dim,coeff), m_C(C.rightBound()), m_decay(q.leftBound())
{
  this->check();
}

GeometricBound::GeometricBound(ScalarType C, ScalarType decay, const FiniteVectorType& u)
  : m_x(u), m_C(C.rightBound()), m_decay(decay.leftBound())
{
  this->check();
}

void GeometricBound::check(){
  if(! (m_decay>1.) )
  {
    throw std::runtime_error("GeometricBound::check - object in an inconsistent state (q<=1)");
  }
}

// --------------------------- operator+ --------------------------------------

GeometricBound operator+(const GeometricBound& x, const GeometricBound& y) {
  return GeometricBound(
      x.getConstant()+y.getConstant(),
      capd::min(x.getGeometricDecay().leftBound(),y.getGeometricDecay().leftBound()),
      x.getExplicitCoefficients() + y.getExplicitCoefficients()
  );
}

GeometricBound& GeometricBound::operator+=(const GeometricBound& y) {
  this->m_C += y.getConstant();
  this->m_decay = capd::min(this->getGeometricDecay().leftBound(),y.getGeometricDecay().leftBound());
  this->m_x += y.getExplicitCoefficients();
  return *this;
}

// --------------------------- operator- --------------------------------------
GeometricBound operator-(const GeometricBound& x, const GeometricBound& y) {
  return GeometricBound(
      x.getConstant()+y.getConstant(),
      capd::min(x.getGeometricDecay().leftBound(),y.getGeometricDecay().leftBound()),
      x.getExplicitCoefficients() - y.getExplicitCoefficients()
  );
}

GeometricBound& GeometricBound::operator-=(const GeometricBound& y) {
  this->m_C += y.getConstant();
  this->m_decay = capd::min(this->getGeometricDecay().leftBound(),y.getGeometricDecay().leftBound());
  this->m_x -= y.getExplicitCoefficients();
  return *this;
}

// ----------------------------------------------------------------------------

GeometricBound operator*(const interval& s, const GeometricBound& x) {
  return GeometricBound(
      x.getConstant()*s,
      x.getGeometricDecay(),
      s*x.getExplicitCoefficients()
  );
}

GeometricBound operator*(const GeometricBound& x, const interval& s) {
  return s*x;
}

GeometricBound& GeometricBound::operator*=(const ScalarType& s) {
  this->m_C *= s;
  this->m_x *= s;
  return *this;
}

// ----------------------------------------------------------------------------
/**
 * This operator realizes the following Matrix by vector multiplication
 * (A*mainVariables,Id*tail)
 */
GeometricBound operator*(const IMatrix& A, const GeometricBound& v) {
  if (A.numberOfRows() > v.dimension() || A.numberOfColumns() > v.dimension())
    throw std::runtime_error("Cannot multiply matrix by GeometricBound - incorrect dimensions");

  GeometricBound::FiniteVectorType x = A*GeometricBound::FiniteVectorType(A.numberOfColumns(),v.getExplicitCoefficients().begin());
  GeometricBound result = v;
  result.projection(A.numberOfRows()) = A*v.projection(A.numberOfColumns());
  return result;
}

// ----------------------------------------------------------------------------

std::ostream& operator << (std::ostream& out, const GeometricBound& x)
{
  out << "{"
      << x.getExplicitCoefficients() << ","
      << x.getConstant() << ","
      << x.getGeometricDecay() << "}";
  return out;
}

// ----------------------------------------------------------------------------

void split(const GeometricBound& X, GeometricBound& x, GeometricBound& dx){
  x.setConstant(0.);
  x.setGeometricDecay(X.getGeometricDecay());
  dx.setConstant(X.getConstant());
  dx.setGeometricDecay(X.getGeometricDecay());
  split(X.getExplicitCoefficients(),x.getExplicitCoefficients(),dx.getExplicitCoefficients());
}

// ----------------------------------------------------------------------------

GeometricBound intersection(const GeometricBound& x, const GeometricBound& y){
  double qX = x.getGeometricDecay().leftBound();
  double qY = y.getGeometricDecay().leftBound();
  double cX = x.getConstant().rightBound();
  double cY = y.getConstant().rightBound();

  // if the same decay, take minimum of constantsx cX and cY
  if(qX==qY)
    return GeometricBound(
        capd::min(cX,cY),
        qX,
        intersection(x.getExplicitCoefficients(),y.getExplicitCoefficients())
    );
  // otherwise, choose the tail with faster decay as the intersection
  if(qX>qY)
      return GeometricBound(
          cX,
          qX,
          intersection(x.getExplicitCoefficients(),y.getExplicitCoefficients())
      );
  return GeometricBound(
      cY,
      qY,
      intersection(x.getExplicitCoefficients(),y.getExplicitCoefficients())
  );
}

} // end of namespace pdes
} // end of namespace capd

/// @}


