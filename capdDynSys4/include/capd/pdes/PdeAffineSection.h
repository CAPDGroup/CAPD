/// @addtogroup pdes
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file PdeAffineSection.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2016 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_PDES_PDE_AFFINE_SECTION_H_
#define _CAPD_PDES_PDE_AFFINE_SECTION_H_

#include <string>
#include "capd/poincare/AffineSection.h"
#include "capd/pdes/PdeAbstractSection.h"

namespace capd{
namespace pdes{

/**
*  TimeMap class provides class that serves as an affine Poincare section.
*  The section is defined by two vectors - new origin of coordinate system and normal vector to hyperplane.
*/

template<typename VectorT, typename MatrixT>
class PdeAffineSection : public PdeAbstractSection<VectorT,MatrixT>
{
public:
  typedef MatrixT MatrixType;
  typedef VectorT VectorType;
  typedef typename MatrixType::RowVectorType FiniteVectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename VectorType::size_type size_type;             ///< integral type used to index containers (vectors, matrices, etc)
  typedef capd::dynset::AbstractSet<FiniteVectorType> Set;   ///< type of abstract base class for all sets

  typedef capd::diffAlgebra::Hessian<ScalarType,FiniteVectorType::csDim,FiniteVectorType::csDim> HessianType;
  typedef capd::diffAlgebra::Jet<MatrixT,0> JetType;

  PdeAffineSection(const VectorType& x, const VectorType& n)
    : m_finiteSection(x.getExplicitCoefficients(),n.getExplicitCoefficients()),
      m_projectionSection(x.getExplicitCoefficients(),n.getExplicitCoefficients())
  {}

  PdeAffineSection(const FiniteVectorType& x, const FiniteVectorType& n)
    : m_finiteSection(x,n), m_projectionSection(x,n)
  {}

  ScalarType operator()(const VectorType& v) const{
    return m_finiteSection(v.getExplicitCoefficients());
  }

  void setNormalVector(const VectorType& n){
    m_finiteSection.setNormalVector(n.getExplicitCoefficients());
  }

  void setOrigin(const VectorType& x){
    m_finiteSection.setOrigin(x.getExplicitCoefficients());
  }
/*
  VectorType getOrigin() const {
    return x;
  }

  VectorType getNormalVector() const {
    return n;
  }
*/
  VectorType gradient(const VectorType& u) const{
    return VectorType(0.,u.getGeometricDecay(),m_finiteSection.gradient(u.getExplicitCoefficients()));
  }

  ScalarType gradientByVector(const VectorType& x, const VectorType& u) const{
    return m_finiteSection.gradient(x.getExplicitCoefficients())*u.getExplicitCoefficients();
  }

  ScalarType evalAt(const Set& s) const{
    return s.evalAt(m_finiteSection);
  }

  const capd::poincare::AbstractSection<MatrixType>& getProjection(size_type dim) const {
    m_projectionSection = capd::poincare::AffineSection<MatrixType>(
      FiniteVectorType(dim,m_finiteSection.getOrigin().begin()),
      FiniteVectorType(dim,m_finiteSection.getNormalVector().begin())
    );
    return this->m_projectionSection;
  }

private:
  capd::poincare::AffineSection<MatrixType> m_finiteSection;
  mutable capd::poincare::AffineSection<MatrixType> m_projectionSection;
}; // end of template AffineSection

}} // namespace capd::poincare

#endif // _CAPD_POINCARE_AFFINE_SECTION_H_

/// @}
