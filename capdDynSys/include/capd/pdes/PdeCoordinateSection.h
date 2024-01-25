/// @addtogroup poincare
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file PdeCoordinateSection.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2016 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_PDES_PDE_COORDINATE_SECTION_H_
#define _CAPD_PDES_PDE_COORDINATE_SECTION_H_

#include <string>
#include <algorithm>
#include "capd/pdes/PdeAbstractSection.h"
#include "capd/poincare/CoordinateSection.h"

namespace capd{
namespace pdes{

/**
*  TimeMap class provides class that serves as Poincare section of the form x_i = c.
*  The section is defined by:
*  - integer - index of variable
*  - constant c that defines affine hyperplane x_i=c
*/

template<typename VectorT, typename MatrixT>
class PdeCoordinateSection : public PdeAbstractSection<VectorT,MatrixT>
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

  /// Constructs section \$ \{ x \in R^D | x_i = c \} (indices start at 0)\$
  PdeCoordinateSection(size_type D,         ///< phase space dimension
                       size_type _i,        ///< index of coordinate which defines section (counted from 0)
                       ScalarType _c = TypeTraits<ScalarType>::zero() ///< value that defines section
  )  : n(D,0.,2), c(_c), m_finiteSection(D,_i,_c)
  {
    setDirection(_i);
  }

  ScalarType operator()(const VectorType& v) const{
    return v[i]-c;
  }

  void setDirection(size_type i){
    if(i>=n.dimension())
       throw std::runtime_error("CoordinateSection::setDirection error - index of variable that defines PoincareSection must be less that dimension.");
    this->n.clear();
    this->i = i;
    this->n[i] = 1.;
    this->m_finiteSection.setDirection(i);
  }

  void setConstant(ScalarType _c){
    this->c = _c;
    this->m_finiteSection.setConstant(_c);
  }

  VectorType getNormalVector() const {
    return n;
  }

  const capd::poincare::AbstractSection<MatrixType>& getProjection(size_type) const {
    return this->m_finiteSection;
  }

  VectorType gradient(const VectorType&) const{
    return n;
  }

  ScalarType gradientByVector(const VectorType& /*x*/, const VectorType& u) const {
    return u[i];
  }

  bool isSpecialSection() const {
    return true;
  }

  ScalarType evalAt(const Set& s) const{
    return ((FiniteVectorType)s)[i]-c;
  }

private:
  capd::poincare::CoordinateSection<MatrixType> m_finiteSection;
  VectorType n;
  ScalarType c;
  size_type i;
}; // end of template PdeCoordinateSection

}} // namespace capd::poincare

#endif // _CAPD_PDES_PDE_COORDINATE_SECTION_H_

/// @}
