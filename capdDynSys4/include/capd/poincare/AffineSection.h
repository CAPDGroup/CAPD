

/////////////////////////////////////////////////////////////////////////////
/// @file AffineSection.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2013 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_AFFINE_SECTION_H_
#define _CAPD_POINCARE_AFFINE_SECTION_H_

#include <string>
#include "capd/poincare/AbstractSection.h"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
/**
*  TimeMap class provides class that serves as an affine Poincare section.
*  The section is defined by two vectors - new origin of coordinate system and normal vector to hyperplane.
*/

template<typename MatrixT>
class AffineSection : public AbstractSection<MatrixT>
{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename VectorType::size_type size_type;             ///< integral type used to index containers (vectors, matrices, etc)
  typedef capd::dynset::AbstractSet<VectorType> Set;   ///< type of abstract base class for all sets
  typedef typename  AbstractSection<MatrixT>::JetType JetType;

  AffineSection(const VectorType& _x, const VectorType& _n) : x(_x), n(_n), c(_n*_x) {}

  ScalarType operator()(const VectorType& v) const{
    return n*v-c;
  }

  void setNormalVector(const VectorType& n){
    this->n = n;
    this->c = n*x;
  }

  void setOrigin(const VectorType& x){
    this->x = x;
    this->c = n*x;
  }

  VectorType getOrigin() const {
    return x;
  }

  VectorType getNormalVector() const {
    return n;
  }

  VectorType gradient(const VectorType&) const{
    return n;
  }

  ScalarType evalAt(const capd::dynset::AbstractSet<VectorType>& s) const{
    return s.evalAt(*this);
  }

private:
  VectorType x, n;
  ScalarType c;
}; // end of template AffineSection
/// @}
}} // namespace capd::poincare

#endif // _CAPD_POINCARE_AFFINE_SECTION_H_


