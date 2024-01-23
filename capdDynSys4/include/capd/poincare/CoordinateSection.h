/////////////////////////////////////////////////////////////////////////////
/// @file CoordinateSection.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2013 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_COORDINATE_SECTION_H_
#define _CAPD_POINCARE_COORDINATE_SECTION_H_

#include <string>
#include <algorithm>
#include "capd/poincare/AbstractSection.h"
#include "capd/autodiff/NodeType.h"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
/**
*  TimeMap class provides class that serves as Poincare section of the form x_i = c.
*  The section is defined by:
*  - integer - index of variable
*  - constant c that defines affine hyperplane x_i=c
*/

template<typename MatrixT>
class CoordinateSection : public AbstractSection<MatrixT>
{
public:
  typedef MatrixT                   MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename VectorType::size_type size_type;             ///< integral type used to index containers (vectors, matrices, etc)
  typedef capd::dynset::AbstractSet<VectorType> Set;   ///< type of abstract base class for all sets
  typedef typename  AbstractSection<MatrixT>::JetType JetType;
  
  /// Constructs section \$ \{ x \in R^D | x_i = c \} (indices start at 0)\$
  CoordinateSection(size_type D,         ///< phase space dimension
                    size_type _i,        ///< index of coordinate which defines section (counted from 0)
                    ScalarType _c = TypeTraits<ScalarType>::zero() ///< value that defines section
  )  : n(D), c(_c) {
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
  }

  void setConstant(ScalarType _c){
    this->c = _c;
  }

  VectorType getNormalVector() const {
    return n;
  }

  VectorType gradient(const VectorType&) const{
    return n;
  }
  ScalarType gradientByVector(const VectorType& /*x*/, const VectorType& u) const {
    return u[this->i];
  }

  bool isSpecialSection() const {
    return true;
  }

  ScalarType evalAt(const capd::dynset::AbstractSet<VectorType>& s) const{
    return ((VectorType)s)[i]-c;
  }

  // -------------------------------------------------------------------

  using AbstractSection<MatrixT>::computeDT; // virtual function
  void computeDT(const JetType& Px, const JetType& vfOnPx, JetType& dT, size_type degree) const
  {
    const size_type dim = Px.dimension();
    // here we apply division lemma. The code is already implemented in evaluation of Acos node for automatic differentiation
    using capd::vectalg::Multiindex;
    const ScalarType* left = Px.begin(i);
    const ScalarType* right = vfOnPx.begin(i);
    ScalarType* result = dT.begin();
    Multiindex a(dim),b(dim),c(dim);
    unsigned j,p;
    const unsigned shift = binomial(dim+degree-1,dim);
    for(j=1;j<degree;++j){
      a.clear();
      a[0]=j;
      const unsigned shiftA = binomial(dim+j-1,dim);
      const unsigned shiftB = binomial(dim+(degree-j)-1,dim);
      do{
        b.clear();
        b[0]=degree-j;
        const unsigned sa = shiftA + a.index(j);
        do{
          p = capd::autodiff::sumAndFindMax(a.begin(),b.begin(),c.begin(),dim);
          const unsigned sb = shiftB + b.index(degree-j);
          result[shift + c.index(degree)] += a[p]*result[sa]*right[sb];
        }while(b.hasNext());
      }while(a.hasNext());
    }

    c.clear();
    c[0] = degree;
    do{
      p = capd::autodiff::findMax(c.begin(),dim);
      const unsigned s = shift + c.index(degree);
      result[s] = -(left[s] + result[s]/double(c[p]))/(*right);
    }while(c.hasNext());
  }

private:
  VectorType n;
  ScalarType c;
  size_type i;
}; // end of template CoordinateSection

/// @}
}} // namespace capd::poincare

#endif // _CAPD_POINCARE_COORDINATE_SECTION_H_

