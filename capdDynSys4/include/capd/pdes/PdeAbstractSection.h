/////////////////////////////////////////////////////////////////////////////
/// @file PdeAbstractSection.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_PDES_ABSTRACT_PDE_SECTION_H_
#define _CAPD_PDES_ABSTRACT_PDE_SECTION_H_

#include <string>
#include "capd/dynset/AbstractSet.h"
#include "capd/diffAlgebra/Hessian.h"
#include "capd/diffAlgebra/Jet.h"
#include "capd/poincare/AbstractSection.h"
#include "capd/pdes/PdeSectionDerivativesEnclosure.h"

namespace capd{

namespace pdes{
/// @addtogroup pdes 
/// @{
template<class VectorT, class MatrixT>
class PdeAbstractSection
{
public:
  typedef MatrixT MatrixType;
  typedef VectorT VectorType;
  typedef typename MatrixType::RowVectorType FiniteVectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename VectorType::size_type size_type;             ///< integral type used to index containers (vectors, matrices, etc)
  typedef capd::dynset::AbstractSet<FiniteVectorType> Set;   ///< type of abstract base class for all sets
  typedef PdeSectionDerivativesEnclosure<VectorType,MatrixType> SectionDerivativesEnclosureType;

  typedef capd::diffAlgebra::Hessian<ScalarType,FiniteVectorType::csDim,FiniteVectorType::csDim> HessianType;
  typedef capd::diffAlgebra::Jet<MatrixT,0> JetType;

  virtual ScalarType operator() (const VectorType& v) const = 0;     ///< evaluates function at a given vector
  virtual VectorType gradient(const VectorType& u) const = 0;        ///< returns gradient of the function computed at vector u
  virtual ScalarType gradientByVector(const VectorType& x, const VectorType& u) const = 0;

  virtual const capd::poincare::AbstractSection<MatrixType>& getProjection(size_type) const = 0;
  /** This is very important function.
      If it returns true, class PoincareMap delegates computation of value of section(set) to the section.
      Otherwise it is assumed that the set has more information to compute value of section(set) in most optimal way.
      This is quite natural as the set knows its own representation.
      This function returns true for instance if the PoincareSection is given by x_i=c, where x_i is i-th coordinate and c is constant.
  */
  virtual bool isSpecialSection() const{
    return false;
  }

  /** This function computes value of section function on a given set.*/
  virtual ScalarType evalAt(const Set& set) const
  {
    throw std::logic_error("AbstractPdeSection::evalAt not implemented");
  }

  /** Simultaneous computation of gradient of return time and derivative of Poincare Map dP.
      @param[in] P - value of Poincare map
      @param[in] D - finite block of solution to first variational equation computed at return time
      @param[out] vfOnP - vector field evaluated on P
  */
  virtual MatrixType computeDP(
        const VectorType& P, const VectorType& vfOnP, const MatrixType& Dxx,
        const FiniteVectorType& Dy, const ScalarType& Dyx,
        FiniteVectorType& Py, ScalarType& Pyx) const
  {

    // first inner product of vector field and section gradient
    VectorType grad = this->gradient(P);
    ScalarType g = grad.getExplicitCoefficients()*vfOnP.getExplicitCoefficients();//this->gradientByVector(P,vfOnP);
    MatrixType DP = Dxx;

    size_type m = Dxx.numberOfRows();
    // Then T_{x_j} for j<= m and Pxx block which is given by an explicit formula
    FiniteVectorType dT(m);
    ScalarType T = 0.; // will contain |T_{x_1}| + ... + |T_{x_m}|
    ScalarType G = 0.; // will contain |grad_1|*V_{1y} + |grad_m|*V_{my}
    for (size_type j = 0; j < m; ++j){
      for(size_type i =0;i<m;++i)
        dT[j] -= grad[i]*Dxx(i+1,j+1);
      dT[j] /= g;
      //ScalarType dT = -(grad.getExplicitCoefficients()*D.column(j))/g;
      T += abs(dT[j]);
      G += (abs(grad[j])*Dy[j]).rightBound();
      for(size_type i = 0; i<m; ++i)
        DP(i+1,j+1) += vfOnP[i]*dT[j];
    }

    // Block Pyx that is j<=m x>m
    Pyx = Dyx + T*abs(vfOnP[m]);

    ScalarType C = abs(g)*G;
    for(size_type j=0;j<=m;++j){
      Py[j] = Dy[j] + C*abs(vfOnP[j]);
    }

    return DP;

  }
}; // end of template AbstractSection

/// @}
}} // namespace capd::poincare

#endif  /* _CAPD_POINCARE_ABSTRACT_SECTION_H_ */

