/// @addtogroup diffAlgebra
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file C2TimeJet.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_C2TIMEJET_H_
#define _CAPD_DIFFALGEBRA_C2TIMEJET_H_

#include "capd/diffAlgebra/C1TimeJet.h"
#include "capd/diffAlgebra/Hessian.h"

namespace capd{
namespace diffAlgebra{

template<typename MatrixT>
class C2TimeJet : public C1TimeJet<MatrixT>
{
  public:
  typedef C2TimeJet SetType;
  typedef MatrixT MatrixType;
  typedef C1TimeJet<MatrixT> BaseCoeff;

  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::RefColumnVectorType RefVectorType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;
  typedef typename VectorType::size_type size_type;
  typedef typename VectorType::iterator iterator;
  typedef typename VectorType::const_iterator const_iterator;

  explicit C2TimeJet(size_type dim)
    : BaseCoeff(dim)
  {
    m_hessian = new HessianType(dim);
  }

  C2TimeJet(VectorType* v, MatrixType* der, HessianType* h)
    : BaseCoeff(v,der), m_hessian(h)
  {}

  ~C2TimeJet(){
    if(this->is_owner) delete m_hessian;
  }

  void clear(){
    BaseCoeff::clear();
    m_hessian->clear();
  }
  void setHessian(const HessianType& hessian){
    *m_hessian = hessian;
  }
  const HessianType& hessian()const{
    return *m_hessian;
  }
  HessianType& hessian(){
    return *m_hessian;
  }

  template<class DynSysT>
  void move(DynSysT& dynsys){
	   *(this->m_value) = dynsys(this->m_currentTime,*(this->m_value),*(this->m_derivative),*(this->m_hessian),*(this->m_derivative),*(this->m_hessian));
  }

  static size_type degree(){
    return 2;
  }
protected:
  HessianType* m_hessian;
}; // the end of class C2Coeff

template<class MatrixT>
struct CoeffTraits<C2TimeJet<MatrixT> >{
	const static bool isC0Jet=false;
	const static bool isC1Jet=false;
	const static bool isC2Jet=true;
	const static bool isCnJet=false;
};

}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_C2TIMEJET_H_

/// @}
