//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file C1TimeJet.h
///
/// @author Tomasz Kapela   @date 2010-04-20
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_C1TIMEJET_H_
#define _CAPD_DIFFALGEBRA_C1TIMEJET_H_

#include "capd/diffAlgebra/C0TimeJet.h"

namespace capd{
namespace diffAlgebra{

template <typename MatrixT>
class C1TimeJet : public capd::diffAlgebra::C0TimeJet<typename MatrixT::RowVectorType>{
public:
  typedef C1TimeJet SetType;
  typedef capd::diffAlgebra::C0TimeJet<typename MatrixT::RowVectorType> BaseCoeff;
  typedef MatrixT MatrixType;

  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::RefColumnVectorType RefVectorType;
  typedef typename MatrixType::RowVectorType VectorType;

  typedef typename VectorType::size_type size_type;
  typedef typename VectorType::iterator iterator;
  typedef typename VectorType::const_iterator const_iterator;

  explicit C1TimeJet(size_type dim)
    : BaseCoeff(dim)
  {
    m_derivative = new MatrixType(dim,dim);
  }

  C1TimeJet(VectorType* v, MatrixType* der)
    : BaseCoeff(v), m_derivative(der)
  {}

  ~C1TimeJet(){
    if(this->is_owner) delete m_derivative;
  }

  void clear(){
    BaseCoeff::clear();
    m_derivative->clear();
  }
  MatrixType& matrix(){
    return *m_derivative;
  }
  const MatrixType& matrix() const{
    return *m_derivative;
  }
  void setMatrix(const MatrixType& m){
    *m_derivative = m;
  }
  friend std::ostream & operator<< (std::ostream & s, const C1TimeJet & coeff){
    s << static_cast<const BaseCoeff&>(coeff);
    s << "\n" << *(coeff.m_derivative);
    return s;
  }

   template<class DynSysT>
   void move(DynSysT& dynsys){
	   *(this->m_value) = dynsys(this->m_currentTime,*(this->m_value),*(this->m_derivative),*(this->m_derivative));
   }

   static size_type degree(){
   	   return 1;
   }
protected:
   MatrixType* m_derivative;
};

template<class MatrixT>
struct CoeffTraits<C1TimeJet<MatrixT> >{
	const static bool isC0Jet=false;
	const static bool isC1Jet=true;
	const static bool isC2Jet=false;
	const static bool isCnJet=false;
};

}} // end of namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_C1TIMEJET_H_
