//////////////////////////////////////////////////////////////////////////////
//   Package:          CAPD

/////////////////////////////////////////////////////////////////////////////
//
/// @file C0TimeJet.h
///
/// @author Tomasz Kapela   @date 2010-04-20
//
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) Tomasz Kapela 2010
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_C0TIMEJET_H_
#define _CAPD_DIFFALGEBRA_C0TIMEJET_H_

#include "capd/diffAlgebra/TimeRange.h"
#include "capd/diffAlgebra/CoeffTraits.h"

namespace capd{
namespace diffAlgebra{

template <typename VectorT>
class C0TimeJet : public capd::diffAlgebra::TimeRange<typename VectorT::ScalarType>{
public:
  typedef C0TimeJet SetType;
  typedef VectorT VectorType;
  typedef typename VectorType::ScalarType ScalarType;
  typedef typename VectorType::size_type size_type;
  typedef typename VectorType::iterator iterator;
  typedef typename VectorType::const_iterator const_iterator;

  explicit C0TimeJet(size_type dim)
    : capd::diffAlgebra::TimeRange<ScalarType>(TypeTraits<ScalarType>::zero() )
  {
    m_value = new VectorType(dim);
    is_owner=true;
  }

  explicit C0TimeJet(VectorType* v)
    : capd::diffAlgebra::TimeRange<ScalarType>(TypeTraits<ScalarType>::zero() ),
      m_value(v), is_owner(false)
  {}

  ~C0TimeJet(){
    if(is_owner) delete m_value;
  }

  void clear(){
    m_value->clear();
  }

  operator VectorType() const{
    return *m_value;
  }

  VectorType& vector(){
    return *m_value;
  }

  const VectorType& vector() const {
    return *m_value;
  }

  void setVector(const VectorType& v){
    *m_value = v;
  }

  friend std::ostream & operator<< (std::ostream & s, const C0TimeJet & coeff){
    s << coeff.getCurrentTime() << "\n" << *(coeff.m_value);
    return s;
  }

  size_type dimension() const{
	  return m_value->dimension();
  }

  template<class DynSysT>
  void move(DynSysT& dynsys){
    *m_value = dynsys(this->m_currentTime,*m_value);
  }

  static size_type degree(){
    return 0;
  }
protected:
   VectorType* m_value;
   bool is_owner;
};

template<class VectorT>
struct CoeffTraits<C0TimeJet<VectorT> >{
	const static bool isC0Jet=true;
	const static bool isC1Jet=false;
	const static bool isC2Jet=false;
	const static bool isCnJet=false;
};

}} // end of namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_C0TIMEJET_H_
