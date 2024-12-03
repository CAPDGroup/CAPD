/////////////////////////////////////////////////////////////////////////////
/// @file C11Rect2.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

/* Author: Daniel Wilczak, 2006 */

#ifndef _CAPD_DYNSET_C11RECT2_H_ 
#define _CAPD_DYNSET_C11RECT2_H_ 

#include "capd/dynset/C1Set.h"
#include "capd/vectalg/Norm.h"
#include "capd/dynsys/C1DynSys.h"
#include "capd/dynset/QRPolicy.h"


namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/**
 *  \f$C^1\f$ doubleton set with reorganization moved by QR decomposition (3rd method).
 *
 *  C^1-Lohner algorithm.
 *  Derivative of the flow moved via QR decomposition (third method)
 *  the set part - QR decomposition (3-rd method)
*/
template <typename MatrixT, typename QRPolicy = FullQRWithPivoting>
class C11Rect2 : public C1Set<MatrixT>{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  C11Rect2(const VectorType& x, double s_factor = 20);
  C11Rect2(const VectorType& x, const VectorType& r0, double sizeFactor = 20);
  C11Rect2(const VectorType& x, const MatrixType& C, const VectorType& r0, double sizeFactor = 20);
  virtual ~C11Rect2(){}

  void move(capd::dynsys::C1DynSys<MatrixType>& c1dynsys);
  void move(capd::dynsys::C1DynSys<MatrixType>& c1dynsys, C11Rect2& result) const;

  std::string show(void) const;
  double setFactor(double sFactor); // sets size factor :  0 = no swapping
  void disableSwapping();
  double getFactor();

  operator VectorType() const;
  operator MatrixType() const;
  C11Rect2 &operator=(const VectorType &v);

  /// This method computes value of functor f at interval vector represented by this set.
  /// It computes the value as an intersection of two intervals:
  /// 1. f(set)
  /// 2. f(center(set)) + Df(set)*(set-center(set))
  /// This set is represented as doubleton X=x+C*r0+B*r. The above reduces to
  /// f(x) + (Df(X)*C)*r0 + (Df(X)*B)*r
  template<class Functional>
  ScalarType evalAt(const Functional& f) const {
    ScalarType r = f(this->m_x);
    VectorType D = f.gradient((VectorType)(*this));
    for(size_type i=0;i<this->m_x.dimension();++i){
      r += (D*this->m_B.column(i))*this->m_r[i];
      r += (D*this->m_C.column(i))*this->m_r0[i];
    }
    if(!intersection(r, f(this->m_currentSet), r)){
      throw std::logic_error("C11Rect2::evalAt - empty intersection!. Report this error to CAPD developers!");
    }    return r;
  }
  protected:

// a representation of the set: x + C*r0 + B*r

  VectorType m_x,m_r,m_r0;
  MatrixType m_B,m_C;

  MatrixType m_D, m_Bjac, m_R, m_Cjac, m_R0; // D + Bjac*R  + Cjac * R0;

  double m_sizeFactor;  // if size of r0 is greater than r *size_factor we  put r in r0 and set r=0
                       // if size_factor = 0 we disable swapping
};
/// @}

// -------------------------------------------------------------------

template<typename MatrixType,typename QRPolicy>
inline void C11Rect2<MatrixType,QRPolicy>::move(capd::dynsys::C1DynSys<MatrixType>& c1dynsys)
{
  move(c1dynsys,*this);
}

template<typename MatrixType,typename QRPolicy>
inline C11Rect2<MatrixType,QRPolicy>::operator typename C11Rect2<MatrixType,QRPolicy>::VectorType(void) const
{
  return m_x + m_C*m_r0 + m_B*m_r;
}

template<typename MatrixType,typename QRPolicy>
inline double C11Rect2<MatrixType,QRPolicy>::setFactor(double sFactor) // sets size factor :  0 = no swapping
{
  double s = m_sizeFactor;
  m_sizeFactor = sFactor;
  return s;
}

template<typename MatrixType,typename QRPolicy>
inline void C11Rect2<MatrixType,QRPolicy>::disableSwapping()
{
  m_sizeFactor = 0.0;
}

template<typename MatrixType,typename QRPolicy>
inline double C11Rect2<MatrixType,QRPolicy>::getFactor() // returns size factor
{
  return m_sizeFactor;
}

}} // the end of the namespace capd::dynset

#endif // _CAPD_DYNSET_C11RECT2_H_ 

