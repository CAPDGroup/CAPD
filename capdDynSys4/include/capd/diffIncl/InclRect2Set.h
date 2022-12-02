/// @addtogroup diffIncl
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file InclRect2Set.h
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2007 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFINCL_INCLRECT2SET_H_ 
#define _CAPD_DIFFINCL_INCLRECT2SET_H_ 

#include "capd/dynset/lib.h"
#include "capd/dynset/C0DoubletonSet.h"
#include "capd/vectalg/Norm.h"

namespace capd{
namespace diffIncl{

///////////////////////////////////////////////////////////////////////////
// InclRect2Set
/// 
/// Set representation for differential inclusions based on capd::dynset::Rect2Set class
///
///   set is represented as: x + C*r0 + B*r   where
///       C*r0 - basic 'Lipschitz part'
///       B*r - QR-decomposition of the remaining errors
///   

template<typename MatrixT>
class InclRect2Set : public capd::dynset::C0DoubletonSet<MatrixT, capd::C0Rect2Policies> {

public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;
  typedef capd::dynset::C0DoubletonSet<MatrixT, capd::C0Rect2Policies> BaseSet;

  // constructors
  explicit InclRect2Set(const VectorType& the_x);
  InclRect2Set(const VectorType& the_x, const VectorType& the_r0);
  InclRect2Set(const VectorType& the_x, const MatrixType& the_C, const VectorType& the_r0);
  InclRect2Set(const VectorType& the_x, const MatrixType& the_C,
               const VectorType& the_r0,
               const VectorType& the_r
  );
  
  
  template<typename DiffIncl>
  void move( DiffIncl& dynsys);
  template<typename DiffIncl>
  void move( DiffIncl & dynsys, InclRect2Set& result) const;
  
  std::vector<VectorType> getCorners() const;

  using BaseSet::get_x;
  using BaseSet::get_r;
  using BaseSet::get_r0;
  using BaseSet::get_B;
  using BaseSet::get_C;
  using BaseSet::operator VectorType;
  using BaseSet::show;
  using BaseSet::affineTransformation;

protected:
  using BaseSet::m_x;
  using BaseSet::m_r;
  using BaseSet::m_r0;
  using BaseSet::m_B;
  using BaseSet::m_C;
};

template<typename MatrixType>
std::vector<typename MatrixType::VectorType> getCorners(const InclRect2Set<MatrixType> & set) ;

// inline definitions
////////////////////////////////////////////////////////////////////////////////
/// Constructors

template<typename MatrixType>
inline InclRect2Set<MatrixType>::InclRect2Set(const VectorType& the_x)
  :  BaseSet(the_x) {
}

template<typename MatrixType>
inline InclRect2Set<MatrixType>::InclRect2Set(const VectorType& the_x,const VectorType& the_r0)
  :  BaseSet(the_x, the_r0) {
}

template<typename MatrixType>
inline InclRect2Set<MatrixType>::InclRect2Set(
      const VectorType& the_x,
      const MatrixType& the_C,
      const VectorType& the_r0
   )
  : BaseSet(the_x, the_C, the_r0){
}

template<typename MatrixType>
inline InclRect2Set<MatrixType>::InclRect2Set(
      const VectorType& the_x,
      const MatrixType &the_C,
      const VectorType& the_r0,
      const VectorType& the_r
   ): BaseSet(the_x, the_C, the_r0, the_r){
}

////////////////////////////////////////////////////////////////////////////////
/// C0set interface overriding
   
////////////////////////////////////////////////////////////////////////////////
template<typename MatrixType>
std::vector<typename MatrixType::VectorType> getCorners(const capd::diffIncl::InclRect2Set<MatrixType> & set) ;

template<typename MatrixType>
template<typename DiffIncl>
void InclRect2Set<MatrixType>::move(DiffIncl & diffIncl) {

  VectorType x = VectorType(*this);

  //VectorType x = m_x + m_C * m_r0 + m_B * m_r; //VectorType(*this);
  //this->m_currentSet = x;
  ScalarType currentTime = this->getCurrentTime();
  //computation of an unperturbed trajectory
  BaseSet::move(diffIncl.getDynamicalSystem());

  //computation of the influence of the perturbations
  VectorType Deltha = diffIncl.perturbations(currentTime, x);

  // Rearrangements
  x = midVector( m_x + Deltha);
  VectorType dr = m_x - x + Deltha;

  m_r = m_r + this->get_invB() * dr;
  m_x = x;

  //this->m_currentSet += m_x + m_C * m_r0 + m_B * m_r;
  this->m_currentSet += Deltha;
}

template<typename MatrixType>
template<typename DiffIncl>
void InclRect2Set<MatrixType>::move(DiffIncl & diffIncl, InclRect2Set<MatrixType>& result) const {
  VectorType x = VectorType(*this);
  ScalarType currentTime = this->getCurrentTime();

  //computation of an unperturbed trajectory
  BaseSet::move(diffIncl.getDynamicalSystem(), result);

  //computation of the influence of the perturbations
  VectorType Deltha = diffIncl.perturbations(currentTime, x);

  // Rearrangements
  x = capd::vectalg::midVector( result.m_x + Deltha);
  VectorType dr = result.m_x + Deltha - x;

  result.m_r = result.m_r + result.get_invB() * dr;
  result.m_x = x;

//      result.m_currentSet = result.m_x + result.m_C * result.m_r0 + result.m_B * result.m_r;
  result.m_currentSet += Deltha;
}

}} // namespace capd::diffIncl

#endif // _CAPD_DIFFINCL_INCLRECT2SET_H_ 

/// @}
