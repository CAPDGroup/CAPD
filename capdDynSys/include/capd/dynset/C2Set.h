/////////////////////////////////////////////////////////////////////////////
/// @file C2Set.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#ifndef _CAPD_DYNSET_C2SET_H_ 
#define _CAPD_DYNSET_C2SET_H_ 

#include <string>

#include "capd/diffAlgebra/TimeRange.h"
#include "capd/dynsys/C2DynSys.h"
#include "capd/diffAlgebra/Hessian.h"
#include "capd/dynset/EnclosureHolder.h"
#include "capd/dynset/SetTraits.h"

namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/**
 * Common interface of all sets that store C2 information (set position and first and second derivatives)
 */
template<typename MatrixT>
class C2Set : public capd::diffAlgebra::TimeRange<typename MatrixT::ScalarType>,
              public capd::dynset::C2EnclosureHolder<MatrixT>{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef capd::diffAlgebra::Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;
  typedef C2Set SetType;
  typedef capd::dynsys::C2DynSys<MatrixType> DynSysType;

  /// constructor, initializes default enclosures and initial time
  C2Set(const VectorType& x, const VectorType& enc, const MatrixType& M, const MatrixType& mEnc, const HessianType& h, const HessianType& hEnc, const ScalarType& t)
      : capd::diffAlgebra::TimeRange<ScalarType>(t),
        capd::dynset::C2EnclosureHolder<MatrixType>(x,enc,M,mEnc,h,hEnc)
  {}
  virtual ~C2Set(){}

  virtual void move(DynSysType& c2dynsys) = 0;

  virtual const ScalarType& operator()(size_type i, size_type j, size_type c) const {
    return this->m_currentHessian(i,j,c);
  }
  virtual VectorType operator()(size_type j, size_type c) const{
    VectorType result(this->m_currentHessian.dimension());
    for(size_type i=0;i<this->m_currentHessian.dimension();i++)
      result[i] = this->m_currentHessian(i,j,c);
    return result;
  }
  virtual MatrixType operator()(size_type i){
    MatrixType result(this->m_currentHessian.dimension(),this->m_currentHessian.dimension());
    for(size_type j=0;j<this->m_currentHessian.dimension();j++)
      for(size_type c=0;c<this->m_currentHessian.dimension();c++)
        result(j+1,c+1) = this->m_currentHessian(i,j,c);

    return result;
  }

  const static size_type degree() { return 2; }
protected:
};


template<class MatrixT>
struct SetTraits< C2Set<MatrixT> >{
	const static bool isC0Set=false;
	const static bool isC1Set=false;
	const static bool isC2Set=true;
	const static bool isCnSet=false;
};

/// @}

}} // the end of the namespace capd::dynset

#endif // _CAPD_DYNSET_C2SET_H_ 

