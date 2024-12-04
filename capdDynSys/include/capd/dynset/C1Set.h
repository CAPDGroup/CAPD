////////////////////////////////////////////////////////////////////////////
/// @file C1Set.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C1SET_H_
#define _CAPD_DYNSET_C1SET_H_

#include <string>

#include "capd/diffAlgebra/TimeRange.h"
#include "capd/dynset/EnclosureHolder.h"
#include "capd/dynsys/C1DynSys.h"
#include "capd/dynset/SetTraits.h"

namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/**
 * Common interface of all sets that store C1 informations (set position and first derivatives)
 */
template<typename MatrixT>
class C1Set : public capd::diffAlgebra::TimeRange<typename MatrixT::ScalarType>,
              public capd::dynset::C1EnclosureHolder<MatrixT>{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef C1Set SetType;                                     ///< defines class of given set (C0, C1, C2, ...)
  typedef capd::dynsys::C1DynSys<MatrixType> DynSysType;     ///< defines minimal regularity of the dynamical system

  /// constructor, initializes default enclosures and initial time
  C1Set(const VectorType& x, const VectorType& enc, const MatrixType& M, const MatrixType& mEnc, const ScalarType& t)
      : capd::diffAlgebra::TimeRange<ScalarType>(t),
        capd::dynset::C1EnclosureHolder<MatrixType>(x,enc,M,mEnc)
  {}

  /// destructor
  virtual ~C1Set(){}

  /// computes image of the set after one step/iterate of the dynamical system
  virtual void move(DynSysType & c1dynsys) = 0;
  const static size_type degree() { return 1; }
};

template<class MatrixT>
struct SetTraits< C1Set<MatrixT> >{
	const static bool isC0Set=false;
	const static bool isC1Set=true;
	const static bool isC2Set=false;
	const static bool isCnSet=false;
};

/// @}
}} // the end of the namespace capd::dynset

#endif // _CAPD_DYNSET_C1SET_H_

