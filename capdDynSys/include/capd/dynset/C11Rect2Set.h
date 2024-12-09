/////////////////////////////////////////////////////////////////////////////
/// @file C11Rect2Set.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

/* Author: Daniel Wilczak, 2006 */

#ifndef _CAPD_DYNSET_C11RECT2SET_H_
#define _CAPD_DYNSET_C11RECT2SET_H_

#include "capd/dynset/C1DoubletonSet.h"
#include "capd/dynset/reorganization/FactorReorganization.h"
#include "capd/dynset/QRPolicy.h"

namespace capd{
namespace dynset{
// @addtogroup dynset
/// @{

/**
 * In C1Rect2Set both  C0 and C1 parts are represented as: x + C*r0 + B*r
 *
 * where:
 * -  x - set center
 * -  C*r0 - 'Lipschitz' part
 * -  B*r  - 'errors' computed via QR method
 */
typedef FactorReorganization<FullQRWithPivoting<> > C11Rect2Policies;

template<typename MatrixT>
class C11Rect2Set : public C1DoubletonSet<MatrixT,C11Rect2Policies>
{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef C1DoubletonSet<MatrixT,C11Rect2Policies> BaseSet;
  typedef typename BaseSet::C0BaseSet C0BaseSet;
  typedef typename BaseSet::C1BaseSet C1BaseSet;
  typedef C1Set<MatrixT> SetType;
  typedef typename C1Set<MatrixT>::DynSysType DynSysType;

  C11Rect2Set(const BaseSet& s) : BaseSet(s) {}

  /// computes image of the set after one step/iterate of the dynamical system and stores it in result
  void move(DynSysType & dynsys, C11Rect2Set& result) const;
  using BaseSet::move;

  std::string name() const { return "C11Rect2Set"; }
  std::string toString() const { return C0BaseSet::toString() +"\n" + C1BaseSet::toString(); }
};

/// @}

}} // the end of the namespace capd::dynset

#endif // _CAPD_DYNSET_C11RECT2SET_H_

