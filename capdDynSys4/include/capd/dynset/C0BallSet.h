
/////////////////////////////////////////////////////////////////////////////
/// @file C0BallSet.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0BALLSET_H_
#define _CAPD_DYNSET_C0BALLSET_H_

#include "capd/dynset/C0Set.h"
#include "capd/vectalg/Norm.h"

namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/**
   the set is represented as: x + Ball(r).

   move: x - center of result; r = Lip(JacPhi) + errors from Phi(x)
      for cascades
      for flows : r = exp(Lip(vectfield)*h) + errors ..
 */
template<typename MatrixT>
class C0BallSet : public C0Set<MatrixT>{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;
  typedef typename MatrixType::size_type size_type;

protected:
  VectorType x;
  ScalarType r;
  NormType *n;

public:
  C0BallSet(const C0BallSet&);
  C0BallSet(const VectorType& x, NormType& n, ScalarType t = TypeTraits<ScalarType>::zero());
  C0BallSet(const VectorType& x, const ScalarType& r, NormType& n, ScalarType t = TypeTraits<ScalarType>::zero());
  ~C0BallSet();
  const C0BallSet& operator=(const C0BallSet&);

  virtual void move(capd::dynsys::DynSys<MatrixType>& dynsys);
  virtual void move(capd::dynsys::DynSys<MatrixType>& dynsys, C0BallSet& result) const;
  virtual std::string show(void) const;
  virtual VectorType affineTransformation(const MatrixType&, const VectorType&) const;
};

/// @}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C0BALLSET_H_


