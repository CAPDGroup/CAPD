/////////////////////////////////////////////////////////////////////////////
/// @file C0FlowballSet.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#ifndef _CAPD_DYNSET_C0FLOWBALLSET_H_
#define _CAPD_DYNSET_C0FLOWBALLSET_H_

#include <stdexcept>
#include "capd/dynset/C0Set.h"
#include "capd/vectalg/Norm.h"

namespace capd{
namespace dynset{

template<typename MatrixT>
class C0FlowballSet : public C0Set<MatrixT>{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;
  typedef typename MatrixType::size_type size_type;

  C0FlowballSet(const VectorType& x, const ScalarType& r, const NormType& aNorm, ScalarType t = TypeTraits<ScalarType>::zero());
  C0FlowballSet(const C0FlowballSet&);
  ~C0FlowballSet(){delete m_n;}

  virtual void move(capd::dynsys::DynSys<MatrixType>& dynsys);
  virtual void move(capd::dynsys::DynSys<MatrixType>& dynsys, C0FlowballSet& result) const;
  virtual std::string show(void) const;
  virtual VectorType affineTransformation(const MatrixType&, const VectorType&) const;

protected:
  VectorType m_x;
  ScalarType m_r;
  NormType *m_n;
};


template<typename MatrixType>
inline C0FlowballSet<MatrixType>::C0FlowballSet(const C0FlowballSet& S)
  : C0Set<MatrixType>(S),
    m_x(S.m_x), m_r(S.m_r), m_n(S.m_n->clone())
{}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C0FLOWBALLSET_H_
