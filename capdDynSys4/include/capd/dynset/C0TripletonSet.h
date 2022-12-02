/////////////////////////////////////////////////////////////////////////////
/// @file C0TripletonSet.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0TRIPLETONSET_H_
#define _CAPD_DYNSET_C0TRIPLETONSET_H_

#include "capd/dynset/C0Set.h"
#include "capd/geomset/CenteredTripletonSet.h"

namespace capd{
namespace dynset{
/// @addtogroup dynset 
/// @{
/**
 * Class \b C0TripletonSet \b represents a subset of R^n in the following form
 *
 * x + C*r0 + intersection(B*r , Q*q )
 *
 * where
 *  x is a point vector
 *  C,B,Q are point matrices, where B and Q are invertible and Q is close to orthogonal
 *  r0,r are interval vectors
 *
 *  Moreover it stores rigorous inverse matrices of B and Q
 *
 */

template<typename MatrixT, typename Policies>
class C0TripletonSet;

template<typename MatrixT, typename Policies>
C0TripletonSet<MatrixT,Policies> operator*(const MatrixT&A, const C0TripletonSet<MatrixT,Policies>& s);

template<typename MatrixT, typename Policies>
class C0TripletonSet : public Policies, public C0Set<MatrixT>, public capd::geomset::CenteredTripletonSet<MatrixT>, public TripletonData<MatrixT>
{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixT::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef C0Set<MatrixT> SetType;
  typedef typename C0Set<MatrixT>::DynSysType DynSysType;
  typedef capd::geomset::CenteredTripletonSet<MatrixT> BaseSet;
  typedef TripletonData<MatrixT> Data;
  typedef Policies Policy;

  C0TripletonSet(const VectorType& x, ScalarType t = TypeTraits<ScalarType>::zero());
  C0TripletonSet(const VectorType& x, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C0TripletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C0TripletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  C0TripletonSet(const VectorType& x, const MatrixType& C,
      const VectorType& r0, const MatrixType& B,
      const VectorType& r,
      ScalarType t = TypeTraits<ScalarType>::zero()
    );
  friend C0TripletonSet operator*<>(const MatrixT&, const C0TripletonSet&);

  /// computes image of the set after one step/iterate of the dynamical system
  void move(DynSysType & dynsys);

  /// this computes next representation of the set given computed one-step enclosure of the form y + jacPhi*deltaX + rem
  static void move(const BaseSet& set, BaseSet& result, VectorType& bound, Data& data);

  /// computes image of the set after one step/iterate of the dynamical system and stores it in result
  void move(DynSysType & dynsys, C0TripletonSet& result) const;

  void setToIdentity(){
    this->m_B.setToIdentity();
    this->m_invB.setToIdentity();
    this->m_r.clear();
    this->m_Q.setToIdentity();
    this->m_invQ.setToIdentity();
    this->m_q.clear();
  }

  /// This method computes value of functor f at interval vector represented by this set.
  /// This set is represented as tripleton X=x+C*r0+intersection(B*r,Q*q). Then f(X) can be computed as
  /// f(x) + (Df(X)*C)*r0 + intersection( (Df(X)*B)*r, (Df(X)*Q)*q )
  template<class Functional>
  ScalarType evalAt(const Functional& f) const {
    ScalarType r = f(this->m_currentSet);
    VectorType gradient = f.gradient(this->m_currentSet);
    if(subset(this->m_x,this->m_currentSet)){
      intersection( r, f(this->m_x) + this->evalAffineFunctional(gradient,this->m_x), r);
    }else{
      VectorType x0 = midVector(this->m_currentSet);
      intersection( r, f(x0) + this->evalAffineFunctional(gradient,x0), r);
    }

    return r;
  }

  using SetType::operator VectorType;
  using SetType::m_currentSet;

  std::string show() const;
  std::string name() const { return "C0TripletonSet"; }

  using BaseSet::get_r0;
  using BaseSet::getElement_r0;
  using BaseSet::get_x;
  using BaseSet::getElement_x;
  using BaseSet::get_B;
  using BaseSet::get_invB;
  using BaseSet::getElement_B;
  using BaseSet::getRow_B;
  using BaseSet::getColumn_B;
  using BaseSet::get_C;
  using BaseSet::getElement_C;
  using BaseSet::getRow_C;
  using BaseSet::getColumn_C;

  // computed rigorous inverses of m_B and m_Q
protected:
  C0TripletonSet(const VectorType& x,
      const MatrixType& C, const VectorType& r0,
      const MatrixType& B, const VectorType& r,
      const MatrixType& Q, const VectorType& q,
      ScalarType t = TypeTraits<ScalarType>::zero()
    );
};

template<typename MatrixT, typename Policies>
C0TripletonSet<MatrixT,Policies> operator*(const MatrixT&A, const C0TripletonSet<MatrixT,Policies>& s)
{
  return C0TripletonSet<MatrixT,Policies>(A*s.m_x,A*s.m_C,s.m_r0,A*s.m_B,s.m_r,A*s.m_Q,s.m_q);
}

/// @}
}} //namespace capd::dynset

#endif // _CAPD_DYNSET_C0TRIPLETONSET_H_
