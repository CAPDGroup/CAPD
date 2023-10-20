/////////////////////////////////////////////////////////////////////////////
/// @file C0DoubletonSetGeometricTail.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2020 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_PDES_C0DOUBLETONSETGEOMETRICTAIL_H_
#define _CAPD_PDES_C0DOUBLETONSETGEOMETRICTAIL_H_

#include "capd/dynset/C0DoubletonSet.hpp"
#include "capd/dynset/C0TripletonSet.hpp"
#include "capd/dynset/lib.h"
#include "capd/pdes/GeometricBound.hpp"
#include "capd/pdes/PdeSolver.h"

namespace capd{
namespace pdes{
/// @addtogroup pdes 
/// @{

template<typename BaseT = capd::dynset::C0DoubletonSet<GeometricBound<capd::interval>::MatrixType,capd::C0Rect2Policies> >
class C0DoubletonSetGeometricTail : public BaseT{
public:
  typedef capd::pdes::GeometricBound<capd::interval> VectorType;
  typedef VectorType::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType FiniteVectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef BaseT FiniteDimensionalBaseSet;
  typedef typename capd::dynset::C0Set<MatrixType>::SetType SetType;
  typedef capd::pdes::PdeSolver<VectorType> DynSysType;
  typedef capd::dynset::DoubletonData<MatrixType> Data;
  typedef typename FiniteDimensionalBaseSet::Policy Policy;

  C0DoubletonSetGeometricTail(const VectorType& x, ScalarType t = TypeTraits<ScalarType>::zero())
    : FiniteDimensionalBaseSet(FiniteVectorType(x.dimension(),x.getExplicitCoefficients().begin()),t), m_currentSeries(x), m_tmp(x)
  {}

  C0DoubletonSetGeometricTail(const VectorType& x, const FiniteVectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero())
    : FiniteDimensionalBaseSet(x.projection(r0.dimension()),r0,t), m_currentSeries(x), m_tmp(x)
  {
    for(size_type i=0;i<r0.dimension();++i)
      m_currentSeries[i] += r0[i];
  }

  C0DoubletonSetGeometricTail(const VectorType& x, const MatrixType& C, const FiniteVectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero())
    : FiniteDimensionalBaseSet(x.projection(r0.dimension()),C,r0,t), m_currentSeries(x), m_tmp(x)
  {
    m_currentSeries.projection(r0.dimension()) += C*r0;
  }

  C0DoubletonSetGeometricTail(const VectorType& x, const MatrixType& C, const FiniteVectorType& r0, const MatrixType& B, const FiniteVectorType& r, ScalarType t = TypeTraits<ScalarType>::zero())
    : FiniteDimensionalBaseSet(x.projection(r0.dimension()),C,r0,B,r,t), m_currentSeries(x), m_tmp(x)
  {
    m_currentSeries.projection(r0.dimension()) += C*r0 + B*r;
  }

  C0DoubletonSetGeometricTail(const VectorType& x, const FiniteDimensionalBaseSet& set)
    : FiniteDimensionalBaseSet(set), m_currentSeries(x), m_tmp(x)
  {}

  void initMove(C0DoubletonSetGeometricTail& result) const;
  void finalizeMove(C0DoubletonSetGeometricTail& result) const;

  /// computes image of the set after one step/iterate of the dynamical system
  void move(DynSysType & dynsys) { move(dynsys,*this); }

  /// computes image of the set after one step/iterate of the dynamical system and stores it in result
  void move(DynSysType & dynsys, C0DoubletonSetGeometricTail& result) const;

  /// This method computes value of functor f at interval vector represented by this set.
  template<class Functional>
  ScalarType evalAt(const Functional& f) const {
    return FiniteDimensionalBaseSet::evalAt(f.getProjection(FiniteDimensionalBaseSet::dimension()));
  }

  VectorType affineTransformation(const MatrixType& A, const VectorType& c) const{
    size_type d = FiniteDimensionalBaseSet::dimension();
    MatrixType M(d,d);
    for(size_type i=1;i<=d;++i)
      for(size_type j=1;j<=d;++j)
    M(i,j) = A(i,j);
    FiniteVectorType u =  FiniteDimensionalBaseSet::affineTransformation(M,FiniteVectorType(d,c.getExplicitCoefficients().begin()));
    VectorType result = m_currentSeries - c;
    for(size_type i=0;i<d;++i)
      result[i] = u[i];
    return result;
  }

  ScalarType evalAffineFunctional(const VectorType& gradient, const VectorType& x0) const{
    size_type d = FiniteDimensionalBaseSet::dimension();
    ScalarType r = FiniteDimensionalBaseSet::evalAffineFunctional(
      FiniteVectorType(d,gradient.getExplicitCoefficients().begin()),
      FiniteVectorType(d,x0.getExplicitCoefficients().begin())
    );    
    for(size_type i=d;i<x0.dimension();++i)
      r += gradient[i]*(this->m_currentSeries[i]-x0[i]);
    return r;
  }
  
  operator VectorType() const { return m_currentSeries; }
  using FiniteDimensionalBaseSet::operator FiniteVectorType;

  std::string name() const { return "C0DoubletonSetGeometricTail"; }

  const VectorType& getCurrentSeries() const { return m_currentSeries; }
  void setLastEnclosure(const VectorType& enc) { m_enc = enc; }
  const VectorType& getLastEnclosure() const { return m_enc; }

  using FiniteDimensionalBaseSet::get_x;
  using FiniteDimensionalBaseSet::getElement_x;
  using FiniteDimensionalBaseSet::affineTransformation;
//protected:
  using FiniteDimensionalBaseSet::m_x;
  using FiniteDimensionalBaseSet::m_r;
  VectorType m_enc, m_currentSeries, m_tmp;
};

// -----------------------------------------------------------------------------------

template<typename BaseT>
void C0DoubletonSetGeometricTail<BaseT>::initMove(C0DoubletonSetGeometricTail& result) const
{
  // the following function can throw an exception leaving output parameters in an inconsistent state
  // do not overwrite parameters of the set until we are sure that they are computed correctly
  if(subset(this->m_x,this->m_currentSet)){
      result.x = this->m_x;
      result.deltaX = this->m_currentSet-this->m_x;
  }else{
    split(this->m_currentSet, result.x, result.deltaX);
  }
  for(size_type i=0;i<result.m_currentSet.dimension();++i){
    ScalarType t;
    if(!intersection(this->m_currentSet[i],this->m_currentSeries[i],t)){
      std::cout
          << "C0DoubletonSetGeometricTail - inconsistent data in set (intersection error): "
          << i << '\n' << result.m_currentSeries[i] << '\n' << result.m_currentSet[i] << std::endl;
      exit(0);
    }
  }
}

// -----------------------------------------------------------------------------------

template<typename BaseT>
void C0DoubletonSetGeometricTail<BaseT>::finalizeMove(C0DoubletonSetGeometricTail& result) const
{
  FiniteDimensionalBaseSet::reorganizeIfNeeded(result.m_B,result.m_invB,result.m_r,result.m_C,result.m_r0);

  swap(result.m_currentSeries,result.m_tmp);
  for(size_type i=0;i<result.m_currentSet.dimension();++i){
    if(!intersection(result.m_currentSeries[i],result.m_currentSet[i],result.m_currentSet[i])){
      std::cout
          << "C0DoubletonSetGeometricTail intersection error: "
          << i << '\n' << result.m_currentSeries[i] << '\n' << result.m_currentSet[i] << std::endl;
      exit(0);
    }
    result.m_currentSeries[i]=result.m_currentSet[i];
  }
}
// -----------------------------------------------------------------------------------

template<typename BaseT>
inline void C0DoubletonSetGeometricTail<BaseT>::move(DynSysType & dynsys, C0DoubletonSetGeometricTail& res) const
{
  C0DoubletonSetGeometricTail result = *this;
  this->initMove(result);

  dynsys.encloseC0Map(result.x,m_currentSeries,result.m_tmp,result.y, result.rem, result.m_enc, result.jacPhi);
  capd::dynset::C0DoubletonSet<MatrixType,Policy>::move(*this,result,result.m_currentSet,result);
  result.setCurrentTime(this->getCurrentTime()+dynsys.getStep());

  this->finalizeMove(result);
  res = result;
}

/// @}
}} // namespace capd::dynset

#endif // _CAPD_PDES_C0DOUBLETONSETGEOMETRICTAIL_H_
