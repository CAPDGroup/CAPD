/////////////////////////////////////////////////////////////////////////////
/// @file C1DoubletonSetGeometricTail.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2020 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_PDES_C1DOUBLETONSETGEOMETRICTAIL_H_
#define _CAPD_PDES_C1DOUBLETONSETGEOMETRICTAIL_H_

#include "capd/dynset/C1DoubletonSet.hpp"
#include "capd/pdes/GeometricBound.h"
#include "capd/pdes/PdeSolver.h"
#include "capd/pdes/C0DoubletonSetGeometricTail.h"

namespace capd{
namespace pdes{
/// @addtogroup pdes 
/// @{

template<typename Policies>
class C1DoubletonSetGeometricTail : public C0DoubletonSetGeometricTail< capd::dynset::C1DoubletonSet<GeometricBound<capd::interval>::MatrixType,Policies> >{
public:
  typedef capd::pdes::GeometricBound<capd::interval> VectorType;
  typedef VectorType::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType FiniteVectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef capd::vectalg::MaxNorm<FiniteVectorType,MatrixType> NormType;
  typedef capd::dynset::C1DoubletonSet<MatrixType,Policies> FiniteDimensionalBaseSet;
  typedef C0DoubletonSetGeometricTail< capd::dynset::C1DoubletonSet<MatrixType,Policies> > BaseSet;
  typedef typename capd::dynset::C1Set<MatrixType>::SetType SetType;
  typedef capd::pdes::PdeSolver<VectorType> DynSysType;
  typedef capd::dynset::DoubletonData<MatrixType> Data;
  typedef Policies Policy;
  typedef std::vector<VectorType> VectorArray;
  typedef typename DynSysType::MatrixArray MatrixArray;

  C1DoubletonSetGeometricTail(const VectorType& x, ScalarType t = TypeTraits<ScalarType>::zero())
    : BaseSet(x,t), m_dyxId(x.dimension(),VectorType(x.dimension())), m_dyx(x.dimension(),VectorType(x.dimension()))
  {
    initC1();
  }

  C1DoubletonSetGeometricTail(const VectorType& x, const FiniteVectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero())
    : BaseSet(x,r0,t), m_dyxId(r0.dimension(),VectorType(x.dimension())), m_dyx(r0.dimension(),VectorType(x.dimension()))
  {
    initC1();
  }

  C1DoubletonSetGeometricTail(const VectorType& x, const MatrixType& C, const FiniteVectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero())
    : BaseSet(x,C,r0,t), m_dyxId(r0.dimension(),VectorType(x.dimension())), m_dyx(r0.dimension(),VectorType(x.dimension()))
  {
    initC1();
  }

  C1DoubletonSetGeometricTail(const VectorType& x, const MatrixType& C, const FiniteVectorType& r0, const MatrixType& B, const FiniteVectorType& r, ScalarType t = TypeTraits<ScalarType>::zero())
    : BaseSet(x,C,r0,B,r,t), m_dyxId(r0.dimension(),VectorType(x.dimension())), m_dyx(r0.dimension(),VectorType(x.dimension()))
  {
    initC1();
  }

  C1DoubletonSetGeometricTail(const VectorType& x, const FiniteDimensionalBaseSet& set)
    : BaseSet(x,set), m_dyxId(x.dimension(),VectorType(x.dimension())), m_dyx(x.dimension(),VectorType(x.dimension()))
  {
    initC1();
  }

  /// computes image of the set after one step/iterate of the dynamical system
  void move(DynSysType & dynsys) { move(dynsys,*this); }
  /// computes image of the set after one step/iterate of the dynamical system and stores it in result
  void move(DynSysType & dynsys, C1DoubletonSetGeometricTail& result) const;

  std::string name() const { return "C1DoubletonSetGeometricTail"; }

  void initC1(){
    for(size_type j=0;j<m_x.dimension();++j)
      m_dyxId[j][j] = 1;
    this->encDyx = this->Dyx = 0.;
    this->Dy = FiniteVectorType(m_x.dimension()+1);
    this->Dy[m_x.dimension()] = 1;
    this->encDy =this->Dy;
    this->m_lastMatrixEnclosure = (MatrixType)(*this);
    this->encDxx = this->m_lastMatrixEnclosure;
    counter = 0;
  }

  inline const VectorType& ind(VectorArray& J, int j) const { return J[j]; }
  inline const VectorType& ind(MatrixArray& J, int j) const { return J[j][0]; }

  /// Here we compute operator norm of f:R^n\to l_\infty
  /// In R^n we assume l_2 norm.
  /// Hence, operator norm is given by
  /// supremum of l_2 norms of rows.
  /// This can be also bouned by supremum of l_1 norms of rows.
  template<class Container>
  ScalarType getVyxNorm(Container& J) const{
    const size_type m = m_dyx.size();
    const size_type M = this->m_currentSeries.dimension();
    ScalarType vyx = 0.;
    /// first check finite number of rows, where the coefficients are given explicitly.
    for(size_type k=m+1;k<=M;++k){
      ScalarType s = 0;
      for(size_type j=0;j<m;++j)
        s += sqr(ind(J,j).getCoefficient(k));
      vyx = intervalHull(sqrt(s),vyx);
    }
    /// check the first row which is bounded uniformly as a geometric series.
    ScalarType s = 0;
    ScalarType q = power(ind(J,0).getGeometricDecay(),-(M+1));
    for(size_type c=0;c<m;++c){
      s += sqr(q*ind(J,c).getConstant());
    }
    return intervalHull(vyx,sqrt(s));
  }

  void setMaxBlockSize(int newMaxBlockSize) { this->maxBlockSize = newMaxBlockSize; }

  using BaseSet::move;
  using BaseSet::affineTransformation;
  using BaseSet::operator VectorType;
  using FiniteDimensionalBaseSet::operator FiniteVectorType;
  using BaseSet::getCurrentSeries;
  using BaseSet::setLastEnclosure;
  using BaseSet::getLastEnclosure;

  using BaseSet::get_x;
  using BaseSet::getElement_x;
//protected:
  using BaseSet::m_x;
  using BaseSet::m_r;
  using BaseSet::m_enc;
  using BaseSet::m_currentSeries;
  using BaseSet::m_tmp;

  VectorArray m_dyxId, m_dyx;

  ScalarType Dxx, Dyx, encDyx;
  FiniteVectorType Dy, encDy;
  MatrixType encDxx;

  bool printLog = false;
  int maxBlockSize = 500;
  int counter;
};

// -----------------------------------------------------------------------------------

template<typename Policies>
void C1DoubletonSetGeometricTail<Policies>::move(DynSysType & dynsys, C1DoubletonSetGeometricTail& res) const
{
  C1DoubletonSetGeometricTail result = *this;
  BaseSet::initMove(result);
  dynsys.encloseC1Map(
    result.x,m_currentSeries,Dy,
    result.m_tmp, result.y, result.rem, result.m_enc,
    result.jacPhi, result.jacRem, result.m_lastMatrixEnclosure,
    result.m_dyxId, result.m_dyx,
    result.Dy, result.encDy
  );

  // move C^0 part
  capd::dynset::C0DoubletonSet<MatrixType,Policy>::move(*this,result,result.m_currentSet,result);

  // move C^1 part
  // first enclosure of the image is given by simply interval evaluation
  result.jacPhi += result.jacRem;
  result.m_currentMatrix = result.jacPhi*this->m_currentMatrix;
  FiniteDimensionalBaseSet::move(*this,result,result.m_currentMatrix,result);

  result.setCurrentTime(this->getCurrentTime()+dynsys.getStep());

  FiniteDimensionalBaseSet::reorganizeIfNeeded(result.m_B,result.m_invB,result.m_r,result.m_C,result.m_r0);
  FiniteDimensionalBaseSet::reorganizeC1IfNeeded(result.m_Bjac,result.m_invBjac,result.m_R,result.m_Cjac,result.m_R0);

  BaseSet::finalizeMove(result);

  // update Dxx,Dyx (finite-width) column-block that is
  // newDxx = Dxx(Id,0)*oldDxx + Dxx(0,oldDyx)
  // newDyx = Dyx(Id,0)*oldDxx + Dyx(0,oldDyx)
  const size_type m = m_dyx.size();
  const size_type M = m_currentSeries.dimension();
  // Loop over columns
  for(size_type i=0;i<m;++i){
    // Firstly, we compute newDxx.
    // Member result.m_currentMatrix already contains Dxx(Id,0)*oldDxx,
    // hence we need to add Dxx(0,oldDyx)
    // Leading coefficients of Dxx(0,oldDyx) will be absorbed to m_R, which is part of newDxx doubleton representation
    auto tmp = result.m_dyx[i].projection(m);
    result.m_R.column(i) += result.m_invBjac*tmp;
    result.m_currentMatrix.column(i) += tmp;

    // Secondly, compute new Dyx.
    // We already have result.m_dyx that cotaints Dyx(0,oldDyx)
    // We update it by Dyx(Id,0)*oldDxx, which splits into two parts
    // - finite number of explicit coefficients in (finite) rows
    // - uniform bound on infinite number of rows
    // Here we process finite number of rows
    // Note: indexing in dyx and dyxId is [column][row]
    for(size_type j=m;j<M;++j){
      for(size_type k=0;k<m;++k)
        result.m_dyx[i][j] += result.m_dyxId[k][j]*this->m_currentMatrix[k][i];
    }
    // Eventualy, we must compute a uniform constant for infinite-number of rows.
    // Here C = constant for Dyx(0,oldDyx)
    ScalarType C = result.m_dyx[i].getConstant();
    // Add to it constant bounding the product Dyx(Id,0)*oldDxx
    for(size_type k=0;k<m;++k)
      C += result.m_dyxId[k].getConstant()*abs(this->m_currentMatrix[k][i]).rightBound();
    result.m_dyx[i].setConstant(C);
  }
  for(size_type i=0;i<m;++i)
    result.m_dyxId[i].setConstant(0.);

  IMaxNorm norm;
  result.Dxx = norm(result.m_currentMatrix);
  result.Dyx = result.getVyxNorm(result.m_dyx);

  // eventualy, compute enclosure of Dyx over time step (needed for Poincare map)
  // newDyx = Dyx(Id,0)*oldDxx + Dyx(0,oldDyx). We sum two norms
  result.encDyx = getVyxNorm(dynsys.getMatrixRemainderCoefficients())*this->Dxx + getVyxNorm(dynsys.getDyxRemainderCoefficients());

  result.encDxx = result.getLastMatrixEnclosure()*this->m_currentMatrix;
  for(size_type j=1;j<=result.encDxx.numberOfColumns();++j){
    for(size_type i=1;i<=result.encDxx.numberOfRows();++i){
      result.encDxx(i,j) += dynsys.getDyxRemainderCoefficients()[j-1][0][i-1];
    }
  }

  if(result.printLog and result.counter++ > result.maxBlockSize)
  {
    result.counter = 0;
    std::ostringstream out;
    out.precision(17);
    out << "time=" << result.m_currentTime << ", " << result.m_currentSet[0] << std::endl;
    out << "maxWidth(Dxx)=" << maxWidth(result.m_currentMatrix) << std::endl;
    out << "Dy=" << result.Dy  << std::endl;
    out << "Dyx=" << result.Dyx  << std::endl << std::endl;
    std::cout << out.str();
  }

  res = result;
}

/// @}
}} // namespace capd::dynset

#endif // _CAPD_PDES_C1DOUBLETONSETGEOMETRICTAIL_H_
