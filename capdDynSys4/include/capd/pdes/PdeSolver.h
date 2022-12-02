/// @addtogroup pdes
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file PdeSolver.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008-2016 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_PDES_PDESOLVER_H_
#define _CAPD_PDES_PDESOLVER_H_

#include <vector>
#include <iostream>
#include <sstream>
#include <tuple>

#include "capd/basicalg/minmax.h"
#include "capd/basicalg/power.h"
#include "capd/intervals/lib.h"
#include "capd/vectalg/lib.h"
#include "capd/dynsys/approveRemainder.h"
#include "capd/pdes/PdeCurve.h"
#include "capd/pdes/DissipativeVectorField.h"

namespace capd {
namespace pdes {

//template<class SeriesT, class StepControlT = capd::dynsys::FixedStepControl<capd::interval> >
template<class SeriesT, class StepControlT = capd::dynsys::ILastTermsStepControl>
class PdeSolver : public PdeCurve<SeriesT>, public capd::dynsys::StepControlInterface<StepControlT,typename SeriesT::ScalarType> {
public:
  typedef SeriesT VectorType;
  typedef StepControlT StepControlType;

  typedef capd::IVector FiniteVectorType;
  typedef PdeCurve<VectorType> CurveType;
  typedef PdeCurve<VectorType> SolutionCurve;
  typedef capd::IMatrix MatrixType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef typename MatrixType::size_type size_type;
  //typedef DissipativeVectorField<VectorType> MapType;
  typedef DissipativeVectorField<VectorType> VectorFieldType;
  typedef typename CurveType::VectorArray VectorArray;
  typedef typename CurveType::MatrixArray MatrixArray;

  PdeSolver(VectorFieldType& f, size_type order);
  ~PdeSolver(){}

  const CurveType& getCurve()  const {
    this->setDomain(0.,rightBound(this->m_step));
    return *this;
  }
  CurveType& getCurve() {
    this->setDomain(0.,rightBound(this->m_step));
    return *this;
  }
  const CurveType& getImplicitCurve() const { return this->m_implicitCurve; }

  template <typename SetType>
  void operator()(SetType& set){
    this->saveCurrentSet(set);
    set.move(*this);
  }

  /// Computes image of the set (in set's representation) and stores it in the result set.
   /// @param[in]  set       C^0 or C^1  set representing initial conditions
   /// @param[out] result    on return contains image of the set
  template <typename SetType>
  void operator()(SetType& set, SetType& result){
    this->saveCurrentSet(set);
    set.move(*this,result);
  }

  VectorFieldType& getVectorField() { return m_vectorField; }

  void encloseC0Map(
      const FiniteVectorType& x0, //< @param[in] an internal point of x, usually center of x
      const VectorType& x,        //< @param[in] set to be moved along by the flow
      VectorType& y,      	      //< @param[out] result - set after one step
      FiniteVectorType& o_phi,    //< @param[out] bound for phi(x0), where phi is a numerical method
      FiniteVectorType& o_rem,    //< @param[out] bound for the error of numerical method over the time step
      VectorType& o_enc,    	    //< @param[out] enclosure of all trajectories starting from x over the time interval (time step of numerical method)
      MatrixType& o_jacPhi  	    //< @param[out] bound for derivative Dphi(x) on finite number of modes
  );

  void encloseC1Map(
      const FiniteVectorType& x0,    //< @param[in] an internal point of x, usually center of x
      const VectorType& x,           //< @param[in] set to be moved along by the flow
      const FiniteVectorType& inDy,  //< @param[in] vector of norms of blocks Dxy: m\times \infinity and Dyy: \infinity\times\infinity
      VectorType& y,      	         //< @param[out] result - set after one step
      FiniteVectorType& o_phi,       //< @param[out] bound for phi(x0), where phi is a numerical method
      FiniteVectorType& o_rem,       //< @param[out] bound for the error of numerical method over the time step
      VectorType& o_enc,    	       //< @param[out] enclosure of all trajectories starting from x over the time interval (time step of numerical method)
      MatrixType& o_jacPhi,  	       //< @param[out] bound for derivative Dphi(x) on finite number of modes
      MatrixType& o_jacRem, 	       //< @param[out] bound for the error of numerical method over a time step for Dphi(x) on finite number of modes
      MatrixType& o_jacEnc, 	       //< @param[out] enclosure for Dxx block over last step
      VectorArray& DyxId,            //< @param[out] solution to variational equation on finite number of leading columns with initial condition (Id,0) )
      VectorArray& Dyx,              //< @param[in,out] solution to variational equation on finite number of leading columns with initial condition (0,Vyx)
      FiniteVectorType& outDy,       //< @param[out] bound on one step norms of blocks Dxy: m\times \infinity and Dyy: \infinity\times\infinity - after time step
      FiniteVectorType& encDy        //< @param[out] bound on one step norms of blocks Dxy: m\times \infinity and Dyy: \infinity\times\infinity - after time step
   );

  void computeImplicitCoefficients(
    const FiniteVectorType& x0,     //< @param[in] an internal point of x, usually center of x
    const VectorType& x,        	  //< @param[in] set to be moved along the trajectories of ODE
    size_type q
  );

  void initRemainderCoefficients(ScalarType /*t*/, const VectorType& /*x*/, unsigned /*degree*/){
    this->m_vectorField.computeODECoefficients(this->getRemainderCoefficients(),this->getOrder()+1);
  }

  void computeTimeStep(const ScalarType& t, const VectorType& x){
    this->m_step = this->isStepChangeAllowed()
        ? this->getStepControl().computeNextTimeStep(*this,t,x.getExplicitCoefficients())
        : capd::min(this->m_fixedTimeStep,this->getMaxStep());
  }

  void setStep(ScalarType h) {
      m_fixedTimeStep = h;
      this->turnOffStepControl();
  }
  void adjustTimeStep(ScalarType h) { this->m_step = h; }
  ScalarType getStep() { return m_step; }

  ScalarType getCoeffNorm(size_type, size_type degree) const;
protected:
  typedef std::tuple<VectorArray&,MatrixArray&,MatrixArray&> C1Data;

  /// this function is used to check, that predicted enclosure is indeed an enclosure for leading modes
  /// used both in C^0 and C^1 computation
  template<class V>
  bool checkRemainderInclusion(const V& T, const FiniteVectorType& rem, ScalarType& beta);

  // functions specific for C^0 algorithm
  void makeSelfConsistentBound(VectorArray& enc)                                                            { this->m_vectorField.makeSelfConsistentBound(enc); }
  void computeODECoefficients(VectorArray& a,size_type p)                                                   { this->m_vectorField.computeODECoefficients(a,p); }
  void scaleRemByStepPower(const VectorArray& enc, ScalarType s)                                            { m_rem = enc[this->getOrder()+1].getExplicitCoefficients()*s; }
  bool checkRemainderInclusion(const VectorArray&, ScalarType& beta)                                        { return this->checkRemainderInclusion(m_rem,predicted_rem,beta); }
  void predictEnclosure(const VectorArray& X, VectorArray& enc, ScalarType R, ScalarType I, ScalarType IP)  { this->predictEnclosure(X,enc,predicted_rem,R,I,IP); }

  // functions specific for C^1 algorithm
  void makeSelfConsistentBound(C1Data enc) {
    this->m_vectorField.makeSelfConsistentBound(std::get<0>(enc),std::get<1>(enc),std::get<2>(enc),this->m_numberOfExplicitColumns);
  }

  void computeODECoefficients(C1Data a,size_type p) {
    this->m_vectorField.computeODECoefficients(std::get<0>(a),std::get<1>(a),p,this->m_numberOfExplicitColumns);
    for(size_type i=0;i<this->m_numberOfExplicitColumns;++i)
      this->m_vectorField.computeODECoefficients(std::get<0>(a),std::get<2>(a)[i],p);
  }

  void scaleRemByStepPower(C1Data enc, ScalarType s) {
    m_rem = std::get<0>(enc)[this->getOrder()+1].getExplicitCoefficients()*s;
    for(size_type i=0;i<this->m_numberOfExplicitColumns;++i){
      m_jacRem.column(i) = std::get<1>(enc)[i][this->getOrder()+1].getExplicitCoefficients()*s;
      m_dyxRem.column(i) = std::get<2>(enc)[i][this->getOrder()+1].getExplicitCoefficients()*s;
    }
  }

  bool checkRemainderInclusion(C1Data, ScalarType& beta) {
    bool success = checkRemainderInclusion(m_rem,predicted_rem,beta);
    for(size_type i=0;i<this->m_numberOfExplicitColumns;++i){
      success = checkRemainderInclusion(m_jacRem.column(i),predicted_jacRem[i],beta) and success;
      success = checkRemainderInclusion(m_dyxRem.column(i),predicted_dyxRem[i],beta) and success;
    }
    return success;
  }

  void predictEnclosure(C1Data X, C1Data enc, ScalarType R, ScalarType I, ScalarType IP) {
    predictEnclosure(std::get<0>(X),std::get<0>(enc),predicted_rem,R,I,IP);
    for(size_type i=0;i<this->m_numberOfExplicitColumns;++i){
      predictEnclosure(std::get<1>(X)[i],std::get<1>(enc)[i],predicted_jacRem[i],R,I,IP);
      predictEnclosure(std::get<2>(X)[i],std::get<2>(enc)[i],predicted_dyxRem[i],R,I,IP);
    }
  }

  /// this function is used to predict enclosure for C^0 part and for finite number of column in C^1 algorithm
  void predictEnclosure(const VectorArray& X, VectorArray& enc, FiniteVectorType& rem, ScalarType R, ScalarType I, ScalarType IP);

  /// this function realizes initial computation usd in both for C^0 and C^1 algorithms
  void initEncloseOneStep(const FiniteVectorType& x0, const VectorType& x);

  /// computes output from one step of both C^0 and C^1 algorithms, that is o_phi, o_jacPhi and o_rem
  void sumTaylorSeriesOfLeadingModes(size_type m, VectorType& x, FiniteVectorType& o_phi, FiniteVectorType& o_rem, MatrixType& o_jacPhi);


  template<class EncType>
  void highOrderEnclosure(const EncType& X, EncType& enc);
  void highOrderEnclosure();

  void setInitialCondition(const FiniteVectorType& x0, const VectorType& x, CurveType& curve);

  // TimeRange is base for all types of sets and nonrigorous CxCoeff
  void saveCurrentSet(const capd::diffAlgebra::TimeRange<ScalarType>& /*set*/){
  }

  void saveCurrentSet(const capd::dynset::C1Set<MatrixType>& set){
    this->setInitMatrix((MatrixType)set);
  }

  CurveType m_implicitCurve;
  VectorFieldType& m_vectorField;
  ScalarType m_step;
  ScalarType m_fixedTimeStep;

  MatrixType m_jacPhi, m_jacRem, m_dyxPhi, m_dyxRem;
  FiniteVectorType m_phi, m_deltaX, m_rem;

  FiniteVectorType predicted_rem;
  std::vector<FiniteVectorType> predicted_jacRem;
  std::vector<FiniteVectorType> predicted_dyxRem;

  int m_numberOfExplicitColumns;
}; // end of class PdeSolver

// ***************************************************************************
template<class SeriesT,class StepControlT>
inline PdeSolver<SeriesT,StepControlT>::PdeSolver(VectorFieldType& vf, size_type order)
  : CurveType(vf.dimension(), order),
    m_vectorField(vf),
    m_implicitCurve(vf.dimension(),order),
    m_jacPhi(vf.dimension(),vf.dimension()),
    m_jacRem(vf.dimension(),vf.dimension()),
    m_dyxPhi(vf.dimension(),vf.dimension()),
    m_dyxRem(vf.dimension(),vf.dimension()),
    m_phi(vf.dimension()),
    m_deltaX(vf.dimension()),
    m_rem(vf.dimension()),
    predicted_rem(vf.dimension()),
    predicted_jacRem(vf.dimension()+1,FiniteVectorType(vf.dimension())),
    predicted_dyxRem(vf.dimension()+1,FiniteVectorType(vf.dimension()))
{
  if (order<1)
    throw std::runtime_error("PdeSolver constructor - order of the method cannot be less than 1.");
}

// ***************************************************************************

template<class SeriesT,class StepControlT>
void PdeSolver<SeriesT,StepControlT>::setInitialCondition(
      const FiniteVectorType& x0, const VectorType& x, CurveType& curve
    )
{
  size_type i;
  const size_type m = x0.dimension();
  const size_type M = this->m_vectorField.dimension();
  if(x.dimension()!=M)
    throw std::runtime_error("PdeSolver::setInitialCondition - unequal working dimensions of solver and initial condition.");

  VectorArray& center = curve.getCoefficientsAtCenter();
  VectorArray& X = curve.getCoefficients();
  MatrixArray& J = curve.getMatrixCoefficients();

  center[0] = x;
  for(i=0;i<m;++i){
    center[0][i] = x0[i];
    m_deltaX[i] = x[i]-x0[i];
  }
  this->m_numberOfExplicitColumns = m;
  for(;i<m_numberOfExplicitColumns;++i)
    center[0][i].split(m_deltaX[i]);

  // set initial condition
  X[0] = x;
  // set Identity as an initial condition for variational equations
  for(size_type j=0;j<m_numberOfExplicitColumns;++j){
    for(size_type c=1;c<=M;++c)
      J[j][0].setCoefficient(c, j+1==c ? 1. : 0.);
    J[j][0].setGeometricDecay(x.getGeometricDecay());
    J[j][0].setConstant(0.);
  }
}

// ***************************************************************************

template<class SeriesT,class StepControlT>
void PdeSolver<SeriesT,StepControlT>::initEncloseOneStep(const FiniteVectorType& x0, const VectorType& x){
  this->setInitialCondition(x0,x,*this);
  this->m_vectorField.computeODECoefficients(this->getCoefficientsAtCenter(),this->getOrder());
  this->m_vectorField.computeODECoefficients(this->getCoefficients(),this->getMatrixCoefficients(),this->getOrder(),this->m_numberOfExplicitColumns);
  this->computeTimeStep(0.,x);
}

// ***************************************************************************

template<class SeriesT,class StepControlT>
void PdeSolver<SeriesT,StepControlT>::sumTaylorSeriesOfLeadingModes(
    size_type m, VectorType& x,
    FiniteVectorType& o_phi, FiniteVectorType& o_rem, MatrixType& o_jacPhi
){
  size_type M = this->m_vectorField.dimension();
  size_type i,j;

  this->sumTaylorSeries(m_phi,this->getCoefficientsAtCenter(),m_step,this->getOrder(),this->m_numberOfExplicitColumns);
  this->sumTaylorSeries(m_jacPhi,this->getMatrixCoefficients(),m_step,this->getOrder(),this->m_numberOfExplicitColumns);
  this->sumTaylorSeries(x,this->getCoefficients(),m_step,this->getOrder(),M);
  x.getExplicitCoefficients() += m_rem;
  for(i=0;i<m;++i){
    o_rem[i] = m_rem[i];
    o_phi[i] = m_phi[i];
    for(j=m;j<this->m_numberOfExplicitColumns;++j)
      o_phi[i] += m_jacPhi(i+1,j+1)*m_deltaX[j];
    ScalarType t = o_phi[i] + o_rem[i];
    for(j=0;j<m;++j){
      o_jacPhi(i+1,j+1) = m_jacPhi(i+1,j+1);
      t += m_jacPhi(i+1,j+1)*m_deltaX[j];
    }
    if(!intersection(t,x[i],x[i])){
      std::ostringstream out;
      out << "Intersection error in sumTaylorSeriesOfLeadingModes. Report this to CAPD developers.\n";
      out << t << std::endl;
      out << x[i] << std::endl;
      out << i << std::endl;
      out << m_deltaX << std::endl;
      throw std::logic_error(out.str());
    }
  }
}

// ***************************************************************************

template<class SeriesT,class StepControlT>
void PdeSolver<SeriesT,StepControlT>::encloseC0Map(
    const FiniteVectorType& x0, const VectorType& x, VectorType& y,
    FiniteVectorType& o_phi, FiniteVectorType& o_rem, VectorType& o_enc,
    MatrixType& o_jacPhi
){
  VectorArray& enc = this->getRemainderCoefficients();

  this->initEncloseOneStep(x0,x);
  this->highOrderEnclosure(this->getCoefficients(),enc);
  this->sumTaylorSeriesOfLeadingModes(x0.dimension(),y,o_phi,o_rem,o_jacPhi);

  this->m_vectorField.updateTail(y,enc,m_step);
  o_enc = enc[0];
}


// ***************************************************************************

template<class SeriesT,class StepControlT>
void PdeSolver<SeriesT,StepControlT>::encloseC1Map(
    const FiniteVectorType& x0, const VectorType& x, const FiniteVectorType& inDy,
    VectorType& y, FiniteVectorType& o_phi, FiniteVectorType& o_rem, VectorType& o_enc,
    MatrixType& o_jacPhi, MatrixType& o_jacRem, MatrixType& o_jacEnc,
    VectorArray& DyxId, VectorArray& Dyx,
    FiniteVectorType& outDy, FiniteVectorType& encDy
){
  VectorArray& X = this->getCoefficients();
  VectorArray& enc = this->getRemainderCoefficients();
  MatrixArray& dyx = this->getDyxCoefficients();
  MatrixArray& jacEnc = this->getMatrixRemainderCoefficients();

  // Step 1a. Set initial condition and compute ODE coefficients
  //          for two IVPs: x(0)=x0, C^0-computation
  //          x(0)=x  and V(0)=(Id,0), C^1-computation
  this->initEncloseOneStep(x0,x);
  // Step 1b. Additional IVP for variational part, necessary to reduce exponential growth of the norm of Dxy block of result.
  //          That is, we compute norm of the product instead of product of separate norms.
  //          0 - on leading Vxx  block
  //          data from Dyx in Vyx block
  //          x(0) = x
  //          V(0) = (0,\pi_{yx}Dyx)
  //          coefficients for C^0 part are already computed in previous step
  for(size_type j=0;j<Dyx.size();++j){
    dyx[j][0] = Dyx[j];
    for(size_type c=1;c<=m_numberOfExplicitColumns;++c)
      dyx[j][0].setCoefficient(c,0.);
    this->m_vectorField.computeODECoefficients(X,dyx[j],this->getOrder());
  }

  // Step 2. - predict remainder for all IVPs.
  //         - compute ODE coefficients and find enclosure
  //         - check inclusion and accept time step on success
  //         - shorten the time step if allowed and accept time step
  //         - refine enclosure if changing time step not allowed
  this->highOrderEnclosure();

  // Step 3. Sum Taylor series and add remainder.
  this->sumTaylorSeriesOfLeadingModes(x0.dimension(),y,o_phi,o_rem,o_jacPhi);

  // Step 4. Copy to o_jacRem a finite-dimensional block of solution to variational equation
  for(size_type j=1;j<=o_jacRem.numberOfColumns();++j){
    for(size_type i=1;i<=o_jacRem.numberOfRows();++i){
      o_jacRem(i,j) = m_jacRem(i,j);
      o_jacEnc(i,j) = this->getMatrixRemainderCoefficients()[j-1][0][i-1];
    }
  }

  // Step 5. Add explicit remainder to finite-dimensional block of leading columns and set jacEnc remainder
  for(size_type i=0;i<Dyx.size();++i){
    this->sumTaylorSeries(Dyx[i],dyx[i],this->m_step,this->getOrder(),x.dimension());
    Dyx[i].getExplicitCoefficients() += m_dyxRem.column(i);
    DyxId[i].setExplicitCoefficients(m_jacPhi.column(i)+m_jacRem.column(i));
  }

  // Step 6. Update all tails using linear differential inequality
  this->m_vectorField.updateTail(y,enc,m_step);
  this->m_vectorField.updateTail(DyxId,Dyx,this->getMatrixRemainderCoefficients(),this->getDyxRemainderCoefficients(),m_step);
  o_enc = enc[0];

  // Step 7. compute norm of Dxy and Dyy blocks over one step using component-wise estimations.
  MatrixType J = this->m_vectorField.blockNorms(enc[0],this->m_numberOfExplicitColumns);
  encDy = matrixAlgorithms::matrixExp(J*ScalarType(0.,this->m_step.rightBound()))*inDy;
  outDy = matrixAlgorithms::matrixExp(J*this->m_step)*inDy;
}

// ***************************************************************************

template<class SeriesT,class StepControlT>
void PdeSolver<SeriesT,StepControlT>::computeImplicitCoefficients(
    const FiniteVectorType& x0, const VectorType& x, size_type p
){
  VectorArray& center = this->m_implicitCurve.getCoefficientsAtCenter();
  VectorArray& X = this->m_implicitCurve.getCoefficients();
  MatrixArray& J = this->m_implicitCurve.getMatrixCoefficients();

  this->setInitialCondition(x0,x,this->m_implicitCurve);
  this->m_vectorField.computeODECoefficients(center,p);
  this->m_vectorField.computeODECoefficients(X,J,p);
}

// ***************************************************************************

template<class SeriesT,class StepControlT>
typename PdeSolver<SeriesT,StepControlT>::ScalarType
PdeSolver<SeriesT,StepControlT>::getCoeffNorm(size_type r, size_type /*degree*/) const
{
  typename TypeTraits<ScalarType>::Real result = 0;
  for(size_type i=0;i<this->dimension();++i){
    result = capd::max(result,rightBound(abs(this->remainderCoefficient(i,r))));
    result = capd::max(result,rightBound(abs(this->coefficient(i,r))));
  }
  return ScalarType(result);
}

// ***************************************************************************

template<class SeriesT,class StepControlT>
void PdeSolver<SeriesT,StepControlT>::predictEnclosure(const VectorArray& X, VectorArray& enc, FiniteVectorType& rem, ScalarType R, ScalarType I, ScalarType IP)
{
  size_type d = this->m_vectorField.firstDissipativeIndex()-1;
  size_type p = this->getOrder();

  // on non-disipative modes we predict for high order enclosure
  this->sumTaylorSeries(enc[0].getExplicitCoefficients(),X,I,p,d);
  size_type n=0;
  for(;n<d;++n){
    rem[n] = enc[p+1][n]==0. ? R : ScalarType(-2,2)*IP*enc[p+1][n];
    enc[0][n] += rem[n];
  }
  // on remaining modes we just copy the initial condition
  for(;n<enc[0].getExplicitCoefficients().dimension();++n)
    enc[0][n] = X[0][n];
  // copy constant and decay
  enc[0].setConstant(X[0].getConstant());
  enc[0].setGeometricDecay(X[0].getGeometricDecay());
}

// ***************************************************************************

template<class SeriesT,class StepControlT>
template<class V>
bool PdeSolver<SeriesT,StepControlT>::checkRemainderInclusion(const V& T, const FiniteVectorType& rem, ScalarType& beta){
  size_type d = this->m_vectorField.firstDissipativeIndex()-1;
  size_type p = this->getOrder();
  bool success = true;
  for(size_type n=0;n<d;++n){
    if(subsetInterior(T[n],rem[n])) continue;
    success = false;
    ScalarType c = abs(rem[n]).rightBound()/abs(T[n]).right();
    c = exp(log(c)/(p+1));
    beta = capd::min(beta,c);
  }
  return success;
}


// ***************************************************************************

template<class SeriesT,class StepControlT>
template<class EncType>
void PdeSolver<SeriesT,StepControlT>::highOrderEnclosure(const EncType& X, EncType& enc)
{
  this->m_step = capd::min(this->m_step,this->getMaxStep());
  ScalarType I = ScalarType(0,1)*this->m_step;
  ScalarType IP = power(I,this->getOrder()+1);
  ScalarType R = this->getAbsoluteTolerance()*ScalarType(-1.,1.);

  int counter = 5*this->m_vectorField.firstDissipativeIndex();
  while(counter){
    this->predictEnclosure(X,enc,R,I,IP);
    // refine constant and decay to assert isolation on dissipative modes
    this->makeSelfConsistentBound(enc);
    this->computeODECoefficients(enc,this->getOrder()+1);
    this->scaleRemByStepPower(enc,IP);

    // check inclusion
    ScalarType beta = 1;
    if(checkRemainderInclusion(enc,beta)) return;
    if(this->isStepChangeAllowed() and beta>0){
      if(beta>0.9){
        m_step = (m_step*beta).leftBound();
        return;
      } else {
        this->m_step *= 0.9;
        I = ScalarType(0,1)*this->m_step;
        IP = power(I,this->getOrder()+1);
      }
    }
    counter--;
  }
  throw capd::dynsys::SolverException<VectorType>("PdeSolver::highOrderEnclosure - cannot find an enclosure",0.,this->getCoefficients()[0],m_step);
}

// ***************************************************************************

template<class SeriesT,class StepControlT>
void PdeSolver<SeriesT,StepControlT>::highOrderEnclosure()
{
  C1Data X = make_tuple(std::ref(this->getCoefficients()),std::ref(this->getMatrixCoefficients()),std::ref(this->getDyxCoefficients()));
  C1Data enc = make_tuple(std::ref(this->getRemainderCoefficients()),std::ref(this->getMatrixRemainderCoefficients()),std::ref(this->getDyxRemainderCoefficients()));
  this->highOrderEnclosure(X,enc);
}

}} // namespace capd::pdes


#endif // _CAPD_PDES_PDESOLVER_H_


/// @}
