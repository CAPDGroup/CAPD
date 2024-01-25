/////////////////////////////////////////////////////////////////////////////
/// @file  PoincareMap_templateMembers.h
///
/// @author Daniel Wilczak
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 200-2013 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_POINCARE_MAP_TEMPLATE_MEMBERS_H_
#define _CAPD_POINCARE_POINCARE_MAP_TEMPLATE_MEMBERS_H_

#include <cassert>
#include "capd/poincare/PoincareMap.h"
#include "capd/poincare/BasicPoincareMap.hpp"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
// -------------------- protected functions ---------------------------------

/*__________________________________________________________________________*/
/// Function returns section sign. It also checks If we crossed section
/// and returned in one step. In this case it throws an exception.
/// @param[in]      theSet      the set after making a step,
/// @return    section sign evaluated on theSet
template <typename SolverT, typename FunT>
template<typename T>
typename PoincareMap<SolverT,FunT>::ScalarType
PoincareMap<SolverT,FunT>::getSign(const T & theSet)
{
  if(this->m_section->isSpecialSection())
    return this->m_section->evalAt(theSet);
  else
    return theSet.evalAt(*this->m_section);
}

/*__________________________________________________________________________*/
/// Function returns section sign. It also checks If we crossed section
/// and returned in one step. In this case it throws an exception.
/// @param[in]      theSet      the set after making a step,
/// @param[in, out] position    position before time step,
/// @param[in]      updatePosition   flag that decides if position has to be updated
/// @param[in]      bound       rough enclosure for trajectory during whole time step
/// @return    section sign
template <typename SolverT, typename FunT>
template<typename T>
typename PoincareMap<SolverT,FunT>::ScalarType
PoincareMap<SolverT,FunT>::getSign(const T & theSet,  VectorType & position, bool updatePosition, const VectorType & bound )
{
  // we check if we cross the section and than return during one step
  // i.e. if during section crossing the section gradient is orthogonal to vector field
  this->checkTransversability(theSet, bound);
  // Now we are sure that either sing is constant during the time step
  // or we crossed the section transversely, so we can compute sign on position after the step.
  if(updatePosition)
    position = (VectorType)theSet;

  return this->getSign(theSet);
}

/*__________________________________________________________________________*/
/// Function returns section sign. It also checks for possible crossing the section
/// and returning in one step. In this case it throws an exception.
/// @param[in]      theSet      the set after making a step,
/// @param[in, out] position    position before time step,
/// @param[in]      updatePosition   flag that decided if position has to be updated
/// @return    section sign
template <typename SolverT, typename FunT>
template<typename T>
typename PoincareMap<SolverT,FunT>::ScalarType
PoincareMap<SolverT,FunT>::getSign( const T & theSet,  VectorType & position, bool updatePosition)
{
  return this->getSign(theSet, position, updatePosition, theSet.getLastEnclosure());
}

/*__________________________________________________________________________*/

template <typename SolverT, typename FunT>
template<typename T>
void PoincareMap<SolverT,FunT>::integrateUntilSectionCrossing(T& theSet, T& nextSet, int n)
// Functions reaches the section. As result:
// - theSet is just before the section,
// - nextSet is on the section or just after it

{
  this->m_solver.getStepControl().init(this->m_solver,theSet.getCurrentTime(),theSet);
  T *prev = &theSet, *next = &nextSet;

  for(int iteration=0;iteration<n;++iteration){
     VectorType setPosition = (VectorType)(*prev);
     this->m_signBeforeSection = getSign(*prev);

    //------------------------------------------------------------------------------------------------
    // LEAVING SECTION AND GOING UP TO THE POINT WHERE SECTION HAS CORRECT SIGN
    // We want to leave section and reach point where the sign of the section function,
    // according to crossingDirection, indicate that next section crossing will be in a good direction.

    while( (m_signBeforeSection.contains(0.0)) ||                  // We are on the section
        !((this->m_crossingDirection == Both) ||         // section sign is not correct
            ((this->m_crossingDirection == MinusPlus) && (m_signBeforeSection < 0.0)) ||
            ((this->m_crossingDirection == PlusMinus) && (m_signBeforeSection > 0.0))
        )
    ){
      this->m_solver(*prev); // we make one step to leave the section and try again
      m_signBeforeSection = this->getSign(*prev, setPosition, true);
    }
    //------------------------------------------------
    // first return to the section

    assert(!m_signBeforeSection.contains(ScalarType(0.0)));

    m_signBeforeSection = m_signBeforeSection.mid();            // we save the sign of the section function, sign does not contain zero
    ScalarType sign = m_signBeforeSection;      // used for checking if we cross the section (if sign is changed)

    // now we try to move the set closer to section than one current time step
    while( (sign>0. and m_signBeforeSection>0.) or (sign<0. and m_signBeforeSection<0.) ) // m_signBeforeSection*sign>0
    {
      this->m_solver(*prev, *next);
      sign = this->getSign(*next, setPosition, true);
      std::swap(prev, next);
    } // end while
  }

  if(&theSet == prev)
    std::swap(theSet,nextSet);
  this->m_lastTimeStep = this->m_solver.getStep();
}

///*__________________________________________________________________________*/
//
template <typename SolverT, typename FunT>
template<typename T>
void PoincareMap<SolverT,FunT>::getCloserToSection(T& theSet, T& nextSet){
  SaveStepControl<SolverT> ssc(this->m_solver);
  const CurveType& c = this->m_solver.getCurve();

  // We know that at t=c.getLeftDomain() the set is before section
  // and at t = c.getRightDomain() it is after the section or on the section
  // Using Interval Newton method we extract from this domain subinterval,
  // when the set cannot be on the section. This will give us a lower bound for the section crossing time.
  ScalarType returnTime = ScalarType(c.getLeftDomain(),c.getRightDomain());

  ScalarType timeBeforeSection = theSet.getCurrentTime();

  // now we resolve for the return time using interval Newton method
  VectorType bound = nextSet.getLastEnclosure();
  ScalarType sign = m_signBeforeSection;

  int i=0;
  for(;i<10;++i)
  {
    VectorType vfOnBound = this->m_solver.getVectorField()(timeBeforeSection+returnTime,bound);
    ScalarType h = returnTime.mid();
    ScalarType valAtH = (*this->m_section)(c(h));
    ScalarType derAtDomain = this->m_section->gradientByVector(bound,vfOnBound);
    ScalarType N = h - valAtH/derAtDomain;

    // There are three cases:
    // 1. N \subset stepMade. This means the set crossed the section in one step
    // 2. N and stepMade are disjoint. This meand that the set did not touch the section. This should not happen.
    // 3. N and stepMade overlap without inclusion. Then possible return time cannot be smaller than N.leftBound().
    // Summarizing, min(N.leftBound(),stepMade.rightBound()) is a lower bound for return time.

    ScalarType newReturnTime;
    if ( !intersection(N,returnTime,newReturnTime) ) // case 2. Take the bound at the end of the step.
    {
      theSet = nextSet;
      return;
    }

    if(returnTime == newReturnTime) // cannot improve
      break;
    returnTime = newReturnTime;
    bound = c(returnTime);
  }
  if(i){
    this->m_solver.setMaxStep(returnTime.leftBound());
    this->m_solver(theSet);
  }
}

/////////////////////////////////////////////////////////////////////////////////////
///
/// Function used to cross the section for all types of sets
///
/// @param[in] theSet                  set just before the section
/// @param[out] oneStepReturnTime      bound for the return time in relative to current time of theSet
/// @param[out] bound                  bound for the intersection of the trajectories with Poincare section
/// @returns        true               if succeed
template <typename SolverT, typename FunT>
template<typename T>
bool PoincareMap<SolverT,FunT>::crossSectionInOneStep(T& theSet, T& nextSet, ScalarType& oneStepReturnTime, VectorType& bound)
{
  SaveStepControl<SolverT> ssc(this->m_solver);
  if(m_signBeforeSection.contains(0.0))
    throw PoincareException<T>("PoincareMap::crossSectionInOneStep exception: initial set is already on the section ", theSet, m_signBeforeSection );
  this->getCloserToSection(theSet,nextSet);
   // make one step
  this->m_solver.setStep(this->m_lastTimeStep);
  this->m_solver(theSet,nextSet);
  ScalarType sign = this->getSign(nextSet);

  if( !(m_signBeforeSection*sign < 0) ) // we did not cross the section in one step
    return false;

  bound = nextSet.getLastEnclosure();
  this->checkTransversability(theSet, bound);

  // now we resolve for the return time using interval Newton method
  const CurveType& c = this->m_solver.getCurve();
  oneStepReturnTime = ScalarType(c.getLeftDomain(),c.getRightDomain());
  ScalarType timeBeforeSection = theSet.getCurrentTime();
  VectorType vfOnBound = this->m_solver.getVectorField()(timeBeforeSection+oneStepReturnTime,bound);

  VectorType x0 = c.getCenter();
  for(int i=0;i<10;++i)
  {
    ScalarType h = oneStepReturnTime.mid();
    // We want to estimate m_section(c(h))
    // taking into account representation of the set.
    // The code
    // ScalarType valAtH = (*this->m_section)(c(h));
    // leads to huge overestimation on computation of valAtH and, in consequence, the return time.
    // Below, we estimate the function as follows.
    // Assume Phi is a numerical method, so that
    // phi_h(X)\subset Phi_h(x0) + DPhi_h(X)*(X-x0) + rem
    // We compute
    // section(phi_h(X)) \subet section(Phi_h(x_0)) + DPhi_h(X)^T*grad(section)*(X-x0) + grad(section)*rem
    // The last term is evaluated using evalAffineFunctional,
    // that for a given set X and vectors g,x computes g*(X-x) using representation of X.
    VectorType valAtCenter = c.valueAtCenter(h);
    VectorType rem = c.remainder(h);
    ScalarType valAtH = (*this->m_section)(valAtCenter) + this->m_section->gradientByVector(valAtCenter + ScalarType(0.,1.)*rem,rem);
    valAtH += theSet.evalAffineFunctional(Transpose(c.oneStepDerivativeOfNumericalMethod(h))*this->m_section->gradient(bound),x0);
    // intersect with previous implementation
    intersection(valAtH,(*this->m_section)(c(h)),valAtH);
    ScalarType derAtDomain = this->m_section->gradientByVector(bound,vfOnBound);
    ScalarType N = h - valAtH/derAtDomain;
    ScalarType newReturnTime;
    if(!intersection(oneStepReturnTime,N,newReturnTime)){
      std::ostringstream out;
      out.precision(17);
      out << "PoincareMap::crossSectionInOneStep error - empty intersection in estimation of the return time. Report this error to CAPD developers!";
      out << "\nreturnTime1=" << oneStepReturnTime;
      out << "\nreturnTime2=" << N;
      out << "\nvalueAtCentre=" << valAtH;
      out << "\nderivative=" << derAtDomain;
      out << "\nvalueAtLeft=" << (*this->m_section)(c(oneStepReturnTime.leftBound()));
      out << "\nvalueAtRight=" << (*this->m_section)(c(oneStepReturnTime.rightBound()));
      out << "\ntimeStepRequested=" << this->m_lastTimeStep;
      out << "\ntimeStepMade=" << this->m_solver.getStep() << std::endl;

      throw std::runtime_error(out.str());
    }
    bound = c(newReturnTime);
    vfOnBound = this->m_solver.getVectorField()(timeBeforeSection+oneStepReturnTime,bound);
    VectorType L = c(oneStepReturnTime.left());
    VectorType R = c(oneStepReturnTime.right());
    for(size_type j=0;j<bound.dimension();++j){
      if(!(vfOnBound[j].contains(0.))){
        intersection(bound[j],intervalHull(L[j],R[j]),bound[j]);
      }
    }

    if(oneStepReturnTime == newReturnTime)
      break;

    oneStepReturnTime = newReturnTime;
  }

  return true;
}

/////////////////////////////////////////////////////////////////////////////////////
///
/// Function used to cross the section for all types of sets
///
/// @param[in,out] theSet              set just before the section, on exit is equal to setAfterTheSection
/// @param[in]     nextSet  set on the section or just after it
/// @returns        set that contains value of Poincare map for given set theSet
template <typename SolverT, typename FunT>
template<typename T>
typename PoincareMap<SolverT,FunT>::VectorType
PoincareMap<SolverT,FunT>::crossSection(T& theSet, T& nextSet)
{
  SaveStepControl<SolverT> ssc(this->m_solver);

  T *prev = &theSet, *next=&nextSet;
  VectorType result = VectorType(theSet);
  this->sectionDerivativesEnclosure.saveEnclosure(theSet); // we save the \frac{d}{dx}\phi

  while(true){
    const VectorType& oneStepBound = next->getLastEnclosure();
    checkTransversability(*prev, oneStepBound);

    ScalarType sign = this->getSign(*next);
    ScalarType check = sign*m_signBeforeSection;
    if(check<0.) break; // the set is on the opisite side if the section

    const VectorType nextPosition = (VectorType)(*next);
    if( check>0. ){ // we did not touch the section yet
      result = nextPosition;
      this->sectionDerivativesEnclosure.saveEnclosure(*next);
    }else{
      this->sectionDerivativesEnclosure.updateEnclosure(*prev,*next);
      ScalarType domain(0.,this->m_solver.getStep().rightBound());
      VectorType vfOnBound = this->m_solver.getVectorField()(prev->getCurrentTime()+domain,oneStepBound);
      for(size_type i=0;i<this->getVectorField().dimension();i++)
      {
        if(!(vfOnBound[i].contains(0.)))
          result[i] = intervalHull(result[i],nextPosition[i]);
        else
          result[i] = intervalHull(result[i],oneStepBound[i]); // bound is an enclosure
      }
    }

    std::swap(prev,next);
    this->m_solver(*prev,*next);
  }

  // Now *next is after the section.
  // We have to resolve for the upper bound on return time.
  const CurveType& c = this->m_solver.getCurve();
  ScalarType stepMade = ScalarType(c.getLeftDomain(),c.getRightDomain());
  ScalarType t0 = prev->getCurrentTime();
  VectorType bound = next->getLastEnclosure();
  VectorType vfOnBound = this->m_solver.getVectorField()(t0+stepMade,bound);

  for(int i=0;i<10;++i)
  {
    ScalarType h = stepMade.mid();
    ScalarType valAtH = (*this->m_section)(c(h));
    ScalarType derAtDomain = this->m_section->gradientByVector(bound,vfOnBound);
    ScalarType N = h - valAtH/derAtDomain;

    ScalarType newReturnTime;
    if ( !intersection(N,stepMade,newReturnTime) ) // case 2. Take the bound at the end of the step.
    {
      std::ostringstream out;
      out.precision(17);
      out << "PoincareMap::crossSection error - empty intersection in estimation of the return time. Report this error to CAPD developers!";
      out << "\nreturnTime1=" << stepMade;
      out << "\nreturnTime2=" << N << std::endl;
      throw std::runtime_error(out.str());
    }

    bound = c(stepMade);
    vfOnBound = this->m_solver.getVectorField()(t0+stepMade,bound);
    VectorType L = c(stepMade.left());
    VectorType R = c(stepMade.right());
    for(size_type j=0;j<bound.dimension();++j){
      if(!(vfOnBound[j].contains(0.))){
        intersection(bound[j],intervalHull(L[j],R[j]),bound[j]);
      }
    }

    if(stepMade == newReturnTime) // cannot improve
      break;
    stepMade = newReturnTime;
  }

  // Eventually, we add bound from the last step. Newton method guarantees that
  // at t=stepMade.rightBound() the set is after section, even if the obtained enclosure touches the section!
  this->m_solver.setStep(stepMade.rightBound());
  this->m_solver(*prev,*next);

  const VectorType& oneStepBound = next->getLastEnclosure();
  checkTransversability(*prev, oneStepBound);
  this->sectionDerivativesEnclosure.updateEnclosure(*prev,*next);
  stepMade = ScalarType(0.,this->m_solver.getStep().rightBound());
  vfOnBound = this->m_solver.getVectorField()(prev->getCurrentTime()+stepMade,oneStepBound);
  const VectorType nextPosition = (VectorType)(*next);

  for(size_type i=0;i<this->getVectorField().dimension();i++)
  {
    if(!(vfOnBound[i].contains(0.)))
      result[i] = intervalHull(result[i],nextPosition[i]);
    else
      result[i] = intervalHull(result[i],oneStepBound[i]); // bound is an enclosure
  }

  if(prev==&theSet)
    std::swap(theSet,nextSet);
  return result;
}

/// Function checks if we crossed section and then returned in one step.
/// In this case it throws an exception of PoincareException<T> type.
template <typename SolverT, typename FunT>
template<typename T>
void PoincareMap<SolverT,FunT>::checkTransversability(
  const T & theSet, const VectorType & bound
){
  ScalarType check = (*this->m_section)(bound);
  const typename SolverT::SolutionCurve& c = this->m_solver.getCurve();
  ScalarType domain(c.getLeftDomain(),c.getRightDomain());

  if(subset(ScalarType(0.0), check)) {  // Is the section crossed?
    ScalarType innerProduct = this->m_section->gradientByVector(bound,this->m_solver.getVectorField()(domain,bound));
    if (innerProduct.contains(0.0)) {  // Is the vector field orthogonal to section gradient?
      throw PoincareException<T>(
              "PoincareMap error: possible nontransversal return to the section ", theSet, theSet,
              bound, this->m_solver.getVectorField()(domain,bound),
              check, this->m_section->gradient(bound), innerProduct
              );
    }
  }

}

/// @}
}} // namespace capd::poincare

#endif // _CAPD_POINCARE_POINCARE_MAP_TEMPLATE_MEMBERS_H_

