

/////////////////////////////////////////////////////////////////////////////
/// @file FirstOrderEnclosure.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef CAPD_DYNSYS_FIRST_ORDER_ENCLOSURE_HPP
#define CAPD_DYNSYS_FIRST_ORDER_ENCLOSURE_HPP

#include <sstream>
#include <string>
#include <stdexcept>

#include "capd/vectalg/iobject.hpp"
#include "capd/dynsys/FirstOrderEnclosure.h"
#include "capd/dynsys/SolverException.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

// the function finds an enclosure for \varphi([0,step],x)
template<typename MapType>
typename MapType::VectorType FirstOrderEnclosure::enclosure(  MapType  & vField,
                                         typename MapType::ScalarType const & currentTime,
                                         typename MapType::MatrixType::RowVectorType const & x,
                                         typename MapType::ScalarType const & step
                                      ) {
  typedef typename MapType::ScalarType ScalarType;
  typedef typename MapType::VectorType VectorType;
  typedef typename TypeTraits<ScalarType>::Real Real;

  ScalarType trialStep = ScalarType(Real(-0.2),Real(1.2))*step;
  int dimension = x.dimension();
  VectorType y(dimension),z(dimension);

  ScalarType h = ScalarType(Real(0.),Real(1.)) * step;
  ScalarType timeRange = currentTime + h;
  typename ScalarType::BoundType multf = 1.5;   // factor to multiply coordinates if inclusion fails

  VectorType val = vField(timeRange,x);
  z = x + trialStep * val + ScalarType(-Real(1.),Real(1.)) * TypeTraits<ScalarType>::epsilon();

  bool found = false;
  int counter=0,
      limit=10 + 2*x.dimension(),    // maximum numbers of attempts to find enclosure
      i;
  while((!found) && (counter<limit)){
    counter++;
    y = x + h * vField(timeRange,z);
    found = true;
    for(i=0;i< dimension;++i){
      if(!(y[i].subsetInterior(z[i]))){
        found = false;
        z[i] = y[i];
        ScalarType s;
        z[i].split(s);
        s = multf*s;
        z[i] += s;
      }
    }
  }

  if(found) return y;
  throw SolverException<VectorType>("Solver error: cannot find enclosure guaranteeing bounds",currentTime,x,step);
}


//###########################################################//

template<typename MapType, typename NormType>
typename MapType::MatrixType FirstOrderEnclosure::jacEnclosure(
                        const MapType& vectorField,
                        const typename MapType::ScalarType& currentTime,
                        const typename MapType::ScalarType& step,
                        const typename MapType::VectorType& enc,
                        const NormType &the_norm,
                        typename MapType::ScalarType* o_logNormOfDerivative
                        )
// the function finds enclosure for Jacobian matrix (variational part)
// source- "C^1-Lohner algorithm" by P. Zgliczynski
{
  typedef typename MapType::MatrixType MatrixType;
  typedef typename MapType::ScalarType ScalarType;
  const static ScalarType I(TypeTraits<ScalarType>::zero().leftBound(),TypeTraits<ScalarType>::one().rightBound());

  int dimension = enc.dimension();
  ScalarType h = I*step;
  MatrixType der = vectorField.derivative(currentTime+h,enc), W(dimension,dimension); // W_3 in paper "C^1 - Lohner algorithm"

  ScalarType l = the_norm(der).rightBound(); // computation of lagarithmic norm
  W = ScalarType(-1,1)*exp(h*l);

  MatrixType result = MatrixType::Identity(dimension) + h*der*W;
  capd::vectalg::intersection(W,result,result);
  if(o_logNormOfDerivative)
    *o_logNormOfDerivative = l;
  return result;
}

//###########################################################//

template<typename MapType>
typename MapType::ScalarType FirstOrderEnclosure::c2Enclosure(
      const MapType& vectorField,
      const typename MapType::ScalarType& step,
      const typename MapType::VectorType& enc,
      typename MapType::MatrixType& jacEnclosure,
      typename MapType::HessianType& hessEnclosure
    )
{
  typedef typename MapType::ScalarType ScalarType;
  typedef typename MapType::VectorType VectorType;
  typedef typename MapType::MatrixType MatrixType;
  typedef typename MapType::HessianType HessianType;

  const static ScalarType I(TypeTraits<ScalarType>::zero().leftBound(),TypeTraits<ScalarType>::one().rightBound());
  ScalarType h = I*step;
  
  int dimension = enc.dimension();
  vectorField.homogenousPolynomial(jacEnclosure);
  ScalarType logNormOfDerivative = capd::vectalg::EuclLNorm<VectorType,MatrixType>()(jacEnclosure).rightBound(); // computation of lagarithmic norm
  MatrixType W(dimension,dimension); // W_3 in paper "C^1 - Lohner algorithm"

  ScalarType hl = exp(h*logNormOfDerivative);
  W = ScalarType(-1,1)*hl;
  jacEnclosure = MatrixType::Identity(dimension) + h*jacEnclosure*W;
  capd::vectalg::intersection(W,jacEnclosure,jacEnclosure);

  int i,j,c;
  HessianType temp(dimension);
  vectorField.homogenousPolynomial(jacEnclosure,temp);

  ScalarType w = 
    logNormOfDerivative.contains(0.0) ? 
        abs(step).rightBound() : (hl-ScalarType(1.))/logNormOfDerivative;
  for(int j=0;j<dimension;++j)
    for(int c=j;c<dimension;++c)
    {
      ScalarType delta =TypeTraits<ScalarType>::zero();
      for(i=0;i<dimension;++i)
        delta += sqr(temp(i,j,c));
      delta = sqrt(nonnegativePart(delta)).rightBound();
      typename ScalarType::BoundType size = (abs(delta * w)).rightBound();
      for(i=0;i<dimension;++i)
        hessEnclosure(i,j,c) = ScalarType(-size,size);
    } // c - loop

  return logNormOfDerivative;
}


/// @}
}} //namespace capd::dynsys

#endif // CAPD_DYNSYS_FIRST_ORDER_ENCLOSURE_HPP


