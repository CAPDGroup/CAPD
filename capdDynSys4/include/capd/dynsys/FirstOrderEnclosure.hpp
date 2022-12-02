

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

#ifndef _CAPD_DYNSYS_FIRST_ORDER_ENCLOSURE_HPP_
#define _CAPD_DYNSYS_FIRST_ORDER_ENCLOSURE_HPP_

#include <sstream>
#include <string>
#include <stdexcept>

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
  MatrixType der = vectorField.derivative(currentTime+h,enc), result(dimension,dimension);

  ScalarType l = the_norm(der).rightBound(); // computation of lagarithmic norm
  ScalarType w = ScalarType(-1,1)*exp(h*l);

  MatrixType W(dimension,dimension); // W_3 in paper "C^1 - Lohner algorithm"
  W = w;

  result = MatrixType::Identity(dimension) + h*der*W;

  int i,j;
  for(i=1;i<=dimension;++i)
    for(j=1;j<=dimension;++j)
    {
      ScalarType d = result(i,j);
      typename ScalarType::BoundType
         l = (w.leftBound() > d.leftBound() ? w.leftBound() : d.leftBound()),
         r = (w.rightBound() < d.rightBound() ? w.rightBound() : d.rightBound());
      result(i,j) = ScalarType(l,r);
    }
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

  int dimension = enc.dimension();
  vectorField.homogenousPolynomial(jacEnclosure);
  ScalarType logNormOfDerivative = capd::vectalg::EuclLNorm<VectorType,MatrixType>()(jacEnclosure).rightBound(); // computation of lagarithmic norm
  ScalarType w = ScalarType(-1,1)*exp(step*logNormOfDerivative);

  MatrixType W(dimension,dimension); // W_3 in paper "C^1 - Lohner algorithm"
  W = w;

  jacEnclosure = MatrixType::Identity(dimension) + step*jacEnclosure*W;

  int i,j,c;
  for(i=1;i<=dimension;++i)
    for(j=1;j<=dimension;++j)
    {
      ScalarType d = jacEnclosure(i,j);
      typename ScalarType::BoundType
         l = (w.leftBound() > d.leftBound() ? w.leftBound() : d.leftBound()),
         r = (w.rightBound() < d.rightBound() ? w.rightBound() : d.rightBound());
      jacEnclosure(i,j) = ScalarType(l,r);
    }

  HessianType temp(dimension);
  vectorField.homogenousPolynomial(jacEnclosure,temp);

  w = (exp(logNormOfDerivative*step)-ScalarType(1.))/logNormOfDerivative;
  for(j=0;j<dimension;++j)
    for(c=j;c<dimension;++c)
    {
      VectorType Bjc(dimension);
      for(i=0;i<dimension;++i)
          Bjc[i] = temp(i,j,c);
      ScalarType delta = Bjc.euclNorm().rightBound();
      typename ScalarType::BoundType size;
      if(logNormOfDerivative.contains(0.0))
        size = abs(delta*abs(step).rightBound()).rightBound();
      else
        size = (abs(delta * w)).rightBound();
      for(i=0;i<dimension;++i)
        hessEnclosure(i,j,c) = ScalarType(-size,size);
    } // c - loop

  return logNormOfDerivative;
}


/// @}
}} //namespace capd::dynsys

#endif // _CAPD_DYNSYS_FIRST_ORDER_ENCLOSURE_HPP_


