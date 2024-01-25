

/////////////////////////////////////////////////////////////////////////////
/// @file FadMap.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_FAD_MAP_H_
#define _CAPD_DYNSYS_FAD_MAP_H_

#include "capd/fadbad/fadiff.h"
#include "capd/fadbad/tadiff.h"
#include "capd/fadbad/differentiate.h"
#include "capd/vectalg/Vector.h"
#include "capd/vectalg/Matrix.h"
#include "capd/diffAlgebra/Hessian.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

// this class is usually used as a Poincare section
// an implementation should provide operator() and gradient methods

template <class VectorT>
class FadFunction
{
public:
  typedef VectorT VectorType;
  typedef typename VectorT::ScalarType ScalarType;

  virtual ScalarType operator()(const VectorType& /*u*/) const =0;
  virtual VectorType gradient(const VectorType& /*u*/) const =0;    
  virtual ~FadFunction() {}
};

// an interface for parameter class for BasicFadTaylor
template<typename Scalar, unsigned D>
class FadMap
{
public:
  typedef Scalar ScalarType;
  typedef capd::vectalg::Matrix<ScalarType,D,D> MatrixType;
  typedef capd::vectalg::Vector<ScalarType,D> VectorType;
  typedef FadFunction<VectorType> FunctionType;
  typedef capd::diffAlgebra::Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;

  virtual ~FadMap() {}

  // Any inherited class must specify the following operator(). 
  // It should compute value of a nonautonomous vector field on a vector u.
  // If the ODE is autonomous, simply do not use the 't' argument in the definition.
  template <class TimeT, typename AVector>
  AVector operator()(const TimeT& /*t*/, const AVector& /*u*/) const;
  
  /// For backward compatibility.
  /// It computes derivative of the vector field for an autonomous vector field.
  /// You can implement it as follow.
  template <typename AVector>
  AVector operator()(const AVector& /*u*/) const;

  /// This operator should compute derivative of the vector field.
  /// One can specify own implementation
  /// or use in the implementation a template function
  /// computeDerivative from file differentiate.h (as in the example of the Lorenz system below)
  /// which performs FAD to compute derivative of a map.
  virtual MatrixType derivative(const ScalarType& /*t*/, const VectorType& /*u*/) const = 0;

  /// computes simultaneously value and derivative of the map for a given vector and time
  virtual VectorType operator()(ScalarType /*t*/, const VectorType&, MatrixType &) const = 0;

  /// For backward compatibility.
  /// It computes derivative of the vector field for an autonomous vector field.
  virtual MatrixType operator[](const VectorType& /*u*/) const = 0;

  /// This function should set parameter value of the vector field
  /// as a Description we may use integers or strings, etc.
  template <typename Description>
  void setParameter(Description /*s*/, Scalar /*value*/){}
  
  /// You must implement this method. It must return the dimension of the phase space.
  unsigned dimension() {return static_cast<unsigned>(D);}

  /// Maximal order of spacial derivatives that this map can compute
  unsigned degree() { return 1u;}
};


/// Sample implementation of FadMap.
/// This class implements the vector field for the Lorenz system.
/// Template parameters are:
/// Scalar: double, interval, MpFloat, MpInterval, etc.
/// D: in this case either 3 or 0. If D=3 then vectors and matrices are allocated
/// on stack and the computations are much faster but we must know they dimension at compile time.
/// This forces separate compilation of all the classes like vectors and matrices for this particular dimension.
/// D=0 means that vectors and matrices are allocated on the storage and they can be of arbitrary
/// dimension specified at runtime.
template<class Scalar, unsigned D>
class LorenzFadMap : public capd::dynsys::FadMap<Scalar,D>
{
public:
  typedef typename capd::dynsys::FadMap<Scalar,D>::HessianType HessianType;
  typedef typename capd::dynsys::FadMap<Scalar,D>::MatrixType MatrixType;
  typedef typename capd::dynsys::FadMap<Scalar,D>::VectorType VectorType;
  typedef FadFunction<VectorType> FunctionType;
  
  LorenzFadMap() : parameters(3)
  {}

  LorenzFadMap(Scalar _s, Scalar _r, Scalar _q)
    : parameters(3)
  {
    parameters[0]=_s;
    parameters[1]=_r;
    parameters[2]=_q;
  }
  /// This operator computes vector field of the Lorenz system.
  template <typename AVector>
  AVector operator()(const AVector& in) const
  {
    AVector out(3);
    out[0] = parameters[0]*(in[1]-in[0]);
    out[1] = in[0]*(parameters[1]-in[2])-in[1];
    out[2] = in[0]*in[1]-parameters[2]*in[2];
    return out;
  }

  /// This operator computes vector field of the Lorenz system.
  /// The system is autonomous, but the operator must have two arguments in order to fit requested signature.
  template <typename TimeT, typename AVector>
  AVector operator()(const TimeT& /*t*/, const AVector& in) const
  {
    return (*this)(in);
  }
  
  MatrixType derivative(const Scalar& t, const VectorType& u) const
  {
    MatrixType der(3,3);
    computeDerivative(*this,t,u,der);
    return der;
  }

  /// computes simultaneously value and derivative of the map for a given vector and time
  VectorType operator()(Scalar t, const VectorType& x, MatrixType& o_der) const
  {
    return computeDerivative(*this,t,x,o_der);
  }

  MatrixType operator[](const VectorType& u) const
  {
    MatrixType der(3,3);
    computeDerivative(*this,TypeTraits<Scalar>::zero(),u,der);
    return der;
  }
  
  void setParameter(int i, Scalar value)
  {
    parameters[i] = value;
  }
  
  unsigned dimension() const {return 3u;}

  std::vector<Scalar> parameters;
  // s = parameters[0]
  // r = parameters[1]
  // q = parameters[2]
};


// This simple class defines a Poincare section for the Lorenz system for use with FadTaylor
template<class Scalar, unsigned D>
class LorenzSection
{
public:
  typedef Scalar ScalarType;
  typedef capd::vectalg::Vector<ScalarType,D> VectorType;

  LorenzSection() : grad(3), parameters(3)
  {
    grad[2] = ScalarType(1);
  }

  LorenzSection(Scalar _s, Scalar _r, Scalar _q)
    : grad(3), parameters(3)
  {
    grad[2] = ScalarType(1);
    parameters[0]=_s;
    parameters[1]=_r;
    parameters[2]=_q;
  }

  // zero of this operator defines the Poincare section
  ScalarType operator()(const VectorType& u) const
  {
    return u[2] - parameters[1] + 1;
  }

  // section must provide a gradient method
  VectorType gradient(const VectorType& /*u*/) const
  {
    return grad;
  }

  VectorType grad;


  void setParameter(unsigned i, Scalar value)
  {
    parameters[i] = value;
  }

  unsigned getDimension() const {return 3;}

  // s = parameters[0]
  // r = parameters[1]
  // q = parameters[2]
  std::vector<Scalar> parameters;
};
/// @}
}} // the end of the namespace capd::dynsys

#endif //define _CAPD_DYNSYS_FAD_MAP_H_


