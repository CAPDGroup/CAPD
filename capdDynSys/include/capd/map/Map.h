/// @addtogroup map
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Map.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_MAP_MAP_H_
#define _CAPD_MAP_MAP_H_

#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>
#include "capd/map/Function.h"
#include "capd/diffAlgebra/Hessian.h"
#include "capd/diffAlgebra/Jet.h"

namespace capd{
/// Functions and Maps
namespace map{

/**
 * \brief This class is used to represent a map \f$ f:R^N->R^M \f$.
 *
 * The class provides methods for:
 * - computation of Taylor coefficients of the map at a given vector up to given degree
 * - computing first order derivative and hessian
 * - jet propagation through the map
 *
 * The map is parsed from a human readable string. Syntax of the formula:
 *\code "[par:a,b,c...;][time:t;]var:x1,x2,...;fun:expression1,expression2,....;" \endcode
 * Sections par and time are optional.
 *
 * - \b var - names of subsequent arguments
 * - \b fun - expressions that define a map. You can use most elementary
 *       functions: sin, cos, exp, log, sqrt, sqr
 *       operators: +,-,*,/,^ (power, integer or not. Exponent cannot depend on variables - we assume gradient of exponent is zero).
 *       constants: (0,1,2,-5,2.5,-0.25, etc)- we recommend usage of representable numbers only.
 *                  If a constant is an interval or a floating point represented with high precision, one should use a parameter instead.
 * - \b par - parameters of the map. Derivatives of them with respect to main variables are assumed to be zero.
 * - \b time - a distinguished parameter with derivative dt/dt=1. Used for nonautonomous ODEs
 *
 * Dimensions N and M are automatically recognized from the number of variables and functions.
 * \note If MatrixT is of fixed dimensions (R,C) then these numbers MUST agree (R,C)=(M,N).
 *
 * \b Example:
 * \code
 *    #include "capd/capdlib.h"
 *    std::string lorenzFormula("par:s,r,q;var:x,y,z;fun:s*(y-x),x*(r-z)-y,x*y-q*z;");
 *    capd::IMap lorenz(lorenzFormula);
 *    lorenz.setParameter("s",10.);
 *    lorenz.setParameter("r",28.);
 *    lorenz.setParameter("q",interval(8.)/interval(3.));
 * \endcode
 *
 * The map can be parsed from a C-like routine. One has to write a global function that defines map and pass it to the constructor.
 * The function that defines a map must have the following signature:
 *
 * \code
 * #include "capd/capdlib.h"
 * void vf(capd::autodiff::Node t,
 *         capd::autodiff::Node in[], int dimIn,
 *         capd::autodiff::Node out[], int dimOut,
 *         capd::autodiff::Node params[], int noParam
 *        );
 * \endcode
 *
 * Here:
 * - \b t - is a time, required for non-autonomous ODEs
 * - \b in - is a C-array of input variables
 * - \b dimIn - is number of input variables
 * - \b out - is a C-array of output variables
 * - \b dimOut - is number of output variables
 * - \b params - is a C-array of parameters
 * - \b noParam - is number of parameters
 *
 * \b Example:
 * \code
 * #include "capd/capdlib.h"
 * using capd::autodiff::Node;
 * void lorenzVectorField(Node t, Node in[], int dimIn, Node out[], int dimOut, Node params[], int noParams)
 * {
 *   out[0] = params[0]*(in[1]-in[0]);
 *   out[1] = in[0]*(params[1]-in[2])-in[1];
 *   out[2] = in[1]*in[2]-params[2]*in[2];
 * }
 *
 * int dimIn=3, dimOut=3, noParam=3;
 * capd::IMap lorenz(lorenzVectorField,dimIn,dimOut,noParam);
 * lorenz.setParameter(0,10.);
 * lorenz.setParameter(1,28.);
 * lorenz.setParameter(2,interval(8.)/interval(3.));
 * \endcode
 *
 * This class can be used as a vector field defining an ODE.
 * It provides methods for computation of Taylor coefficients of solutions to ODEs and associated variational equations.
 * Although the user can define own ODE solver based on computed coefficients by class Map, we recommend to use class capd::dynsys::BasicTaylor and inherited from it.
 * These classes provide suitable interface for integration of ODEs and variational equations.
 *
 * <b>About implementation.</b>
 *
 * Directed acyclic graph (DAG) that encodes a map is stored as well indexed array. Thus processing of it is fast.
 * Moreover, this representation is well suited for:
 * - GPU computing (simple array allocation on the device)
 * - source translation techniques
 * During parsing process the expressions are optimized. Thus subexpressions that appear many times in these formulas are computed only once. For example in the vector field:
 * \code
 *    std::string twistedOScillator("var:x,y;fun:y*(x^2+y^2),-x*(x^2+y^2);");
 * \endcode
 * subexpression x^2+y^2 is recognized and computed only once.
 *
 * \remark The parser does not perform factorization of expressions.
 */
template<typename MatrixT>
class Map : public BasicFunction<typename MatrixT::ScalarType>
{
public:
  typedef MatrixT MatrixType;                                     ///< public type: Matrix used in the computations (template parameter)
  typedef typename MatrixType::RowVectorType VectorType;          ///< public type: Vector that represent arguments of the map (dimension is equal to number of variables)
  typedef typename MatrixType::ColumnVectorType ImageVectorType;  ///< public type: Vector that represent value of the map (dimension is equal to number of functions)
  typedef typename MatrixType::RefColumnVectorType RefColumnVectorType;  ///< public type: reference vector that represent value of the map (dimension is equal to number of functions)
  typedef typename VectorType::ScalarType ScalarType;
  typedef Function<VectorType> FunctionType;
  typedef capd::autodiff::DagIndexer<ScalarType> DAG;
  typedef capd::autodiff::Node NodeType;

  typedef BasicFunction<ScalarType> BaseFunction;
  typedef capd::diffAlgebra::Hessian<ScalarType,ImageVectorType::csDim,VectorType::csDim> HessianType;
  typedef capd::diffAlgebra::Jet<MatrixT,0> JetType;
  typedef typename BaseFunction::size_type size_type;

  Map();  ///< creates an univariate function R->R, f(x)=0
  Map(const std::string &, size_type degree=1); ///< parses expression from given string. Allocates memory for jet propagation of degree 'degree'
  Map(const char *, size_type degree=1);        ///< parses expression from given string. Allocates memory for jet propagation of degree 'degree'
  Map(const Map& f);                            ///< copying constructor

  /**
   * parses expression from given routine. Allocates memory for jet propagation of degree 'degree'
   * @param f  routine that defines function, its signature is: void (NodeType time, NodeType in[], int dimIn, NodeType out[], int dimOut, NodeType param[], int noParam)
   * @param dimIn  number of input variables (dimension of domain)
   * @param dimOut  number of output variables (dimension of codomain)
   * @param noParam  number of parameters
   * @param degree  maximal degree of derivatives that can be computed
  */
  template<typename Function>
  Map(Function f, int dimIn, int dimOut, int noParam, size_type degree=1);

  Map& operator=(const char *);                ///< parses expression from a given string and reallocates DAG, does not change already specified degree
  Map& operator=(const std::string &);         ///< parses expression from a given string and reallocates DAG, does not change already specified degree
  Map& operator=(const Map& f);                ///< assignment from an object

  /**
   * parses expression from given routine. Allocates memory for jet propagation of degree 'degree'
   * @param f - routine that defines function, its signature is: void (NodeType time, NodeType in[], int dimIn, NodeType out[], int dimOut, NodeType param[], int noParam)
   * @param dimIn - number of input variables (dimension of domain)
   * @param dimOut - number of output variables (dimension of codomain)
   * @param noParam - number of parameters
  */
  template<typename Function>
  void reset(Function f, int dimIn, int dimOut, int noParam, size_type degree);

  /**
   * parses expression from given routine. It does not change actual maximal degree of jets that moght be propagated by this map.
   * @param f - routine that defines function, its signature is: void (NodeType time, NodeType in[], int dimIn, NodeType out[], int dimOut, NodeType param[], int noParam)
   * @param dimIn - number of input variables (dimension of domain)
   * @param dimOut - number of output variables (dimension of codomain)
   * @param noParam - number of parameters
  */
  template<typename Function>
  void reset(Function f, int dimIn, int dimOut, int noParam);

  ImageVectorType operator()(const VectorType& u) const;                             ///< evaluates map at a given vector.
  ImageVectorType operator()(ScalarType t, const VectorType& u) const;               ///< evaluates at a given vector and for given time
  ImageVectorType operator()(const VectorType& u, MatrixType& out_derivative) const;                ///< computes simultaneously value and derivative of the map for a given vector
  ImageVectorType operator()(ScalarType t, const VectorType& u, MatrixType& out_derivative) const; ///< computes simultaneously value and derivative of the map for a given vector and time
  MatrixType operator[](const VectorType& u) const;                                 ///< computes derivative of the map for a given vector
  MatrixType derivative(const VectorType& u) const;                                 ///< computes derivative of the map for a given vector
  MatrixType derivative(ScalarType t, const VectorType& u) const;                   ///< computes derivative of the map for a given vector and time

  /**
   *  Computes derivative of the map at a vector 'u' that was send to any vector or function evaluating value.
   *  In this case computation can be speeded up.
   */
  void homogenousPolynomial(MatrixType& o_der) const;
  /**
   *  Computes derivative and hessian of the map at a vector 'u' that was send to any vector or function evaluating value.
   *  In this case computation can be speeded up.
   */
  void homogenousPolynomial(const MatrixType& o_der, HessianType& o_hessian) const;

  /**
   *  Computes jet to degree 'degree' of the map at a vector 'u', at which the map was previously evaluated.
   *  In this case computation can be speeded up.
   */
  template<class JetT>
  void homogenousPolynomial(JetT& x, size_type degree) const;

  ImageVectorType operator()(const VectorType& x, MatrixType& out_df, HessianType& out_hf) const;                ///< computes value, first and second order derivatives for a given vector
  ImageVectorType operator()(ScalarType t, const VectorType& x, MatrixType& out_df, HessianType& out_hf) const;  ///< computes value, first and second order derivatives for a given vector and time

  /**
   * computes propagation of the jet x through the map and returns jet of normalized partial derivatives of the function \f$ f\circ x \f$.
   * \note if \f$R^N->R^M\f$ with \f$N\neq M\f$ then input and output jets have different dimensions.
   */
  JetType operator()(const JetType& x) const;

  /**
   * computes propagation of the jet x through the map and returns jet of normalized partial derivatives of the function \f$ f\circ x \f$.
   * \note if \f$R^N->R^M\f$ with \f$N\neq M\f$ then input and output jets have different dimensions.
   */
  JetType operator()(ScalarType t, const JetType& x) const;

  ImageVectorType operator()(const VectorType& x, JetType& out_jet) const;
  ImageVectorType operator()(ScalarType t, const VectorType& x, JetType& out_jet) const;


  /**
   * iterative computation of Taylor coefficients up to given order for a solution to ODE represented by this map.
   * @param order - order of Taylor method (degree of polynomial approximation to the solution)
   * @param[in] coeffs[0] - initial condition for ODE must be the first vector in the table
   * @param[out] coeffs[1],...,coeffs[order] - computed coefficients of polynomial approximation
   */
  void computeODECoefficients(VectorType coeffs[], size_type order) const;

  /**
   * iterative computation of Taylor coefficients up to given order for a solution to ODE represented by this map together with first order variational equations
   * @param order - order of Taylor method (degree of polynomial approximation to the solution)
   * @param[in] coeffs[0] - initial condition for ODE must be the first vector in the table
   * @param[out] coeffs[1],...,coeffs[order] - computed coefficients of polynomial approximation
   * @param[in] dCoeffs[0] - initial condition for variational equation
   * @param[out] dCoeffs[1],...,dCoeffs[order] - computed coefficients of polynomial approximation to the solution to variational equation
   */
  void computeODECoefficients(VectorType coeffs[], MatrixType dCoeffs[], size_type order) const;

  /**
   * iterative computation of Taylor coefficients up to given order for a solution to ODE represented by this map together with first and second order variational equations
   * @param order - order of Taylor method (degree of polynomial approximation to the solution)
   * @param[in] coeffs[0] - initial condition for ODE must be the first vector in the table
   * @param[out] coeffs[1],...,coeffs[order] - computed coefficients of polynomial approximation
   * @param[in] dCoeffs[0] - initial condition for variational equation
   * @param[out] dCoeffs[1],...,dCoeffs[order] - computed coefficients of polynomial approximation to the solution to variational equation
   * @param[in] hCoeffs[0] - initial condition for second order variational equation
   * @param[out] hCoeffs[1],...,hCoeffs[order] - computed coefficients of polynomial approximation to the solution to second order variational equation
   * \note Second order variational equations are written in the base of normalized derivatives (Taylor coefficients). To obtain derivatives of solution wrt to initial condition one must multiply them by suitable factorials.
   */
  void computeODECoefficients(VectorType coeffs[], MatrixType dCoeffs[], HessianType hCoeffs[], size_type order) const; ///< iterative computation of Taylor coefficients up to given order for a solution to ODE, first and second order variational equations

  /**
   * Iterative computation of Taylor coefficients up to given order for a solution to ODE and its variational equations.
   * @param degree - The maximal order of derivatives with respect to initial conditions is specified by this argument
   * @param order - order of Taylor method (degree of polynomial approximation to the solution)
   * @param[in,out] coeffs[0] - jets of initial conditions for the solution and all variational equations,
   *                coeffs[1],...,coeffs[order] - computed coefficients  of polynomial approximation to the solution all variational equations
   * \note coeffs[i] can be a jet of degree bigger than the argument 'degree'.
   */
  void computeODECoefficients(JetType coeffs[], size_type degree, size_type order) const;

  // using BaseFunction::setParameters;
  void setParameters(const VectorType& values);   ///< simultaneously sets values of many parameters. It is assumed that given vector contains values of subsequent parameters.
  size_type degree() const { return m_degree; }    ///< returns maximal possible order of space-derivative that the map can compute. Specified by the user in the constructor (default value is 1).
  void setDegree(size_type degree);                ///< set maximal order of space-derivatives that the map can compute.

  /**
   * This method generates code of new class in which DAG is constructed as static array, dimensions are fixed giving the compiler chance
   * for most aggressive optimization.
   *
   * This method is used to speed up your computation.
   * DAG that represents an expression is allocated dynamically and thus its processing might be not that fast as in hand optimized code.
   *
   * @param[in] className - name of the class that will be generated. The code will be generated into file realativePath/className.h.
   * @param[in] relativePath - path to folder where generated code will be save
   * @param[in] userNamespace - must be specified - the code will be put into this namespace.
   *
   * \attention The method does not check whether arguments are correct identifiers in C++.
   */
  void codeTranslation(const char className[], const char userNamespace[], const char relativePath[]=".", size_type maxTimeOrder=32) const;
protected:
  void generateHeaderFile(const char className[], const char userNamespace[], const char relativePath[], size_type maxTimeOrder=32) const;
  void generateHppFile(const char className[], const char userNamespace[], const char relativePath[], size_type maxTimeOrder=32) const;
  void generateCppFile(const char className[], const char userNamespace[], const char relativePath[], size_type maxTimeOrder=32) const;

  void evalAndCopyResult(JetType& c) const;

  size_type m_degree;                  ///< maximal order of derivatives that the map can compute

  void checkDegree(size_type degree) const;
  void checkOrder(size_type order) const;
}; // the end of class Map


// ------------------------------------------------------------------------------

template<typename MatrixT>
template<typename Function>
inline
Map<MatrixT>::Map(Function f, int dimIn, int dimOut, int noParam, size_type degree)
  : BaseFunction(f,dimIn,dimOut,noParam), m_degree(degree)
{
  this->realloc(this->m_degree);
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
template<typename Function>
inline
void Map<MatrixT>::reset(Function f, int dimIn, int dimOut, int noParam, size_type degree)
{
  this->clean();
  BaseFunction::reset(f,dimIn,dimOut,noParam);
  this->realloc(degree);
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
template<typename Function>
inline
void Map<MatrixT>::reset(Function f, int dimIn, int dimOut, int noParam)
{
  this->reset(f,dimIn,dimOut,noParam,this->m_degree);
}

// -----------------------------------------------------------------------

template<typename MatrixT>
inline typename Map<MatrixT>::ImageVectorType
Map<MatrixT>::operator()(ScalarType t, const VectorType& u, MatrixType& der) const
{
  this->setCurrentTime(t);
  return this->operator()(u,der);
}

// -----------------------------------------------------------------------

template<typename MatrixT>
inline MatrixT Map<MatrixT>::operator[](const VectorType& v) const
{
  return derivative(v);
}

// -----------------------------------------------------------------------

template<typename MatrixT>
inline MatrixT Map<MatrixT>::derivative(const VectorType& v) const
{
  MatrixType der(this->imageDimension(),this->dimension());
  this->operator()(v,der);
  return der;
}

// -----------------------------------------------------------------------

template<typename MatrixT>
inline MatrixT Map<MatrixT>::derivative(ScalarType t, const VectorType& v) const
{
  MatrixType der(this->imageDimension(),this->dimension());
  this->operator()(t,v,der);
  return der;
}

// -----------------------------------------------------------------------

template<typename MatrixT>
inline typename Map<MatrixT>::ImageVectorType
Map<MatrixT>::operator()(ScalarType t, const VectorType& v, MatrixType& der, HessianType& h) const
{
  this->setCurrentTime(t);
  return this->operator()(v,der,h);
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
inline typename Map<MatrixT>::ImageVectorType
Map<MatrixT>::operator()(ScalarType t, const VectorType& v) const
{
  this->setCurrentTime(t);
  return this->operator()(v);
}

// -----------------------------------------------------------------------

template<typename MatrixT>
inline typename Map<MatrixT>::JetType
Map<MatrixT>::operator()(ScalarType t, const JetType& c) const
{
  this->setCurrentTime(t);
  return (*this)(c);
}

template<typename MatrixT>
inline typename Map<MatrixT>::ImageVectorType
Map<MatrixT>::operator()(ScalarType t, const VectorType& x, JetType& c) const{
  this->setCurrentTime(t);
  return (*this)(x,c);
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
inline void Map<MatrixT>::checkOrder(size_type order) const
{
  if(order>this->getOrder())
  {
    std::ostringstream s;
    s << "exception in class Map\n";
    s << "Cannot compute " << order << "-th Taylor coefficient\n";
    s << "The order is " << this->getOrder();
    throw std::out_of_range (s.str());
  }
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
inline void Map<MatrixT>::checkDegree(size_type degree) const
{
  if(degree>this->degree())
  {
    std::ostringstream s;
    s << "exception in class Map\n";
    s << "Cannot compute jet of degree " << degree<< "\n";
    s << "The maximal degree is " << this->degree();
    throw std::out_of_range (s.str());
  }
}

// ------------------------------------------------------------------------------

template<typename MatrixT>
template<class JetT>
void Map<MatrixT>::homogenousPolynomial(JetT& x, size_type degree) const
{
  this->checkDegree(degree);

  using namespace capd::autodiff;
  size_type i;

  for(i=0;i<this->dimension();++i)
  {
    typename JetType::const_iterator b = x.begin(i,degree-1), e=x.end(i,degree-1);
    typename DAG::iterator p = this->m_dag.begin(i,degree-1);
    while(b!=e)
    {
      *p = *b;
      b++;
      p += this->getOrder()+1;
    }
  }

  if(degree>1)
      this->evalHomogenousPolynomial(degree-1);
  this->evalHomogenousPolynomial(degree);

  for(i=0;i<this->dimension();++i)
  {
    typename JetType::iterator b = x.begin(i,degree), e = x.end(i,degree);
    typename DAG::iterator p = this->m_dag.begin(this->m_pos[i],degree);
    while(b!=e)
    {
      *b = *p;
      b++;
      p += this->getOrder()+1;
    }
  }
}


}} // the end of the namespace capd::map

#endif // _CAPD_MAP_MAP_H_

/// @}
