/////////////////////////////////////////////////////////////////////////////
/// @file NonlinearSection.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2013 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_NONLINEAR_SECTION_H_
#define _CAPD_POINCARE_NONLINEAR_SECTION_H_

#include <string>
#include "capd/map/Function.h"
#include "capd/poincare/AbstractSection.h"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
/**
*  TimeMap class provides class that serves as Poincare section of the form x_i = c.
*  The section is defined by:
*  - integer - index of variable
*  - constant c that defines affine hyperplane x_i=c
*/

template<typename MatrixT>
class NonlinearSection : public AbstractSection<MatrixT>, public capd::map::Function<typename MatrixT::RowVectorType>
{
public:
  typedef MatrixT                   MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename VectorType::size_type size_type;             ///< integral type used to index containers (vectors, matrices, etc)
  typedef capd::dynset::AbstractSet<VectorType> Set;   ///< type of abstract base class for all sets
  typedef typename  AbstractSection<MatrixT>::JetType JetType;

  typedef capd::map::Function<VectorType> BaseFunction;

  ///< parses expression from given string. Allocates memory for ODE solver of order 'order'
  NonlinearSection(const std::string& s) : BaseFunction(s) {}
  ///< parses expression from given string. Allocates memory for ODE solver of order 'order'
  NonlinearSection(const char* s) : BaseFunction(s) {}

  /**
   * parses expression from given routine. Allocates memory for jet propagation of degree 'degree'
   * @param f - routine that defines function, its signature is: void (NodeType time, NodeType in[], int dimIn, NodeType out[], int dimOut, NodeType param[], int noParam)
   * @param dimIn - number of input variables (dimension of domain)
   * @param dimOut - number of output variables (dimension of codomain)
   * @param noParam - number of parameters
  */
  template<typename Functional>
  NonlinearSection(Functional f, int dimIn, int noParam, size_type degree=1)
    : BaseFunction(f,dimIn,noParam,degree)
  {}

  ///< substitution from string. Parse expression and compute new DAG.
  NonlinearSection& operator=(const std::string& s){
    BaseFunction::operator=(s);
    return *this;
  }
  ///< substitution from string. Parse expression and compute new DAG.
  NonlinearSection& operator=(const char* s){
    BaseFunction::operator=(s);
    return *this;
  }

  ScalarType operator()(const VectorType& v) const{
    return BaseFunction::operator()(v);
  }

  VectorType gradient(const VectorType& v) const{
    return BaseFunction::gradient(v);
  }

  ScalarType evalAt(const capd::dynset::AbstractSet<VectorType>& s) const{
    return s.evalAt(*this);
  }

  using AbstractSection<MatrixT>::computeDT; // virtual function
  void computeDT(const JetType& /*Px*/, const JetType& /*vfOnPx*/, JetType& /*dT*/, size_type) const
  {
     throw std::runtime_error("NonlinearSection::computeDT is not implemented for jets");
  }

}; // end of template NonlinearSection

/// @}
}} // namespace capd::poincare

#endif // _CAPD_POINCARE_NONLINEAR_SECTION_H_

