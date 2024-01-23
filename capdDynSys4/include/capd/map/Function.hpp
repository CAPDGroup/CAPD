/// @addtogroup map
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Function.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_MAP_FUNCTION_HPP_
#define _CAPD_MAP_FUNCTION_HPP_

#include <algorithm>
#include <sstream>
#include <string>
#include <stdexcept>

#include "capd/map/BasicFunction.hpp"
#include "capd/map/Function.h"
#include "capd/map/Parser.h"
#include "capd/autodiff/eval.hpp"

namespace capd{
namespace map{

// ------------------------------------------------------------------------------

template<typename VectorType>
Function<VectorType>::Function(void)
{
  this->realloc(1);
}

// ------------------------------------------------------------------------------

template<typename VectorType>
Function<VectorType>::Function(const std::string &s) : BaseFunction(s)
{
  this->check(s);
  this->realloc(1);
}

// ------------------------------------------------------------------------------

template<typename VectorType>
Function<VectorType>::Function(const char *s) : BaseFunction(s)
{
  this->check(s);
  this->realloc(1);
}

// ------------------------------------------------------------------------------

template<typename VectorType>
Function<VectorType>::Function(const Function& f) : BaseFunction(f)
{
  this->realloc(1);
}

// ------------------------------------------------------------------------------

template<typename VectorType>
Function<VectorType>& Function<VectorType>::operator=(const std::string &s)
{
  try{
    BaseFunction::operator=(s);
    this->check(s);
    this->realloc(1);
  } catch(std::runtime_error &r)
  {
    BaseFunction::createDefault();
    this->createDefault();
    std::string re = "exception in Function &Function::operator=(const string &)\n";
    re += r.what();
    throw std::runtime_error(re);
  }
  return *this;
}

// ------------------------------------------------------------------------------

template<typename VectorType>
Function<VectorType>& Function<VectorType>::operator=(const char *s)
{
  return this->operator=(std::string(s));
}

// ------------------------------------------------------------------------------

template<typename VectorType>
Function<VectorType>& Function<VectorType>::operator=(const Function& f)
{
  if(&f==this) return *this;
  BaseFunction::operator=(f);
  return *this;
}

// ------------------------------------------------------------------------------

template<typename VectorType>
typename Function<VectorType>::ScalarType
   Function<VectorType>::operator()(const ScalarType& v) const
{
  using namespace capd::autodiff;
  this->m_dag(VarNo(0),CoeffNo(0)) = v;
  this->evalHomogenousPolynomial();
  return this->m_dag(VarNo(this->m_pos[0]),CoeffNo(0));
}

// ------------------------------------------------------------------------------

template<typename VectorType>
typename Function<VectorType>::ScalarType
Function<VectorType>::operator()(const VectorType &v) const
{
  using namespace capd::autodiff;
  unsigned d = v.dimension();
  if(d > this->m_var.size())
  {
    std::ostringstream message;
    message << "Exception in Function::operator(const VectorType &)\n";
    message << "VectorType has " << v.dimension() << " coordinates\n";
    message << "Function object has " << this->m_var.size() << " variables \n";

    throw std::runtime_error(message.str());
  }
  for(unsigned i=0;i<d;++i)
    this->m_dag(VarNo(i),CoeffNo(0)) = v[i];
  this->evalHomogenousPolynomial();
  return this->m_dag(VarNo(this->m_pos[0]),CoeffNo(0));
}

// -----------------------------------------------------------------------

template<typename VectorType>
VectorType Function<VectorType>::gradient(VectorType u) const
{
  using namespace capd::autodiff;
  this->setArgument(u);
  this->eval(1,CoeffNo(0));

  for(size_type i=0;i<u.dimension();++i)
    u[i] = this->m_dag(VarNo(this->m_pos[0]),DerNo(i),CoeffNo(0));
  return u;
}

// ------------------------------------------------------------------------------

template<typename VectorType>
void Function<VectorType>::setParameters(const VectorType & values)
{
  BaseFunction::setParameters(values.begin(),values.dimension());
}

// -----------------------------------------------------------------------

template<typename VectorType>
void Function<VectorType>::check(const std::string& s)
{
  if(this->m_pos.size()!=1)
  {
    std::ostringstream out;
    out << "Cannot parse a scalar valued function from expression: " << s;
    throw std::runtime_error(out.str());
  }
}


}} // namespace capd::map

#endif // _CAPD_MAP_FUNCTION_HPP_

/// @}
