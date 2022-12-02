/// @addtogroup map
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file BasicFunction.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_MAP_BASICFUNCTION_HPP_
#define _CAPD_MAP_BASICFUNCTION_HPP_

#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <algorithm>

#include "capd/autodiff/DagIndexer.hpp"
#include "capd/map/BasicFunction.h"

namespace capd{
namespace map{

// -----------------------------------------------------------------------------

template<typename Scalar>
BasicFunction<Scalar>::BasicFunction()
{
  this->createDefault();
}

// -----------------------------------------------------------------------------

template<typename Scalar>
BasicFunction<Scalar>::BasicFunction(const std::string&s)
{
  this->createFromText(s);
}

// -----------------------------------------------------------------------------

template<typename Scalar>
BasicFunction<Scalar>::BasicFunction(const char* s)
{
  this->createFromText(s);
}

// -----------------------------------------------------------------------------

template<typename Scalar>
BasicFunction<Scalar>::BasicFunction(const BasicFunction& f)
{
  this->copyObject(f);
}

// -----------------------------------------------------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::operator=(const std::string& s)
{
  this->clean();
  this->createFromText(s);
}

// -----------------------------------------------------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::operator=(const char* s)
{
  this->clean();
  this->createFromText(s);
}

// -----------------------------------------------------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::operator=(const BasicFunction& f)
{
  this->clean();
  this->copyObject(f);
}

//--------------------- createDefault ----------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::createDefault()
{
  this->createFromText("var:x;fun:0;");
}

// -----------------------------------------------------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::createFromText(std::string s)
{
  this->clean();
  capd::map::removeWhiteSpaces(s);
  this->m_indexOfFirstParam = capd::map::parseVariables(s,this->m_var);
  parseMap(this->m_indexOfFirstParam,s,this->m_var,this->m_fullGraph,this->m_pos);
  this->m_evalPath.clear();
  for(unsigned i=0;i<this->m_fullGraph.size();++i)
    if( this->m_fullGraph[i].op!=capd::autodiff::NODE_NULL and
        this->m_fullGraph[i].op!=capd::autodiff::NODE_PARAM and
        this->m_fullGraph[i].op!=capd::autodiff::NODE_CONST and
        this->m_fullGraph[i].op!=capd::autodiff::NODE_TIME and
        this->m_fullGraph[i].op!=capd::autodiff::NODE_COS and
        this->m_fullGraph[i].op!=capd::autodiff::NODE_VAR
     ) this->m_evalPath.push_back(capd::autodiff::MyNode(this->m_fullGraph[i]));
  capd::autodiff::Int4ToAbstractNode(this->m_evalPath,this->m_nodes,this->m_dag);
}

// -----------------------------------------------------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::copyObject(const BasicFunction& f)
{
  this->m_fullGraph = f.m_fullGraph;
  this->m_evalPath = f.m_evalPath;
  this->m_pos = f.m_pos;
  this->m_var = f.m_var;
  this->m_dag = f.m_dag;
  this->m_indexOfFirstParam = f.m_indexOfFirstParam;
  capd::autodiff::Int4ToAbstractNode(this->m_evalPath,this->m_nodes,this->m_dag);
}

// --------------------- clean ---------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::clean()
{
  this->m_var.clear();
  this->m_pos.clear();
  this->m_fullGraph.clear();
  this->m_evalPath.clear();
  this->deleteNodes();
}

// -----------------------------------------------------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::deleteNodes()
{
  for(unsigned i=0;i<this->m_nodes.size();++i)
    delete this->m_nodes[i];
  this->m_nodes.clear();
}

// -----------------------------------------------------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::setCurrentTime(const ScalarType& a_time) const
{
  using namespace capd::autodiff;
  this->m_dag(VarNo(this->dimension()),CoeffNo(0)) = a_time;
}

// -----------------------------------------------------------------------------

template<typename Scalar>
const Scalar& BasicFunction<Scalar>::getCurrentTime() const
{
  using namespace capd::autodiff;
  return this->m_dag(VarNo(this->dimension()),CoeffNo(0));
}

// -----------------------------------------------------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::differentiateTime() const
{
  using namespace capd::autodiff;
  this->m_dag(VarNo(this->dimension()),CoeffNo(1)) = ScalarType(1.);
}

// -----------------------------------------------------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::setOrder(size_type order)
{
  if(this->getOrder() == order) return;
  this->m_dag.setOrder(order);
  for(unsigned i=0;i<this->m_evalPath.size();++i)
  {
    this->m_nodes[i]->left = this->m_dag.coefficients() + this->m_dag.timeJetSize()*this->m_evalPath[i].left;
    this->m_nodes[i]->right = this->m_dag.coefficients() + this->m_dag.timeJetSize()*this->m_evalPath[i].right;
    this->m_nodes[i]->result = this->m_dag.coefficients() + this->m_dag.timeJetSize()*this->m_evalPath[i].result;
  }
}

/* _________________________________________________________________ */

template<typename Scalar>
void BasicFunction<Scalar>::setParameter(const std::string &name, const ScalarType& value)
{
  using namespace capd::autodiff;
  for(unsigned i=this->m_indexOfFirstParam;i<this->m_var.size();++i)
  {
    if(this->m_var[i]!=name) continue;
    this->m_dag(VarNo(i),CoeffNo(0)) = value;
    return;
  }
}

// ------------------------------------------------------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::setParameter(size_type d, const Scalar& value)
{
  using namespace capd::autodiff;
  this->m_dag(VarNo(this->m_indexOfFirstParam+d+1),CoeffNo(0)) = value;
}


// ------------------------------------------------------------------------------

template<typename Scalar>
typename BasicFunction<Scalar>::ScalarType
BasicFunction<Scalar>::getParameter(size_type d) const
{
  using namespace capd::autodiff;
  const size_type vars = this->m_var.size();
  if(this->m_indexOfFirstParam + d +1< vars)
    return this->m_dag(VarNo(this->m_indexOfFirstParam+d+1),CoeffNo(0));
  throw std::runtime_error("BasicFunction::getParameter: incorrect index of parameter");
}

/* _________________________________________________________________ */

template<typename Scalar>
typename BasicFunction<Scalar>::ScalarType
BasicFunction<Scalar>::getParameter(const std::string &name) const
{
  using namespace capd::autodiff;
  for(unsigned i=this->m_indexOfFirstParam;i<this->m_var.size();++i)
  {
    if(this->m_var[i]!=name) continue;
    return this->m_dag(VarNo(i),CoeffNo(0));
  }
  throw std::runtime_error("BasicFunction::getParameter: udfeined name of parameter.");
}

// ------------------------------------------------------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::setParameters(const Scalar* values, size_type d)
{
  using namespace capd::autodiff;
  const size_type vars = this->m_var.size();
  if(this->m_indexOfFirstParam +1 + d != vars){
    throw std::runtime_error("BasicFunction::setParameters : incorrect number of parameters! ");
  }

  size_type cnt= this->m_indexOfFirstParam+1;
  size_type i=0;
  while(cnt < vars and i<d)
  {
      this->m_dag(VarNo(cnt),CoeffNo(0)) = values[i];
      ++cnt;
      ++i;
  }
}

// ------------------------------------------------------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::realloc(size_type degree)
{
  using namespace capd::autodiff;
  std::vector<Scalar> paramValues(this->m_fullGraph.size());

  bool copyParams = (this->m_fullGraph.size() == (unsigned)this->m_dag.numberOfNodes());
  if(copyParams)
    for(unsigned i=0;i<this->m_fullGraph.size();++i)
      paramValues[i] = this->m_dag(VarNo(i),CoeffNo(0));

  this->m_dag.resize(
      this->dimension(),
      this->imageDimension(),
      degree,
      this->m_fullGraph.size(),
      this->getOrder()
      );

 for(unsigned i=0;i<this->m_fullGraph.size();++i)
  {
    if(this->m_fullGraph[i].op == capd::autodiff::NODE_CONST)
      this->m_dag(VarNo(i),CoeffNo(0)) = this->m_fullGraph[i].val;
    else if(copyParams)
      this->m_dag(VarNo(i),CoeffNo(0)) = paramValues[i];
  }
 this->deleteNodes();
 capd::autodiff::Int4ToAbstractNode(this->m_evalPath,this->m_nodes,this->m_dag);
}

// -----------------------------------------------------------------------------

template<typename Scalar>
void BasicFunction<Scalar>::applyC1Mask() const
{
  using namespace capd::autodiff;
  if(this->getMask()){
    for(unsigned i=0;i<this->dimension();++i)
      if(!this->getMask(i))
        this->m_dag(VarNo(i),DerNo(i),CoeffNo(0)) = TypeTraits<ScalarType>::zero();
  }
}
}} // namespace capd::map

#endif // _CAPD_MAP_BASICFUNCTION_HPP_

/// @}
