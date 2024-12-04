/////////////////////////////////////////////////////////////////////////////
/// @file DagIndexer.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_DAGINDEXER_HPP_
#define _CAPD_AUTODIFF_DAGINDEXER_HPP_

#include <algorithm>
#include "capd/autodiff/DagIndexer.h"
#include "capd/vectalg/Container.hpp"

namespace capd{
namespace autodiff{

template<class Scalar>
DagIndexer<Scalar>::DagIndexer(Dim domain, Dim image, Degree degree, size_type nodes, Order order)
  : m_order(order), m_indexArray(1,domain,degree,MultiindexData()), m_mask(0)
{
  this->allocate(domain,image,degree,nodes,order);
}

// -----------------------------------------------------------------------------

template<class Scalar>
DagIndexer<Scalar>::DagIndexer(const DagIndexer& dag)
  : m_order(dag.m_order), m_indexArray(dag.m_indexArray), m_mask(0)
{
  this->allocateAndCopy(dag);
}

// -----------------------------------------------------------------------------

template<class Scalar>
DagIndexer<Scalar>& DagIndexer<Scalar>::operator=(const DagIndexer& dag)
{
  if(this == &dag) return *this;
  delete[] m_coefficients;
  resetMask();
  this->allocateAndCopy(dag);
  return *this;
}

// -----------------------------------------------------------------------------

template<class Scalar>
DagIndexer<Scalar>::~DagIndexer()
{
  delete[] m_coefficients;
  if(this->m_mask)
    delete[] m_mask;
}

// -----------------------------------------------------------------------------

template<class Scalar>
void DagIndexer<Scalar>::setOrder(Order order)
{
  size_type newTimeJetSize = (order+1)*this->jetSize();
  size_type blockSize = newTimeJetSize*m_numberOfNodes;
  ScalarType* coeff = new ScalarType[blockSize];
  std::fill(coeff,coeff+blockSize,TypeTraits<ScalarType>::zero());

  for(size_type i=0;i<this->m_numberOfNodes;++i)
    getC0Coeff(coeff,VarNo(i),JetSize(newTimeJetSize),CoeffNo(0)) = getC0Coeff(m_coefficients,VarNo(i),JetSize(m_timeJetSize),CoeffNo(0));
  delete[] m_coefficients;
  m_coefficients = coeff;

  if(m_mask){
    bool* new_mask = new bool[newTimeJetSize];
    for(size_type i=0;i<this->jetSize();++i){
      for(size_type j=0;j<=order;++j)
        new_mask[i*(order+1)+j] = m_mask[i*(m_order+1)];
    }
    delete[] m_mask;
    m_mask = new_mask;
  }

  this->m_order = Order(order);
  this->m_timeJetSize = newTimeJetSize;
  this->createIndexArray();
}

// -----------------------------------------------------------------------------

template<class Scalar>
void DagIndexer<Scalar>::resize(Dim domain, Dim image, Degree degree, size_type nodes, Order order)
{
  delete[] m_coefficients;
  resetMask();
  this->allocate(domain,image,degree,nodes,order);
}

// -----------------------------------------------------------------------------

template<class Scalar>
void DagIndexer<Scalar>::createIndexArray(){
  m_indexArray = IndexArray(1,m_domainDimension,m_degree,MultiindexData());
  m_indexArray[0] = MultiindexData(capd::vectalg::Multiindex(m_domainDimension),m_order);
  for(size_type i=1;i<=m_degree;++i)
  {
    capd::vectalg::Multipointer a = m_indexArray.first(i);
    do{
      capd::vectalg::Multiindex mi(m_domainDimension,a);
      m_indexArray(0,a) = MultiindexData(mi,m_order);
    }while(m_indexArray.hasNext(a));
  }
}

// -----------------------------------------------------------------------------

template<class Scalar>
void DagIndexer<Scalar>::allocate(Dim domain, Dim image, Degree degree, size_type nodes, Order order)
{
  this->m_domainDimension = domain;
  this->m_imageDimension = image;
  this->m_degree = degree;
  this->m_numberOfNodes = nodes;
  this->m_order = Order(order);
  this->m_timeJetSize = (this->m_order+1)*this->jetSize();
  size_type blockSize = this->m_numberOfNodes*this->m_timeJetSize;
  this->m_coefficients = new ScalarType[blockSize];
  std::fill(this->m_coefficients,this->m_coefficients+blockSize,TypeTraits<ScalarType>::zero());
  this->createIndexArray();
}

// -----------------------------------------------------------------------------

template<class Scalar>
void DagIndexer<Scalar>::allocateAndCopy(const DagIndexer& dag)
{
  this->m_domainDimension = dag.m_domainDimension;
  this->m_imageDimension = dag.m_imageDimension;
  this->m_degree = dag.m_degree;
  this->m_numberOfNodes = dag.m_numberOfNodes;
  this->m_order = dag.m_order;
  this->m_timeJetSize = dag.m_timeJetSize;
  this->m_indexArray = dag.m_indexArray;

  size_type blockSize = this->m_numberOfNodes*this->m_timeJetSize;
  this->m_coefficients = new ScalarType[blockSize];
  std::copy(dag.coefficients(),dag.coefficients()+blockSize,this->m_coefficients);

  if(dag.m_mask){
    this->m_mask = new bool[this->m_timeJetSize];
    std::copy(dag.m_mask,dag.m_mask+dag.m_timeJetSize,this->m_mask);
  }
}

template<class Scalar>
void DagIndexer<Scalar>::addMultiindexToMask(const capd::vectalg::Multiindex& mi){
  this->add(mi);
  this->fillByZeroes();
}

template<class Scalar>
void DagIndexer<Scalar>::add(const capd::vectalg::Multiindex& mi){
  const MultiindexData& m = this->getIndexArray()(0,mi);
  const MultiindexData::ConvolutionPairs& p = m.getConvolutionPairs(this->m_order);
  for(MultiindexData::ConvolutionPairs::const_iterator i=p.begin();i!=p.end();++i){
    this->m_mask[i->first] = true;
    this->m_mask[i->second] = true;
  }
}

template<class Scalar>
void DagIndexer<Scalar>::fillByZeroes(){
  for(size_type i=0;i<this->m_numberOfNodes;++i)
    std::fill(this->m_coefficients+i*this->timeJetSize()+1,this->m_coefficients+(i+1)*this->timeJetSize(),TypeTraits<ScalarType>::zero());
}

template<class Scalar>
void DagIndexer<Scalar>::resetMask(){
  if(m_mask){
    delete[] m_mask;
    m_mask = 0;
  }
}

}}

#endif
