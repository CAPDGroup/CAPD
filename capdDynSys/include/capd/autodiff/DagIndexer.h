/////////////////////////////////////////////////////////////////////////////
/// @file DagIndexer.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_DAGINDEXER_H_
#define _CAPD_AUTODIFF_DAGINDEXER_H_

#include <vector>
#include "capd/basicalg/TypeTraits.h"
#include "capd/basicalg/factrial.h"
#include "capd/vectalg/ColumnVector.h"
#include "capd/diffAlgebra/CnContainer.h"

#include "c2_jet_indices_buffer.h"
#include "c3_jet_indices_buffer.h"

namespace capd{
namespace autodiff{

  inline int sumAndFindMax(const int a[], const int b[], int c[], const int n){
  int p=0,m=0;
  for(int i=0;i<n;++i)
  {
    c[i]=a[i]+b[i];
    if(c[i]>m){
      m=c[i];
      p=i;
    }
  }
  return p;
}

inline int findMax(const int c[], const int n){
  int p=0;
  for(int i=0;i<n;++i)
  {
    if(c[i]>c[p]){
      p=i;
    }
  }
  return p;
}

using capd::vectalg::__size_type;
using capd::vectalg::__difference_type;

// The code is written in almost pure C.
// Therefore there is a lot of integral arguments that can be easily incorrectly used.
// The following classes are wrappers for integral that assure
// correct order when passing arguments to functions.

// These classes should not be used in real computations, rather on development stage.
// Comment out the following line in order to use debug mode.
//#define _Dag_Indexer_Debug_Mode_

inline __size_type index(__size_type dim, __size_type j, __size_type c)
{
  return C2_Jet_Indices_Buffer<__size_type, 20>::compute_index(dim, j, c);
}

inline __size_type index(__size_type dim, __size_type j, __size_type c, __size_type k)
{
  return C3_Jet_Indices_Buffer<__size_type, 20>::compute_index(dim, j, c, k);
}

#ifdef _Dag_Indexer_Debug_Mode_

struct JetSize
{
  explicit JetSize(__size_type i) : m_i(i) {}

  inline operator __size_type (void) { return m_i; }
  __size_type m_i;
};

struct Order
{
  explicit Order(__size_type i) : m_i(i) {}

  inline operator __size_type (void) { return m_i; }
  __size_type m_i;
};

struct VarNo
{
  explicit VarNo(__size_type i) : m_i(i) {}

  inline operator __size_type (void) { return m_i; }
  __size_type m_i;
};

struct CoeffNo
{
  explicit CoeffNo(__size_type i) : m_i(i) {}

  inline operator __size_type (void) { return m_i; }
  __size_type m_i;
};

struct DerNo
{
  explicit DerNo(__size_type i) : m_i(i) {}

  inline operator __size_type (void) { return m_i; }
  __size_type m_i;
};

template<class ScalarType>
inline ScalarType& getC0Coeff(ScalarType* data, VarNo varNo, JetSize jetSize, CoeffNo coeffNo)
{
  return data[varNo*jetSize+coeffNo];
}

template<class ScalarType>
inline ScalarType& getC1Coeff(ScalarType* data, VarNo varNo, DerNo derNo, JetSize jetSize, Order order, CoeffNo coeffNo)
{
  return data[varNo*jetSize + (derNo+1)*(order+1) + coeffNo];
}

template<class ScalarType>
inline ScalarType& getC2Coeff(ScalarType* data, unsigned dim, VarNo varNo, DerNo j, DerNo c, JetSize jetSize, Order order, CoeffNo coeffNo)
{
  return data[varNo*jetSize + index(dim,j,c)*(order+1) + coeffNo];
}

template<class ScalarType>
inline ScalarType& getC3Coeff(ScalarType* data, unsigned dim, VarNo varNo, DerNo i, DerNo j, DerNo c, JetSize jetSize, Order order, CoeffNo coeffNo)
{
  return data[varNo*jetSize + index(dim,i,j,c)*(order+1) + coeffNo];
}

#else

  typedef __size_type VarNo;
  typedef __size_type DerNo;
  typedef __size_type CoeffNo;
  typedef __size_type JetSize;
  typedef __size_type Order;
  typedef __size_type Degree;
  typedef __size_type Dim;

  #define getC0Coeff(data,varNo,jetSize,coeffNo) data[varNo*jetSize+coeffNo]
  #define getC1Coeff(data,varNo,derNo,jetSize,order,coeffNo) data[varNo*jetSize + (derNo+1)*(order+1) + coeffNo]
  #define getC2Coeff(data,dim,varNo,j,c,jetSize,order,coeffNo) data[varNo*jetSize + index(dim,j,c)*(order+1) + coeffNo]
  #define getC3Coeff(data,dim,varNo,j,c,k,jetSize,order,coeffNo) data[varNo*jetSize + index(dim,j,c,k)*(order+1) + coeffNo]

#endif

/// Stores information about decomposition of a Multiinex 'z' into possible sums of x+y=z
/// Used to optimizs convolutions.
/// All the data here is redundant and precomputed to avoid extra runtime computation.
struct MultiindexData{
  typedef __size_type size_type;
  MultiindexData() : p(0),k(1), index(0) {}

  MultiindexData(capd::vectalg::Multiindex k, size_type order) :k(k){
    const size_type dim = k.dimension();
    const size_type deg = k.module();
    this-> p = findMax(k.begin(),dim);
    this->index = totalIndex(k,order+1,dim,deg);
    this->convolution.resize(order+1);
    this->convolutionFromEpToK.resize(order+1);
    capd::vectalg::Multiindex a(dim), b=k;
    if(this->index==0){
      for(unsigned coeffNo=0;coeffNo<=order;++coeffNo)
        for(unsigned j=0;j<=coeffNo;++j)
          convolution[coeffNo].push_back(IndexPair(j,coeffNo-j));
    } else {
      do{
        size_type ia = totalIndex(a,order+1,dim,deg);
        size_type ib = totalIndex(b,order+1,dim,deg);
        for(unsigned coeffNo=0;coeffNo<=order;++coeffNo)
          for(unsigned j=0;j<=coeffNo;++j){
            convolution[coeffNo].push_back(IndexPair(ia+j,ib+coeffNo-j));
            if(a[p]>0)
             convolutionFromEpToK[coeffNo].push_back(IndexPair(ia+j,ib+coeffNo-j));
          }
      }while(k.hasNext(a.begin(),b.begin()));
    }
  }

  static size_type totalIndex(const capd::vectalg::Multiindex& a, size_type order, size_type dim, size_type deg){
    size_type ma = a.module();
    return order*(a.index(deg) + (ma>0 ? binomial(dim+ma-1,dim) : 0));
  }
  size_type p; /// largest index in multiindex
  capd::vectalg::Multiindex k;
  size_type index; /// redundant data - index of k
  typedef std::pair<size_type,size_type> IndexPair;
  typedef std::vector< IndexPair > ConvolutionPairs;
  std::vector< ConvolutionPairs > convolution;
  std::vector< ConvolutionPairs > convolutionFromEpToK;

  const ConvolutionPairs& getConvolutionPairs(size_type coeffNo) const { return convolution[coeffNo]; }
  const ConvolutionPairs& getConvolutionPairsFromEpToK(size_type coeffNo) const { return convolutionFromEpToK[coeffNo]; }
};

template<class ScalarT>
class DagIndexer
{
public:

  typedef ScalarT ScalarType;
  typedef __size_type size_type;
  typedef capd::vectalg::ColumnVector<ScalarType,0> RefVectorType;
  typedef ScalarType* iterator;
  typedef const ScalarType* const_iterator;
  typedef capd::diffAlgebra::CnContainer<MultiindexData,0,0,0> IndexArray;

  DagIndexer(Dim domain=1, Dim image=1, Degree degree=1, size_type nodes=1, Order order=0);
  DagIndexer(const DagIndexer& dag);
  ~DagIndexer();

  DagIndexer& operator=(const DagIndexer& dag);
/*
  template<class iterator>
  inline void setVector(iterator b, iterator e)
  {
    ScalarType *c = m_coefficients;
    while(b!=e)
    {
      *c = *b;
      ++b;
      c += jetSize();
    }
  }
*/

  ScalarType& operator()(VarNo varNo,CoeffNo coeffNo) {  return getC0Coeff(m_coefficients,varNo,JetSize(m_timeJetSize),coeffNo); }
  ScalarType& operator()(VarNo varNo, DerNo derNo, CoeffNo coeffNo) { return getC1Coeff(m_coefficients,varNo,derNo,JetSize(m_timeJetSize),m_order,coeffNo); }
  ScalarType& operator()(VarNo varNo, DerNo j, DerNo c, CoeffNo coeffNo) { return getC2Coeff(m_coefficients,m_domainDimension,varNo,j,c,JetSize(m_timeJetSize),m_order,coeffNo); }
  ScalarType& operator()(VarNo varNo, DerNo i, DerNo j, DerNo c, CoeffNo coeffNo) { return getC3Coeff(m_coefficients,m_domainDimension,varNo,i,j,c,JetSize(m_timeJetSize),m_order,coeffNo); }

  Dim domainDimension() const { return m_domainDimension; }
  Dim imageDimension()  const { return m_imageDimension; }
  Dim degree()          const { return m_degree; }
  JetSize jetSize()          const { return JetSize(binomial(m_domainDimension+m_degree,m_degree)); }
  JetSize timeJetSize()      const { return JetSize(m_timeJetSize); }

  /**
   * This method defines a mask for computation of partial derivatives of the function represented by the instance.
   * Each element of the range [b,e) should be a valid Multiindex. The user can specify which partial derivatives he/she needs tp compute.
   * Dependent derivatives are added to the list automatically and those independent are not evaluated which significantly speeds up the computation.
   *
   * Example:
   * setMask({Multiindex({1,1,0}),Multiindex({2,0,0})});
   *
   * Here we request derivatives dx1dx2 and d^2x1. They depend on first order derivatives dx1 and dx2 which will be added automatically.
   *
   * @param [b,e) - iterator range of Multiindxes
   */
  template<class Iterator>
  void setMask(Iterator b, Iterator e);
  const bool* getMask() const { return m_mask; }
  bool getMask(size_type j) const { return getC1Coeff(m_mask,VarNo(0),DerNo(j),JetSize(m_timeJetSize),m_order,CoeffNo(0)); }
  bool getMask(size_type j, size_type c) const { return getC2Coeff(m_mask,m_domainDimension,VarNo(0),DerNo(j),DerNo(c),JetSize(m_timeJetSize),m_order,CoeffNo(0)); }
  void addMultiindexToMask(const capd::vectalg::Multiindex& i);
  void resetMask();

  ScalarType* coefficients()             { return m_coefficients;}
  const ScalarType* coefficients() const { return m_coefficients;}
  Order getOrder()              const { return m_order; }
  void setOrder(Order order);
  void resize(Dim domain, Dim image, Degree degree, size_type nodes, Order order);
  size_type numberOfNodes() const {return m_numberOfNodes;} ///< returns total number of nodes in DAG representing expression
  iterator begin();               ///< iterator selection. Returns iterator to the first element in container
  iterator end();                 ///< iterator selection. Returns iterator to the first element in container
  iterator begin(size_type i);          ///< iterator selection. Returns iterator to the first coefficient of the i-th component
  iterator end(size_type i);            ///< iterator selection. Returns iterator to an element after the last element the i-th component
  iterator begin(size_type i, size_type d);   ///< iterator selection. Returns iterator to the first coefficient of the i-th component of the homogeneous part of degree 'd'
  iterator end(size_type i, size_type d);     ///< iterator selection. Returns iterator to an element after the last coefficient of the i-th component of the homogeneous part of degree 'd'

  const_iterator begin() const;                 ///< iterator selection. Returns iterator to the first element in container
  const_iterator end() const;                   ///< iterator selection. Returns iterator to the first element in container
  const_iterator begin(size_type i) const;            ///< iterator selection. Returns iterator to the first coefficient of the i-th component
  const_iterator end(size_type i) const;              ///< iterator selection. Returns iterator to an element after the last element the i-th component
  const_iterator begin(size_type i, size_type d) const;     ///< iterator selection. Returns iterator to the first coefficient of the i-th component of the homogeneous part of degree 'd'
  const_iterator end(size_type i, size_type d) const;       ///< iterator selection. Returns iterator to an element after the last coefficient of the i-th component of the homogeneous part of degree 'd'

  const IndexArray& getIndexArray() const{
    return this->m_indexArray;
  }
private:
  void add(const capd::vectalg::Multiindex& i);
  void fillByZeroes();

  /// allocates memory when all parameters are known. All coefficients are set to zero.
  void allocate(Dim domain, Dim image, Degree degree, size_type nodes, Order order);

  /// allocates memory and copies data from an existing object.
  void allocateAndCopy(const DagIndexer& dag);

  /// precomputes arrays of indices for convolutions
  void createIndexArray();

  ScalarType* m_coefficients; ///< pointer to allocated memory
  Dim m_domainDimension;      ///< total dimension of the domain (with time, parameters, etc)
  Dim m_imageDimension;       ///< dimension of the counterdomain
  Degree m_degree;            ///< degree of jet (polynomial of space variables)

  size_type m_numberOfNodes;  ///< total number of nodes in DAG
  Order m_order;              ///< order of polynomial with respect to distinguished time variable

  JetSize m_timeJetSize;      ///< size of chunk of memory needed for (m_order+1) jets.
  IndexArray m_indexArray;    ///< array of precomputed pairs of indexes for computation of convolutions
  bool* m_mask;               ///< a pointer to mask of derivatives
};

// ------------------- iterator selections --------------------------

template<class T>
inline typename DagIndexer<T>::iterator DagIndexer<T>::begin(){
  return this->m_coefficients;
}

template<class T>
inline typename DagIndexer<T>::iterator DagIndexer<T>::end(){
  return this->m_coefficients + this->m_numberOfNodes*this->timeJetSize();
}

template<class T>
inline typename DagIndexer<T>::iterator DagIndexer<T>::begin(size_type i){
  return this->m_coefficients + i*this->timeJetSize();
}

template<class T>
inline typename DagIndexer<T>::iterator DagIndexer<T>::end(size_type i){
  return this->m_coefficients + (i+1)*this->timeJetSize();
}

template<class T>
inline typename DagIndexer<T>::iterator DagIndexer<T>::begin(size_type i, size_type d){
  return d==0
         ? this->begin(i)
         : this->begin(i) + (this->getOrder()+1)*binomial(this->m_domainDimension+d-1,d-1);
}

template<class T>
inline typename DagIndexer<T>::iterator DagIndexer<T>::end(size_type i, size_type d){
  return this->begin(i) + (this->getOrder()+1)*binomial(this->m_domainDimension+d,d);
}

// ------------------- const_iterator selections --------------------------

template<class T>
inline typename DagIndexer<T>::const_iterator DagIndexer<T>::begin() const{
  return this->m_coefficients;
}

template<class T>
inline typename DagIndexer<T>::const_iterator DagIndexer<T>::end() const {
  return this->m_coefficients + this->m_numberOfNodes*this->timeJetSize();
}

template<class T>
inline typename DagIndexer<T>::const_iterator DagIndexer<T>::begin(size_type i) const{
  return this->m_coefficients + i*this->timeJetSize();
}

template<class T>
inline typename DagIndexer<T>::const_iterator DagIndexer<T>::end(size_type i) const{
  return this->m_coefficients + (i+1)*this->timeJetSize();
}

template<class T>
inline typename DagIndexer<T>::const_iterator DagIndexer<T>::begin(size_type i, size_type d) const{
  return d==0
         ? this->begin(i)
         : this->begin(i) + (this->getOrder()+1)*binomial(this->m_domainDimension+d-1,d-1);
}

template<class T>
inline typename DagIndexer<T>::const_iterator DagIndexer<T>::end(size_type i, size_type d) const{
  return this->begin(i) + (this->getOrder()+1)*binomial(this->m_domainDimension+d,d);
}

template<class T>
template<class Iterator>
void DagIndexer<T>::setMask(Iterator b, Iterator e){
  if(this->m_mask)
    delete[] m_mask;
  this->m_mask = new bool[this->m_timeJetSize];
  std::fill(this->m_mask,this->m_mask+this->m_order+1,true);
  std::fill(this->m_mask+this->m_order+1,this->m_mask+this->m_timeJetSize,false);

  while(b!=e){
    this->add(*b);
    ++b;
  }
  this->fillByZeroes();
}

}} // namespace capd::autodiff

#endif
