/////////////////////////////////////////////////////////////////////////////
/// @file C0GraphicalSet.h
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#ifndef _CAPD_DYNSET_C0GRAPHICALSET_H_
#define _CAPD_DYNSET_C0GRAPHICALSET_H_

#include "capd/dynset/C0Set.h"
#include "capd/dynsys/DynSys.h"


namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

///////////////////////////////////////////////////////////////////////////////////
/// C0GraphicalSet is an envelope class for any class derived from C0Set.
/// It adds a possibility of an additional Output after each 'move' of the original set.
///
///  Output only needs to implement function
///    void show(C0Set & set)
///  which can e.g. draws on a screen or logs the set position to a file.
//////////////////////////////////////////////////////////////////////////////////
  
template<typename BaseSetT, typename OutputT>
class C0GraphicalSet : public BaseSetT{
public:
  typedef BaseSetT BaseSet;
  typedef typename BaseSet::MatrixType MatrixType;
  typedef typename BaseSet::VectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;
  typedef typename BaseSet::SetType SetType;
  typedef typename C0Set<MatrixType>::DynSysType DynSysType;
  typedef OutputT Output;

  C0GraphicalSet(const BaseSet & set, Output * output = 0)
    : BaseSet(set),  
      m_output(output) {
  }
  
  C0GraphicalSet(const BaseSet & set, Output & output)
    : BaseSet(set),  
      m_output(&output) {
  }
  C0GraphicalSet(const C0GraphicalSet &c) : BaseSet(c){
    m_output = c.m_output;
  }
  
  void move(DynSysType & dynsys){
    BaseSet::move(dynsys);
    if(m_output) m_output->show(*this);
  }

  void move(DynSysType & dynsys, C0GraphicalSet & result){
    //ResultType* set = dynamic_cast<ResultType*>(m_set);
    BaseSet::move(dynsys,result);
    if(m_output) m_output->show(result);
  }

//  VectorType affineTransformation(const MatrixType& A, const VectorType& v) const {
//    return BaseSet::affineTransformation(A, v);
//  }

//  std::string show(void) const{
//    return m_set->show();
//  }

//  operator VectorType() const{
//    return static_cast<VectorType>(*m_set);
//  }

  C0GraphicalSet &operator = (const VectorType &v){
    this->BaseSet::operator=(v);
    return *this;
  }

//  C0GraphicalSet &operator=(const C0GraphicalSet &S){
//    m_set = S.m_set;
//    m_output = S.m_output;
//    return *this;
//  }
  
//  SetType & getSet(){
//    return *m_set;
//  }
  
  Output & getOutput(){
    return *m_output();
  }
  Output & setOutput(Output & output){
    Output * old_output = m_output;
    m_output = output;
    return *old_output;
  }
protected:
  Output * m_output;
};
/// @}

}} // the end of the namespace capd::dynset

#endif // _CAPD_DYNSET_C0GRAPHICALSET_H_

