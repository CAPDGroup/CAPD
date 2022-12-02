/////////////////////////////////////////////////////////////////////////////
/// @file C11Rect2Set.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#ifndef _CAPD_DYNSET_C11RECT2SET_HPP_
#define _CAPD_DYNSET_C11RECT2SET_HPP_

#include <sstream>
#include "capd/basicalg/minmax.h"
#include "capd/dynset/C11Rect2Set.h"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"

namespace capd{
namespace dynset{


template<typename MatrixType>
void C11Rect2Set<MatrixType>::move(DynSysType& dynsys, C11Rect2Set& result) const
{
  // important: here we assume that both m_r and m_r0 contains zero
  // this is assured by each constructor and each step of this algorithm

  int dim = this->m_x.dimension();
  VectorType y(dim), rem(dim), enc(dim);
  MatrixType jacPhi(dim,dim), jacEnc(dim,dim), jacRem(dim,dim);
  MatrixType B(dim,dim), Q(dim,dim);

  VectorType xx = VectorType(*this);

  // the following function can throw an exception leaving output parameters in an inconsistent state
  // do not overwrite parameters of the set until we are sure that they are computed correctly
  dynsys.encloseC1Map(this->getCurrentTime(),
                    this->m_x, xx,          // input parameters
                    y, rem, enc,            // C^0 output
                    jacPhi, jacRem, jacEnc  // C^1 output
                    );

  jacPhi += jacRem;
  result.m_x = y + dynsys.Remainder(this->getCurrentTime(),this->m_x,rem);
  result.m_C = jacPhi * this->m_C;
  B = result.m_B = jacPhi * this->m_B;

  // ---------- C^0 - part -----------------

  // here we compute enclosure of the image after one iteration of the map/flow
  result.m_currentSet = result.m_x + result.m_C * this->m_r0 + result.m_B*this->m_r;

  // here we compute representation for the new set
  // xx is unnecessary now
  split(result.m_x, xx);
  split(result.m_C, Q);
  xx += Q * this->m_r0;

  // we assume that Policies provides algorithms for computation
  // of B, its inverse invB and updates
  this->computeBinvB(result.m_B,result.m_invB,this->m_r);

  // eventually we compute new representation of r
  result.m_r = (result.m_invB * B) * this->m_r + result.m_invB * xx;

  // ---------- C^1 - part -----------------
  result.m_D = jacPhi*this->m_D;
  B = result.m_Bjac = jacPhi*this->m_Bjac;
  result.m_Cjac = jacPhi*this->m_Cjac;

  // here we compute enclosure of the image after one iteration of the map/flow
  result.m_currentMatrix = jacPhi*this->m_currentMatrix;
  intersection(result.m_currentMatrix,result.m_D+result.m_Cjac*this->m_R0 + result.m_Bjac*this->m_R,result.m_currentMatrix);

  // here we compute representation for the new set
  // jacRem is unnecessary now
  split(result.m_D, jacRem);
  split(result.m_Cjac, Q);
  jacRem += Q * this->m_R0;

  // we assume that Policies provides algorithms for computation
  // of B, its inverse invB and updates
  this->computeBinvB(result.m_Bjac,result.m_invBjac,this->m_R);

  // eventually we compute new representation of r
  result.m_R = (result.m_invBjac * B) * this->m_R + result.m_invBjac * jacRem;

  if(&result != this){
    result.m_r0 = this->m_r0;
    result.m_R0 = this->m_R0;
  }

  result.setCurrentTime(this->getCurrentTime()+dynsys.getStep());
  result.setLastEnclosure(enc);
  result.setLastMatrixEnclosure(jacEnc);
  //this->Policies::reorganizeC0IfNeeded((C0BaseSet&)result);
  //this->Policies::reorganizeC1IfNeeded((C1BaseSet&)result);
}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C11RECT2SET_HPP_
