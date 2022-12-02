/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file DynSys.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_MOVE_H_
#define _CAPD_DYNSYS_MOVE_H_

#include "capd/vectalg/Norm.h"
#include "capd/basicalg/TypeTraits.h"
#include "capd/dynset/SetTraits.h"
#include "capd/diffAlgebra/CoeffTraits.h"

namespace capd{
namespace dynsys{

template<class T, class SetT,
  bool isSet =  capd::dynset::SetTraits<typename SetT::SetType>::isC0Set
             or capd::dynset::SetTraits<typename SetT::SetType>::isC1Set
           >
struct C1SetMove{
	  static void move(SetT& set, T& solver){
		  set.move(solver);
	  }
    static void move(SetT& set, SetT& result, T& solver){
		  set.move(solver, result);
	  }
};

// this excludes C2Set, CnSet as well as attempts of sending IVector as an argument to Taylor(HOE) solver.
template<class T, class SetT>
struct C1SetMove<T,SetT,false>{
	  static void move(SetT&, T&){
		  throw std::logic_error("Logic error: This solver can integrate C0Sets and C1Sets, only.\n");
	  }
    static void move(SetT&, SetT&, T&){
		  throw std::logic_error("Logic error: This solver can integrate C0Sets and C1Sets, only.\n");
	  }
};

template<class T, class SetT,
  bool isSet =  capd::dynset::SetTraits<typename SetT::SetType>::isC0Set
             or capd::dynset::SetTraits<typename SetT::SetType>::isC1Set
             or capd::dynset::SetTraits<typename SetT::SetType>::isC2Set
           >
struct C2SetMove{
	  static void move(SetT& set, T& solver){
		  set.move(solver);
	  }
    static void move(SetT& set, SetT& result, T& solver){
		  set.move(solver, result);
	  }
};

// this excludes CnSet as well as attempts of sending IVector as an argument to C2Taylor solver.
template<class T, class SetT>
struct C2SetMove<T,SetT,false>{
	  static void move(SetT& , T&){
		  throw std::logic_error("Logic error: This solver can integrate C0Sets, C1Sets and C2Sets, only.\n");
	  }
    static void move(SetT&, SetT&, T&){
		  throw std::logic_error("Logic error: This solver can integrate C0Sets, C1Sets and C2Sets, only.\n");
	  }
};

template<class T, class SetT,
  bool isSet =  capd::dynset::SetTraits<typename SetT::SetType>::isC0Set
             or capd::dynset::SetTraits<typename SetT::SetType>::isC1Set
             or capd::dynset::SetTraits<typename SetT::SetType>::isC2Set
             or capd::dynset::SetTraits<typename SetT::SetType>::isCnSet
           >
struct CnSetMove{
	  static void move(SetT& set, T& solver){
		  set.move(solver);
	  }
    static void move(SetT& set, SetT& result, T& solver){
		  set.move(solver, result);
	  }
};

// this excludes attempts of sending IVector as an argument to C2Taylor solver.
template<class T, class SetT>
struct CnSetMove<T,SetT,false>{
	  static void move(SetT&, T&){
		  throw std::logic_error("Logic error: This solver can integrate sets only. Do not use it for nonrigorous computations.\n");
	  }
    static void move(SetT&, SetT&, T&){
		  throw std::logic_error("Logic error: This solver can integrate sets only. Do not use it for nonrigorous computations.\n");
	  }
};

// #########################################

template<class T, class JetT, bool isC0JetOrC2Jet = capd::diffAlgebra::CoeffTraits<JetT>::isC0Jet or capd::diffAlgebra::CoeffTraits<JetT>::isC1Jet >
struct C1JetMove{
	  static void move(JetT& Jet, T& solver){
		  Jet.move(solver);
	  }
};

// this excludes C2Jet, CnJet
template<class T, class JetT>
struct C1JetMove<T,JetT,false>{
	  static void move(JetT&, T& ){
		  throw std::logic_error("Logic error: This solver can integrate C0Jets and C1Jets, only.\n");
	  }
};

template<class T, class JetT, bool isC0JetOrC2Jet = capd::diffAlgebra::CoeffTraits<JetT>::isC0Jet or capd::diffAlgebra::CoeffTraits<JetT>::isC1Jet or capd::diffAlgebra::CoeffTraits<JetT>::isC2Jet >
struct C2JetMove{
	  static void move(JetT& Jet, T& solver){
		  Jet.move(solver);
	  }
};

// this excludes CnJet
template<class T, class JetT>
struct C2JetMove<T,JetT,false>{
	  static void move(JetT&, T&){
		  throw std::logic_error("Logic error: This solver can integrate C0Jets, C1Jets and C2Jets, only.\n");
	  }
};

template<class T, class JetT, bool isC0JetOrC2Jet = capd::diffAlgebra::CoeffTraits<JetT>::isC0Jet or capd::diffAlgebra::CoeffTraits<JetT>::isC1Jet or capd::diffAlgebra::CoeffTraits<JetT>::isC2Jet or capd::diffAlgebra::CoeffTraits<JetT>::isCnJet >
struct CnJetMove{
	  static void move(JetT& Jet, T& solver){
		  Jet.move(solver);
	  }
};

template<class T, class JetT>
struct CnJetMove<T,JetT,false>{
	  static void move(JetT&, T&){
		  throw std::logic_error("Logic error: This solver can integrate Jets. Do not use it for rigorous computations.\n");
	  }
};

}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_DYNSYS_H_

/// @}
