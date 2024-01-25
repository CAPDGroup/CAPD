/////////////////////////////////////////////////////////////////////////////
/// @file SolverPool.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef __CAPD_THREADING_SOLVERPOOL_H__
#define __CAPD_THREADING_SOLVERPOOL_H__

#include "capd/threading/SolverFactory.h"

namespace capd{
namespace threading{
/// @addtogroup threading 
/// @{

template<class Solver>
struct SolverPool{
  SolverPool(unsigned poolSize, const SolverFactory<Solver>& factory){
    for(unsigned i=0;i<poolSize;++i)
      _solver.push_back(factory.createSolver());
  }
  ~SolverPool(){
    for(Solver* s : _solver) delete s;
  }

  void setOrder(int order){
    for(Solver* s : _solver) s->setOrder(order);
  }
  void setTolerance(double absTolerance, double relTolerance){
    for(Solver* s : _solver) {
      s->setAbsoluteTolerance(absTolerance);
      s->setRelativeTolerance(relTolerance);
    }
  }
  Solver& getSolver(int id) { return *_solver[id]; }
  std::vector<Solver*> _solver;
};

/// @}
}} // namespace capd::threading

#endif
