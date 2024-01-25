/////////////////////////////////////////////////////////////////////////////
/// @file SolverFactory.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef __CAPD_THREADING_SOLVERFACTORY_H__
#define __CAPD_THREADING_SOLVERFACTORY_H__
#include <string>

namespace capd{
namespace threading{
/// @addtogroup threading 
/// @{

/**
 This is an interface of an abstract factory wich creates new instances of a specific solver;
 */
template<class Solver>
struct SolverFactory{
  virtual Solver* createSolver() const = 0;
};

template<class Solver>
struct DefaultSolverFactory : SolverFactory<Solver>{
  virtual Solver* createSolver() const { return new Solver(); }
};

template<class S, class VF>
struct BaseTPMap{
  VF vectorField;
  S solver;
  BaseTPMap(const VF& vf, int order)
    : vectorField(vf),
      solver(vectorField,order)
  {}

  void setOrder(int order){
    solver.setOrder(order);
  }
  void setAbsoluteTolerance(double tol){
    solver.setAbsoluteTolerance(tol);
  }
  void setRelativeTolerance(double tol){
    solver.setRelativeTolerance(tol);
  }
};

template<class TM, class VectorFieldType = typename TM::VectorFieldType>
struct TMap : public BaseTPMap<typename TM::Solver,VectorFieldType>{
  TM tm;
  using BaseTPMap<typename TM::Solver,VectorFieldType>::solver;
  using BaseTPMap<typename TM::Solver,VectorFieldType>::vectorField;

  TMap(const typename TM::VectorFieldType& vf, int order)
    : BaseTPMap<typename TM::Solver,VectorFieldType>(vf,order), tm(solver)
  {}

};

template<class PM, class Section, class VectorFieldType = typename PM::VectorFieldType>
struct PMap : public BaseTPMap<typename PM::Solver,VectorFieldType> {
  using BaseTPMap<typename PM::Solver,VectorFieldType>::solver;
  using BaseTPMap<typename PM::Solver,VectorFieldType>::vectorField;
  Section section;
  PM pm;

  PMap(const VectorFieldType& vf, int order, Section s, typename PM::CrossingDirection dir)
    : BaseTPMap<typename PM::Solver,VectorFieldType>(vf,order),
      section(s),
      pm(solver,section,dir)
  {}
};

template<class PM, class Section>
struct PMapFactory : public SolverFactory< PMap<PM,Section> >{
  typedef PMap<PM,Section> Solver;

  typename PM::VectorFieldType vectorField;
  int order;
  Section section;
  typename PM::CrossingDirection dir;

  PMapFactory(std::string vf, int order, int degree, Section s, typename PM::CrossingDirection dir)
    : vectorField(vf,degree),
      order(order),
      section(s),
      dir(dir)
  {}

  template<class F>
  PMapFactory(F vf, int order, int dim, int params, int degree, Section s, typename PM::CrossingDirection dir)
    : vectorField(vf,dim,dim,params,degree),
      order(order),
      section(s),
      dir(dir)
  {}

  Solver* createSolver() const {
    return new Solver(vectorField,order,section,dir);
  }
};

template<class TM>
struct TMapFactory : public SolverFactory< TMap<TM> >{
  typedef TMap<TM> Solver;

  typename TM::VectorFieldType vectorField;
  int order;

  TMapFactory(std::string vf, int order, int degree)
    : vectorField(vf,degree),
      order(order)
  {}

  template<class F>
  TMapFactory(F vf, int order, int dim, int params, int degree)
    : vectorField(vf,dim,dim,params,degree),
      order(order)
  {}

  Solver* createSolver() const {
    return new Solver(vectorField,order);
  }
};

/// @}
}} // namespace capd::threading

#endif
