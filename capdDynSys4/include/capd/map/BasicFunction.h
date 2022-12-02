/// @addtogroup map
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file BasicFunction.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_MAP_BASICFUNCTION_H_
#define _CAPD_MAP_BASICFUNCTION_H_

#include <string>
#include <map>
#include <vector>
#include "capd/autodiff/DagIndexer.h"
#include "capd/map/Parser.h"
#include "capd/autodiff/NodeType.h"
#include "capd/autodiff/eval.hpp"

namespace capd{
namespace map{

/**
 * This class is a basic for further protected inheritance to classes Function, Map  ...
 * It provides constructors and methods common for all inherited classes. In particular
 * <ul>
 * <li>parsing expressions from a string</li>
 * <li>parsing expressions from routines and computing DAGs</li>
 * <li>memory management</li>
 * <li>setting parameters in expressions</li>
 * <li>evaluation of DAG</li>
 * </ul>
 *\attention The class is for inheritance only and is mainly written to avoid code redundancy. It is for developers of the library only.
 * Users of the CAPD library should not use this class directly in own programs.
 */
template<typename Scalar>
class BasicFunction
{
public:
  typedef Scalar ScalarType;
  typedef capd::autodiff::DagIndexer<ScalarType> DAG;
  typedef typename DAG::size_type size_type;
  typedef capd::autodiff::Node NodeType;

  size_type dimension() const {return m_indexOfFirstParam;}   ///< returns number of main variables (does not count time and parameters)
  size_type imageDimension() const {return m_pos.size(); }    ///< returns number of functions (dimension of counterdomain)

  void setOrder(size_type);                                   ///< allocates necessary memory for ODE solver of order 'order'
  size_type getOrder() const { return m_dag.getOrder(); }     ///< returns current order of expansion with respect to the time variable

  void setParameter(size_type d, const Scalar& value);                  ///< sets new value of a parameter of index d.
  void setParameter(const std::string &name, const ScalarType& value); ///< sets new value of a parameter. If parameter not found there is no effect.
  void setParameters(const Scalar* values, size_type d);                ///< simultaneously sets values of many parameters. It is assumed that given vector contains values of subsequent parameters.

  ScalarType getParameter(size_type d) const;
  ScalarType getParameter(const std::string& name) const;
  void setCurrentTime(const ScalarType& a_time) const;  ///< sets actual value of variable that represents time in an ODE.
  const ScalarType& getCurrentTime() const;             ///< returns actual value of variable that represents time in an ODE
  void differentiateTime() const;                       ///< sets first derivative of time with respect to time equal to 1

  /**
   * The iterator range [b,e) should contain a range of Multiinideces the user requires to compute.
   * The method automatically adds all the depending partial derivatives to this collection and defines
   * a mask for computation of partial derivtives.
   * @param [b,e) iterator range which contains collection of multiindices
   * @warning The method causes undefined behavior if a multiindex in the collection exceeds limits of the map (like dimension, maximal allowed degree).
   */
  template<class Iterator>
  void setMask(Iterator b, Iterator e){
    using namespace capd::autodiff;
    this->m_dag.setMask(b,e);
    this->differentiateTime();
  }
  /**
   * Adds new multiindex (along with dependencies) to the \b existing mask.
   * @param mi multiindex to be added to the mask
   * @warning causes undefined behavior if the mask has not been set before call to this method.
   */
  void addMultiindexToMask(const capd::vectalg::Multiindex& mi){
    this->m_dag.addMultiindexToMask(mi);
    this->differentiateTime();
  }
  /**
   * Resets the mask of derivatives. In consequence, full jet of derivatives will be computed after call to any method that computes derivative, hessian or jet.
   */
  void resetMask(){
    this->m_dag.resetMask();
  }

  const bool* getMask() const {
    return this->m_dag.getMask();
  }

  bool getMask(size_type i) const {
    return this->m_dag.getMask(i);
  }
  bool getMask(size_type i, size_type j) const {
    return this->m_dag.getMask(i,j);
  }
protected:

  BasicFunction();  ///< creates an univariate constant function f(x)=0
  BasicFunction(const std::string& s);   ///< parses expression from given string.
  BasicFunction(const char* s);          ///< parses expression from given string.
  BasicFunction(const BasicFunction& f); ///< copying constructor

  /** parses expression from given routine. Allocates memory for jet propagation of degree 'degree'
   * @param f - routine that defines function, its signature is: void (NodeType time, NodeType in[], int dimIn, NodeType out[], int dimOut, NodeType param[], int noParam)
   * @param dimIn - number of input variables (dimension of domain)
   * @param dimOut - number of output variables (dimension of codomain)
   * @param noParam - number of parameters
   */
  template<typename Function>
  BasicFunction(Function f, int dimIn, int dimOut, int noParam){
    this->reset(f,dimIn,dimOut,noParam);
  }

  /** parses expression from given routine. Allocates memory for jet propagation of degree 'degree'
   * @param f - routine that defines function, its signature is: void (NodeType time, NodeType in[], int dimIn, NodeType out[], int dimOut, NodeType param[], int noParam)
   * @param dimIn - number of input variables (dimension of domain)
   * @param dimOut - number of output variables (dimension of codomain)
   * @param noParam - number of parameters
   */
  template<typename Function>
  void reset(Function f, int dimIn, int dimOut, int noParam);

  /// virtual destructor - releases memory allocated for DAG representing function, map or vector field
  virtual ~BasicFunction() { this->deleteNodes(); }

  void operator=(const std::string&);       ///< substitution from string. Parses expression and compute new DAG.
  void operator=(const char* );             ///< substitution from string. Parses expression and compute new DAG.
  void operator=(const BasicFunction&);     ///< assignment from another object

  void createDefault();                     ///< used to create a default object without expression. Creates a univariate constant function f(x)=0
  void createFromText(std::string s);       ///< parses expression from given string. Allocates memory for ODE solver of order 'order'
  void copyObject(const BasicFunction& f);  ///< an auxiliary function used in copying constructor and assignment operator
  void clean();                             ///< resets all allocated data
  void realloc(size_type degree);           ///< reallocates memory for computation of derivatives up to order 'degree'

  template<class V>
  void setArgument(const V& v) const;
  void applyC1Mask() const;

  void evalHomogenousPolynomial() const {
    for(unsigned i=0;i<this->m_nodes.size();++i)
      this->m_nodes[i]->evalC0HomogenousPolynomial();
  }

  void evalHomogenousPolynomial(size_type degree, size_type coeffNo=0) const {
    if(this->getMask()==0){
      for(unsigned i=0;i<this->m_nodes.size();++i)
        this->m_nodes[i]->evalHomogenousPolynomial(degree,coeffNo);
    } else {
      for(unsigned i=0;i<this->m_nodes.size();++i)
        this->m_nodes[i]->evalHomogenousPolynomial(degree,coeffNo,this->getMask());
    }
  }

  void eval(size_type coeffNo) const {
    for(unsigned i=0;i<this->m_nodes.size();++i)
      this->m_nodes[i]->evalC0(coeffNo);
  }

  void eval(size_type degree, size_type coeffNo) const {
    if(!this->getMask()){
      for(unsigned i=0;i<this->m_nodes.size();++i)
        this->m_nodes[i]->eval(degree,coeffNo);
    } else {
      for(unsigned i=0;i<this->m_nodes.size();++i)
        this->m_nodes[i]->eval(degree,coeffNo,this->getMask());
    }
  }

  void deleteNodes();

  std::vector<capd::autodiff::AbstractNode<ScalarType>*> m_nodes;
  std::vector<capd::autodiff::Node> m_fullGraph;  ///< graph representing the expression
  std::vector<capd::autodiff::MyNode> m_evalPath; ///< reduced graph - only nodes with nontrivial evaluations are left
  std::vector<int> m_pos;                         ///< indices of roots of expressions for each component
  std::vector<std::string> m_var;                 ///< variables, time and parameters
  mutable DAG m_dag;                              ///< data structure that stores all the coefficients. Provides suitable indexing and evaluation of expression.
  size_type m_indexOfFirstParam;                   ///< if equal to m_var.size() then no parameters specified
};


// -----------------------------------------------------------------------------

/// Assignment of new function parsed from a routine.
template<typename Scalar>
template<typename Function>
void BasicFunction<Scalar>::reset(Function f, int dimIn, int dimOut, int noParam)
{
  using namespace autodiff;
  Node* var = new Node[dimIn+noParam+1];
  Node* out = new Node[dimOut];
  int i;
  // set variables
  for(i=0;i<dimIn;++i)
  {
    var[i] = Node(NODE_NULL,NODE_NULL,i,NODE_VAR);
    this->m_fullGraph.push_back(var[i]);
  }
  // set time
  var[i].isTimeDependentOnly = true;
  var[i].isConst = false;
  var[i].op = NODE_TIME;
  var[i].result = i;
  this->m_fullGraph.push_back(var[i]);

  this->m_indexOfFirstParam = dimIn;
  // set parameters
  for(;i<dimIn+noParam+1;++i)
  {
    var[i].isTimeDependentOnly = true;
    var[i].isConst = true;
    var[i].op = NODE_PARAM;
    var[i].result = i;
    this->m_fullGraph.push_back(var[i]);
  }

  // eval expression
  Node::dag = &this->m_fullGraph;
  f(var[dimIn],var,dimIn,out,dimOut,var+dimIn+1,noParam);
  for(i=0;i<dimOut;++i)
    this->m_pos.push_back(out[i].result);

  optimizeDAG(this->m_fullGraph,this->m_pos);

  for(i=0;i<(int)this->m_fullGraph.size();++i)
    if( this->m_fullGraph[i].op!=NODE_NULL and
        this->m_fullGraph[i].op!=NODE_PARAM and
        this->m_fullGraph[i].op!=NODE_CONST and
        this->m_fullGraph[i].op!=NODE_TIME and
        this->m_fullGraph[i].op!=NODE_COS and
        this->m_fullGraph[i].op!=NODE_VAR
     ) this->m_evalPath.push_back(MyNode(this->m_fullGraph[i]));
  Int4ToAbstractNode(this->m_evalPath,this->m_nodes,this->m_dag);

  delete[] var;
  delete[] out;
  this->m_var.resize(dimIn+noParam+1);
}

// ------------------------------------------------------------------------------

template<typename Scalar>
template<class V>
void BasicFunction<Scalar>::setArgument(const V& v) const
{
  using namespace capd::autodiff;
  for(size_type i=0;i<this->dimension();++i)
  {
    this->m_dag(VarNo(i),CoeffNo(0)) = v[i];
    this->m_dag(VarNo(i),DerNo(i),CoeffNo(0)) = TypeTraits<ScalarType>::one();
    for(size_type j=i+1;j<this->dimension();++j){
      this->m_dag(VarNo(i),DerNo(j),CoeffNo(0)) = TypeTraits<ScalarType>::zero();
      this->m_dag(VarNo(j),DerNo(i),CoeffNo(0)) = TypeTraits<ScalarType>::zero();
    }
  }
}
}} // the end of the namespace capd::map

#endif // _CAPD_MAP_BASICFUNCTION_H_

/// @}
