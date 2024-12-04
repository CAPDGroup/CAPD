/////////////////////////////////////////////////////////////////////////////
/// @file eval.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_HPP_
#define _CAPD_AUTODIFF_EVAL_HPP_

#define CAPD_MAKE_NODE(NodeName,ClassName) case NodeName : { p = new ClassName##Node<T>;  break; }

#include <stdexcept>
#include <sstream>

#include "capd/autodiff/DagIndexer.h"
#include "capd/autodiff/EvalAcos.h"
#include "capd/autodiff/EvalAdd.h"
#include "capd/autodiff/EvalAsin.h"
#include "capd/autodiff/EvalAtan.h"
#include "capd/autodiff/EvalDiv.h"
#include "capd/autodiff/EvalExp.h"
#include "capd/autodiff/EvalLog.h"
#include "capd/autodiff/EvalMul.h"
#include "capd/autodiff/EvalNaturalPow.h"
#include "capd/autodiff/EvalHalfIntPow.h"
#include "capd/autodiff/EvalNegIntPow.h"
#include "capd/autodiff/EvalOneMinusSqr.h"
#include "capd/autodiff/EvalPow.h"
#include "capd/autodiff/EvalQuarticPow.h"
#include "capd/autodiff/EvalSinCos.h"
#include "capd/autodiff/EvalSqr.h"
#include "capd/autodiff/EvalSqrt.h"
#include "capd/autodiff/EvalSub.h"
#include "capd/autodiff/EvalThirdPow.h"
#include "capd/autodiff/EvalUnaryMinus.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{
template<class T>
void Int4ToAbstractNode(const std::vector<MyNode>& node, std::vector<AbstractNode<T>* >& out, capd::autodiff::DagIndexer<T>& dag)
{
  const JetSize jetSize = dag.timeJetSize();
  T* data = dag.coefficients();
  out.resize(node.size());
  for(unsigned i=0;i<node.size();++i)
  {
    AbstractNode<T>* p;
    T* left = data + jetSize*node[i].left;
    T* right = data + jetSize*node[i].right;
    T* result = data + jetSize*node[i].result;

    switch(node[i].op)
    {
    CAPD_MAKE_NODE(NODE_ADD,Add);
    CAPD_MAKE_NODE(NODE_CONST_PLUS_VAR,ConstPlusVar);
    CAPD_MAKE_NODE(NODE_CONST_PLUS_CONST,ConstPlusConst);
    CAPD_MAKE_NODE(NODE_CONST_PLUS_TIME,ConstPlusTime);
    CAPD_MAKE_NODE(NODE_CONST_PLUS_FUNTIME,ConstPlusFunTime);
    CAPD_MAKE_NODE(NODE_TIME_PLUS_VAR,TimePlusVar);
    CAPD_MAKE_NODE(NODE_TIME_PLUS_FUNTIME,TimePlusFunTime);
    CAPD_MAKE_NODE(NODE_FUNTIME_PLUS_VAR,FunTimePlusVar);
    CAPD_MAKE_NODE(NODE_FUNTIME_PLUS_FUNTIME,FunTimePlusFunTime);

    CAPD_MAKE_NODE(NODE_SUB,Sub);
    CAPD_MAKE_NODE(NODE_CONST_MINUS_CONST,ConstMinusConst);
    CAPD_MAKE_NODE(NODE_CONST_MINUS_VAR,ConstMinusVar);
    CAPD_MAKE_NODE(NODE_CONST_MINUS_TIME,ConstMinusTime);
    CAPD_MAKE_NODE(NODE_CONST_MINUS_FUNTIME,ConstMinusFunTime);
    CAPD_MAKE_NODE(NODE_TIME_MINUS_CONST,TimeMinusConst);
    CAPD_MAKE_NODE(NODE_TIME_MINUS_FUNTIME,TimeMinusFunTime);
    CAPD_MAKE_NODE(NODE_TIME_MINUS_VAR,TimeMinusVar);
    CAPD_MAKE_NODE(NODE_FUNTIME_MINUS_CONST,FunTimeMinusConst);
    CAPD_MAKE_NODE(NODE_FUNTIME_MINUS_TIME,FunTimeMinusTime);
    CAPD_MAKE_NODE(NODE_FUNTIME_MINUS_FUNTIME,FunTimeMinusFunTime);
    CAPD_MAKE_NODE(NODE_FUNTIME_MINUS_VAR,FunTimeMinusVar);
    CAPD_MAKE_NODE(NODE_VAR_MINUS_CONST,VarMinusConst);
    CAPD_MAKE_NODE(NODE_VAR_MINUS_TIME,VarMinusTime);
    CAPD_MAKE_NODE(NODE_VAR_MINUS_FUNTIME,VarMinusFunTime);

    CAPD_MAKE_NODE(NODE_UNARY_MINUS,UnaryMinus);
    CAPD_MAKE_NODE(NODE_UNARY_MINUS_CONST,UnaryMinusConst);
    CAPD_MAKE_NODE(NODE_UNARY_MINUS_TIME,UnaryMinusTime);
    CAPD_MAKE_NODE(NODE_UNARY_MINUS_FUNTIME,UnaryMinusFunTime);

    CAPD_MAKE_NODE(NODE_MUL,Mul);
    CAPD_MAKE_NODE(NODE_MUL_CONST_BY_VAR,MulConstByVar);
    CAPD_MAKE_NODE(NODE_MUL_CONST_BY_CONST,MulConstByConst);
    CAPD_MAKE_NODE(NODE_MUL_CONST_BY_TIME,MulConstByTime);
    CAPD_MAKE_NODE(NODE_MUL_CONST_BY_FUNTIME,MulConstByFunTime);
    CAPD_MAKE_NODE(NODE_MUL_TIME_BY_VAR,MulTimeByVar);
    CAPD_MAKE_NODE(NODE_MUL_TIME_BY_FUNTIME,MulTimeByFunTime);
    CAPD_MAKE_NODE(NODE_MUL_FUNTIME_BY_VAR,MulFunTimeByVar);
    CAPD_MAKE_NODE(NODE_MUL_FUNTIME_BY_FUNTIME,MulFunTimeByFunTime);

    CAPD_MAKE_NODE(NODE_DIV,Div);
    CAPD_MAKE_NODE(NODE_DIV_VAR_BY_CONST,DivVarByConst);
    CAPD_MAKE_NODE(NODE_DIV_VAR_BY_TIME,DivVarByTime);
    CAPD_MAKE_NODE(NODE_DIV_VAR_BY_FUNTIME,DivVarByFunTime);
    CAPD_MAKE_NODE(NODE_DIV_TIME_BY_CONST,DivTimeByConst);
    CAPD_MAKE_NODE(NODE_DIV_FUNTIME_BY_CONST,DivFunTimeByConst);
    CAPD_MAKE_NODE(NODE_DIV_FUNTIME_BY_TIME,DivFunTimeByTime);
    CAPD_MAKE_NODE(NODE_DIV_FUNTIME_BY_FUNTIME,DivFunTimeByFunTime);
    CAPD_MAKE_NODE(NODE_DIV_CONST_BY_CONST,DivConstByConst);

    CAPD_MAKE_NODE(NODE_SQR,Sqr);
    CAPD_MAKE_NODE(NODE_SQR_CONST,SqrConst);
    CAPD_MAKE_NODE(NODE_SQR_TIME,SqrTime);
    CAPD_MAKE_NODE(NODE_SQR_FUNTIME,SqrFunTime);

    CAPD_MAKE_NODE(NODE_CUBE,Cube);
    CAPD_MAKE_NODE(NODE_CUBE_CONST,CubeConst);
    CAPD_MAKE_NODE(NODE_CUBE_TIME,CubeTime);
    CAPD_MAKE_NODE(NODE_CUBE_FUNTIME,CubeFunTime);

    CAPD_MAKE_NODE(NODE_QUARTIC,Quartic);
    CAPD_MAKE_NODE(NODE_QUARTIC_CONST,QuarticConst);
    CAPD_MAKE_NODE(NODE_QUARTIC_TIME,QuarticTime);
    CAPD_MAKE_NODE(NODE_QUARTIC_FUNTIME,QuarticFunTime);

    CAPD_MAKE_NODE(NODE_SQRT,Sqrt);
    CAPD_MAKE_NODE(NODE_SQRT_CONST,SqrtConst);
    CAPD_MAKE_NODE(NODE_SQRT_TIME,SqrtTime);
    CAPD_MAKE_NODE(NODE_SQRT_FUNTIME,SqrtFunTime);

    CAPD_MAKE_NODE(NODE_POW,Pow);
    CAPD_MAKE_NODE(NODE_POW_CONST,PowConst);
    CAPD_MAKE_NODE(NODE_POW_TIME,PowTime);
    CAPD_MAKE_NODE(NODE_POW_FUNTIME,PowFunTime);

    CAPD_MAKE_NODE(NODE_NATURAL_POW,NaturalPow);
    CAPD_MAKE_NODE(NODE_NATURAL_POW_CONST,NaturalPowConst);
    CAPD_MAKE_NODE(NODE_NATURAL_POW_TIME,NaturalPowTime);
    CAPD_MAKE_NODE(NODE_NATURAL_POW_FUNTIME,NaturalPowFunTime);

    CAPD_MAKE_NODE(NODE_NEG_INT_POW,NegIntPow);
    CAPD_MAKE_NODE(NODE_NEG_INT_POW_CONST,NegIntPowConst);
    CAPD_MAKE_NODE(NODE_NEG_INT_POW_TIME,NegIntPowTime);
    CAPD_MAKE_NODE(NODE_NEG_INT_POW_FUNTIME,NegIntPowFunTime);

    CAPD_MAKE_NODE(NODE_HALF_INT_POW,HalfIntPow);
    CAPD_MAKE_NODE(NODE_HALF_INT_POW_CONST,HalfIntPowConst);
    CAPD_MAKE_NODE(NODE_HALF_INT_POW_TIME,HalfIntPowTime);
    CAPD_MAKE_NODE(NODE_HALF_INT_POW_FUNTIME,HalfIntPowFunTime);

    CAPD_MAKE_NODE(NODE_EXP,Exp);
    CAPD_MAKE_NODE(NODE_EXP_CONST,ExpConst);
    CAPD_MAKE_NODE(NODE_EXP_TIME,ExpTime);
    CAPD_MAKE_NODE(NODE_EXP_FUNTIME,ExpFunTime);
    CAPD_MAKE_NODE(NODE_LOG,Log);
    CAPD_MAKE_NODE(NODE_LOG_CONST,LogConst);
    CAPD_MAKE_NODE(NODE_LOG_TIME,LogTime);
    CAPD_MAKE_NODE(NODE_LOG_FUNTIME,LogFunTime);

    CAPD_MAKE_NODE(NODE_SIN,Sin);
    CAPD_MAKE_NODE(NODE_SIN_CONST,SinConst);
    CAPD_MAKE_NODE(NODE_SIN_TIME,SinTime);
    CAPD_MAKE_NODE(NODE_SIN_FUNTIME,SinFunTime);

    CAPD_MAKE_NODE(NODE_ONE_MINUS_SQR,OneMinusSqr);
    CAPD_MAKE_NODE(NODE_ONE_MINUS_SQR_CONST,OneMinusSqrConst);
    CAPD_MAKE_NODE(NODE_ONE_MINUS_SQR_TIME,OneMinusSqrTime);
    CAPD_MAKE_NODE(NODE_ONE_MINUS_SQR_FUNTIME,OneMinusSqrFunTime);

    CAPD_MAKE_NODE(NODE_ATAN,Atan);
    CAPD_MAKE_NODE(NODE_ATAN_CONST,AtanConst);
    CAPD_MAKE_NODE(NODE_ATAN_TIME,AtanTime);
    CAPD_MAKE_NODE(NODE_ATAN_FUNTIME,AtanFunTime);

    CAPD_MAKE_NODE(NODE_ASIN,Asin);
    CAPD_MAKE_NODE(NODE_ASIN_CONST,AsinConst);
    CAPD_MAKE_NODE(NODE_ASIN_TIME,AsinTime);
    CAPD_MAKE_NODE(NODE_ASIN_FUNTIME,AsinFunTime);

    CAPD_MAKE_NODE(NODE_ACOS,Acos);
    CAPD_MAKE_NODE(NODE_ACOS_CONST,AcosConst);
    CAPD_MAKE_NODE(NODE_ACOS_TIME,AcosTime);
    CAPD_MAKE_NODE(NODE_ACOS_FUNTIME,AcosFunTime);

    default:
      std::ostringstream out;
      out << "Implementation error! Node no. " << node[i].op << " is missed in switch-case block in function evaluating DAG.\nPlease report this bug to developers!\n ";
      throw std::logic_error(out.str());
    }
    p->left = left;
    p->right = right;
    p->result = result;
    p->setDag(&dag);
    out[i] = p;
  }
}
///@}
}} // namespace capd::autodiff

#endif
