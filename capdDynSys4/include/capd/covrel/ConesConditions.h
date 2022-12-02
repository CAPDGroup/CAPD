/////////////////////////////////////////////////////////////////////////////
/// @file ConesConditions.h
///
/// @author Tomasz Kapela  @date 2009-08-03
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2016 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_COVREL_CONESCONDITIONS_H_
#define _CAPD_COVREL_CONESCONDITIONS_H_
#include "capd/auxil/logger_deprecated.h"

#define LOGLN(x)  #x << (x) << "\n"

namespace capd {
namespace covrel {
/// @addtogroup covrel 
/// @{
template <typename MatrixType>
bool isPositiveDefined(const MatrixType & A){
  if((A.numberOfColumns()!=2)||(A.numberOfRows()!=2))
    throw "Implemented only for 2x2 matrices ";
  typename MatrixType::ScalarType det = A[0][0]*A[1][1]-A[1][0]*A[0][1];
  return ((A[0][0]>0) && ((A[0][0]*A[1][1]-A[1][0]*A[0][1])>0));
}
/**
 * Function checks cones condition.
 *
 * Definition of cones condition: \f$ \forall x,y\in |N1|  Q2(f_c(x),f_c(y)) > Q1(x,y) \f$,
 * where f_c is map f expressed in coordinates of N1 and N2
 * Cones condition is satisfied if matrix
 * \f$ V = Df_c^T * Q2 * Df_c - Q1 \f$
 * is positive defined.
 *
 * @param f  map R^n -> R^n
 * @param N1 first h-set
 * @param N2 second h-set
 * @param Q1 quadratic for associated with N1
 * @param Q2 quadratic for associated with N2
 */
template <typename Map, typename HSet2D, typename QuadraticForm>
bool checkConesCondition(Map & f, const HSet2D & N1, const HSet2D & N2,
                         const QuadraticForm & Q1, const QuadraticForm & Q2){
   typedef typename Map::VectorType VectorType;
   typedef typename Map::MatrixType MatrixType;
   VectorType N = N1.center()+N1.coordinateSystem()*N1.box();  // support of N1
   MatrixType dFc = N2.get_I_invB()*f[N]*N1.get_I_B();
   MatrixType V = Transpose(dFc)*Q2.getQ()*dFc - Q1.getQ();
   bool result = isPositiveDefined(V);
  // if(!result)
     capd::sbug  << std::fixed << "checking cones conditions:\n"
           //     << " N1 : " << N1.show() << "\n Q1 :" << Q1.getQ()
           //     << "\n N2 : " << N2.show() << "\n Q2 :" << Q2.getQ()
                << LOGLN(dFc) << LOGLN(V) << LOGLN(result)
            //    << std::setprecision(10) << LOGLN(f[N]) << LOGLN(N) << LOGLN(diam(N))
                << "\n";

   return result;
}

/**
 * Function checks cones condition.
 *
 * Definition of cones condition: \f$ \forall x,y\in |N1|  Q2(f_c(x),f_c(y)) > Q1(x,y) \f$,
 * where f_c is map f expressed in coordinates of N1 and N2
 * Cones condition is satisfied if matrix
 * \f$ V = Df_c^T * Q2 * Df_c - Q1 \f$
 * is positive defined.
 *
 * @param f  map R^n -> R^n
 * @param N1 first h-set with connes
 * @param N2 second h-set with connes
 * @param param  parameters of the computation
 */
template <typename Map, typename HSetWithCones>
bool checkConesCondition(
  Map & f, 
  const HSetWithCones & N1, 
  const HSetWithCones & N2,
  const capd::covrel::CheckCoveringRelation2DParameters & param = capd::covrel::CheckCoveringRelation2DParameters(1, 1, 1, 1)
){
   typedef typename Map::VectorType VectorType;
   typedef typename Map::MatrixType MatrixType;
   typedef capd::covrel::GridSet<MatrixType> GridSet;
   GridSet grid(N1.get_x().dimension());
	 N1.gridSet(grid, param.xGrid, param.yGrid); 
   capd::sbug  << std::fixed << "checking cones conditions:\n";
    
   for (typename GridSet::const_iterator i = grid.begin(); i != grid.end(); ++i) {
     VectorType N = *i+ grid.C * grid.r;
   //VectorType N = N1.center()+N1.coordinateSystem()*N1.box();  // support of N1
    MatrixType dFc = N2.get_I_invB()*f[N]*N1.get_I_B();
    MatrixType V = capd::vectalg::transpose(dFc)*N2.Q().getQ()*dFc - N1.Q().getQ();
    bool result = isPositiveDefined(V);
  // if(!result)
     capd::sbug  << LOGLN(N) 
           //     << " N1 : " << N1.show() << "\n Q1 :" << Q1.getQ()
           //     << "\n N2 : " << N2.show() << "\n Q2 :" << Q2.getQ()
                << LOGLN(dFc) << LOGLN(V) << LOGLN(result)
            //    << std::setprecision(10) << LOGLN(f[N]) << LOGLN(N) << LOGLN(diam(N))
                << "\n";
     if(!result)
       return false;
   }
   return true;
}

/**
 *  Checks covering relation \f$ N_1 \Rightarrow N_2 \f$ and corresponding cones condition.
 *
 *  \see checkConesCondition
 */
template <typename C0SetType, typename Map, typename DynSys, typename HSetWithCones>
bool checkCoveringRelationAndConesCondition(
    Map & f,                      ///<  map R^n -> R^n used to check cones condition
    DynSys & ds,                  ///<  dynamical system used to move set in checking covering relation
    const HSetWithCones & N1,     ///<  first h-set with cones
    const HSetWithCones & N2,     ///<  second h-set with cones
    const capd::covrel::CheckCoveringRelation2DParameters & param = capd::covrel::CheckCoveringRelation2DParameters(1, 1, 1, 1)
){
   return (capd::covrel::checkCoveringRelation2D<C0SetType>(ds, N1, N2, param)) and
          (checkConesCondition(f, N1, N2, param));
}

/**
 *  Checks covering relation \f$ N_1  \Rightarrow  N_2 \f$ and corresponding cones condition.
 *
 *  \see checkConesCondition
 */
template <typename C0SetType, typename Map, typename DynSys, typename HSetWithCones>
bool checkCoveringRelationAndConesCondition(
    Map & f,                      ///<  map R^n -> R^n used to check cones condition
    const HSetWithCones & N1,     ///<  first h-set with cones
    const HSetWithCones & N2,     ///<  second h-set with cones
    const capd::covrel::CheckCoveringRelation2DParameters & param = capd::covrel::CheckCoveringRelation2DParameters(1, 1, 1, 1)
){
   capd::dynsys::DynSysMap<Map> ds(f);
   return checkCoveringRelationAndConesCondition(f, ds, N1, N2, param);
}
/// @}
}} // end of namespace capd::covrel
#undef LOGLN
#endif /* CONESCONDITIONS_H_ */
