#ifndef CAPD_DYNSYS_BASIC_SOLVER_MASK_TEST_DATA_H
#define CAPD_DYNSYS_BASIC_SOLVER_MASK_TEST_DATA_H

#include <boost/test/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/vectalg/Matrix.hpp"
#include "capd/diffAlgebra/Jet.hpp"
#include "capd/diffAlgebra/Hessian.hpp"

using namespace capd;
typedef capd::vectalg::Matrix<int,0,0> BMatrix;
typedef capd::diffAlgebra::Hessian<bool,0,0> BHessian;
typedef capd::diffAlgebra::Jet<BMatrix,0> BJet;

template<class M, class O>
void check(M m, O b, O e){
  for(;b!=e;++m,++b){
    if(not *m)
      BOOST_CHECK_EQUAL(*b,0.);
  }
}

// compare numerical solution to exact solution - C^1 part
inline void check(const BMatrix& mask, const DMatrix& der){
  check(mask.begin(),der.begin(),der.end());
}

// compare numerical solution to exact solution - C^2 part
inline void check(const BHessian& mask, const DHessian& h){
  check(mask.begin(),h.begin(),h.end());
}


#endif /* CAPD_DYNSYS_BASIC_SOLVER_MASK_TEST_DATA_H */

