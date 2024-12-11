/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.5.4)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2014 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* CVS $Id: fi_lib.hpp,v 1.19 2014/01/30 17:23:54 cxsc Exp $ */

#include "interval.hpp"
#include "real.hpp"
#include "float.h"
#include <iostream>
#include <cstdlib>

#define erf_kettenbruch 1

namespace fi_lib {
  using namespace cxsc;
  //using cxsc::interval;
  //using cxsc::real;
  
  bool NANTEST(real);
  real q_rtrg(real x, long int k);
  
  interval j_erf(interval);
  interval j_erfc(interval);
  interval j_acos(interval);
  interval j_acsh(interval);
  interval j_acot(interval);
  interval j_acth(interval);
  interval j_asin(interval);
  interval j_asnh(interval);
  interval j_atan(interval);
  interval j_atnh(interval);
  interval j_cos(interval);
  interval j_cosh(interval);
  interval j_cot(interval);
  interval j_coth(interval);
  interval j_exp(interval);
  interval j_ex10(interval);
  interval j_exp2(interval);
  interval j_expm(interval);
  interval j_log(interval);
  interval j_lg10(interval);
  interval j_lg1p(interval);
  interval j_log2(interval);
  interval j_sin(interval);
  interval j_sinh(interval);
  interval j_sqr(interval);
  interval j_sqrt(interval);
  interval j_tan(interval);
  interval j_tanh(interval);
  
  // internal
  real erf_intv(const real&);
  real erfc_intv(const real&);

  
  // real functions
  
  real q_expx2(real);
  real q_erf(real);
  real q_erfc(real);
  real q_acos(real);
  real q_acsh(real);
  real q_acth(real);
  real q_pred(real);  
  real q_asin(real);
  real q_succ(real);
  real q_asnh(real); 
  real q_atan(real);
  real q_atnh(real);
  real q_cos1(real, long int);
  real q_ex10(real);
  real q_exp2(real);
  real q_lg10(real);
  real q_log2(real);
  real q_sin(real x);
  real q_sin1(real, long int);
  real q_sqr(real);
  real q_abortnan(int, real*, int);
  real q_abortr1(int, real*, int);
  interval q_abortr2(int, real*, real*, int);
  real q_abortdivd(int, real*);
  interval q_abortdivi(int, real*, real*);
  real q_atn1(real);
  real q_p1l1(int, real, real, real);
  real q_p2l1(real);
  real q_log1(real);
  real q_l1p1(real);
  real q_acot(real);
  real q_cos(real);
  real q_cot(real);
  real q_tan(real);
  real q_sqrt(real);
  real q_expm(real);
  real q_p2ex(real);
  real q_p1ex(real);
  real q_p1lg(int, real, real, real);
  real q_p2lg(real);
  real q_log(real);
  real q_lg1p(real);
  real q_exp(real);
  real q_sinh(real);
  real q_cosh(real);
  real q_coth(real);
  real q_tanh(real);   
  real q_expx2(real);
  real q_ep1(real);
  real q_p1e1(real);
  real q_p2e1(real);
  real q_epm1(real);
  real q_rtrg(real, long int);
  real q_r2tr(real, long int);
  real q_cth1(real);
  real scandown();
  real scanup();
  interval scanInterval();
  void printup(real);
  void printdown(real);
  void printInterval(interval);
  real q_abs(real);
  real q_min(real, real);
  real q_max(real, real);
  real q_mid(interval);
  interval eq_ii(interval);
  interval eq_id(real);
  interval add_ii(interval, interval);
  interval add_id(interval, real);
  interval add_di(real, interval);
  interval sub_ii(interval, interval);
  interval sub_id(interval, real);
  interval sub_di(real, interval);
  interval mul_ii(interval, interval);
  interval mul_id(interval, real);
  interval mul_di(real, interval);
  interval div_ii(interval,interval);
  interval div_di(real,interval);
  interval div_id(interval, real);
  int in_di(real,interval);
  int in_ii(interval,interval);
  int ieq_ii(interval,interval);
  int is_ii(interval,interval);
  int ig_ii(interval,interval);
  int ise_ii(interval,interval);
  int ige_ii(interval,interval);
  int dis_ii(interval,interval);
  interval hull(interval, interval);
  interval intsec(interval, interval);
  real q_diam(interval);
  real q_comp(int, real, int);
  real q_cmps(real, int);
  int q_sign(real);
  real q_mant(real);
  real q_mnts(real);
  int q_expo(real);
}

#include "fi_lib_consts.hpp"





