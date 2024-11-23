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

/* CVS $Id: rmath.inl,v 1.31 2014/01/30 17:23:48 cxsc Exp $ */

#include "rtsrmath.h"
#include "xscclass.hpp"
#define CXSC_INCLUDE
#undef LINT_ARGS
#include <fi_lib.hpp>
#undef CXSC_INCLUDE

//#include "fi_lib.hpp"

namespace cxsc{
using namespace fi_lib;
inline real t_std_fct_call(int (fct)(const ExtReal *,ExtReal *),const real &arg)
{
   real erg(arg);
   int rnd=t_grnd();
   t_srnd(NEAR);
   ExtReal a,r;
   t_ltoe((LongReal*)&erg,&a);
   fct(&a,&r);
   t_etol(&r,(LongReal*)&erg);
   t_srnd(rnd);
   return erg;
}

inline real t_std_fct_call(int (fct)(const ExtReal *,const ExtReal *,ExtReal *),const real & arg1, const real & arg2)
{
   real erg(arg1),expo(arg2); // expo semioptimal.. :I
   int rnd=t_grnd();
   t_srnd(NEAR);
   ExtReal a,b,r;
   t_ltoe((LongReal*)&erg,&a);
   t_ltoe((LongReal*)&expo,&b);
   fct(&a,&b,&r);
   t_etol(&r,(LongReal*)&erg);
   t_srnd(rnd);
   return erg;
}

inline real sqr(const real &arg) throw() { return (arg*arg); }
inline real sqrt(const real & arg)       { return q_sqrt(*(double *)&arg); } 
// { return t_std_fct_call(t_sqte,arg); }
inline real sqrt(const real & arg,int n) { return pow(arg,1.0/n); }
// inline real sqrtm1(const real & arg)     { return t_std_fct_call(t_sqme,arg); }
inline real sqrtm1(const real & arg) { 
   real erg(arg);
   int rnd=t_grnd();
   t_srnd(NEAR);
   ExtReal a,r;
   t_ltoe((LongReal*)&erg,&a);
   t_sqme(&a,&r);
   t_etol(&r,(LongReal*)&erg);
   t_srnd(rnd);
   return erg;
}

inline real sin(const real & arg) throw() { return q_sin(*(double*)&arg); } // { return t_std_fct_call(t_sine,arg); }
inline real cos(const real & arg) throw() { return q_cos(*(double*)&arg); } // { return t_std_fct_call(t_cose,arg); }
inline real tan(const real & arg) throw() { return q_tan(*(double*)&arg); } //{ return t_std_fct_call(t_tane,arg); }
inline real cot(const real & arg) throw() { return q_cot(*(double*)&arg); } //{ return t_std_fct_call(t_cote,arg); }

inline real asin(const real & arg)       { return q_asin(*(double*)&arg); } // { return t_std_fct_call(t_asne,arg); }
inline real acos(const real & arg)       { return q_acos(*(double*)&arg); } // { return t_std_fct_call(t_acse,arg); }
inline real atan(const real & arg)       { return q_atan(*(double*)&arg); } // { return t_std_fct_call(t_atne,arg); }
inline real acot(const real & arg)       { return q_acot(*(double*)&arg); } // { return t_std_fct_call(t_acte,arg); }

inline real expm1(const real & arg)throw() { return q_expm(*(double*)&arg); } // { return t_std_fct_call(t_exme,arg); }
inline real lnp1(const real & arg)       { return q_lg1p(*(double*)&arg); } // { return t_std_fct_call(t_lnpe,arg); }

inline real exp(const real & arg) throw() { return q_exp(*(double*)&arg); } // { return t_std_fct_call(t_expe,arg); }
inline real ln(const real & arg)         { return q_log(*(double*)&arg); } // { return t_std_fct_call(t_lnee,arg); }
inline real log2(const real & arg)       { return q_log2(*(double*)&arg); } // { return t_std_fct_call(t_lnee,arg); }
inline real log10(const real & arg)      { return q_lg10(*(double*)&arg); } // { return t_std_fct_call(t_lnee,arg); }

inline real sinh(const real & arg) throw() { return q_sinh(*(double*)&arg); } // { return t_std_fct_call(t_snhe,arg); }
inline real cosh(const real & arg) throw() { return q_cosh(*(double*)&arg); } // { return t_std_fct_call(t_cshe,arg); }
inline real tanh(const real & arg) throw() { return q_tanh(*(double*)&arg); } // { return t_std_fct_call(t_tnhe,arg); }
inline real coth(const real & arg) throw() { return q_coth(*(double*)&arg); } // { return t_std_fct_call(t_cthe,arg); }

inline real asinh(const real & arg) { return q_asnh(*(double*)&arg); } // { return t_std_fct_call(t_ashe,arg); }
inline real acosh(const real & arg) { return q_acsh(*(double*)&arg); } // { return t_std_fct_call(t_ache,arg); }
inline real atanh(const real & arg) { return q_atnh(*(double*)&arg); } // { return t_std_fct_call(t_anhe,arg); }
inline real acoth(const real & arg) { return q_acth(*(double*)&arg); } // { return t_std_fct_call(t_athe,arg); }

/*!
\param arg The value for which to compute the value of the error function
\return The computed result of the error function

In mathematics, the error function (also called the Gauss error function) is a non-elementary function which occurs in probability, statistics and partial differential equations.

When the results of a series of measurements are described by a normal distribution with standard deviation s and expected value 0, then
\f$ \mbox{erf} \left( \frac{a}{\sigma \sqrt{2}} \right) \f$ is the probability that the error of a single measurement lies between \f$ -a \f$ and \f$ +a \f$ .

The error and complementary error functions occur, for example, in solutions of the heat equation when boundary conditions are given by the Heaviside step function.
*/
inline real erf(const real & arg)   { return q_erf(*(double*)&arg);  } // { return t_std_fct_call(t_athe,arg); }
/*!
\param arg The value for which to compute the value of the complementary error function
\return The computed result of the complementary error function

\sa erf(const real & arg)
*/
inline real erfc(const real & arg)  { return q_erfc(*(double*)&arg); } // { return t_std_fct_call(t_athe,arg); }

//inline real pow(const real & arg,const real &expo) { return t_std_fct_call(t_powe,arg,expo); }
inline real pow(const real & arg,const real &expo) { 
   real erg(arg),expohelp(expo); // expo semioptimal.. :I
   int rnd=t_grnd();
   t_srnd(NEAR);
   ExtReal a,b,r;
   t_ltoe((LongReal*)&erg,&a);
   t_ltoe((LongReal*)&expohelp,&b);
   t_powe(&a,&b,&r);
   t_etol(&r,(LongReal*)&erg);
   t_srnd(rnd);
   return erg;
}

inline real power(const real & arg,const int n)    { return pow(arg,n); }

} // namespace cxsc

