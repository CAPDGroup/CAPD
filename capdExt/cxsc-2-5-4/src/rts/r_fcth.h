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

/* CVS $Id: r_fcth.h,v 1.21 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_fcth.h                              */
/*                                                              */
/*      Description     : function prototypes                   */
/*                                                              */
/*                   include files for hardware                 */
/*                   r_cnsu =b=                                 */
/****************************************************************/

#ifdef IEEE_HARDWARE
/* need sigcontext at r_ftrp */
#if SUN4_OS4_C+SUN4_GNU_C
#include <signal.h>
#include <floatingpoint.h>
#include <sys/ieeefp.h>
#endif /* sun */
#if SUN4_CPP_C
#include <signal.h>
#include <floatingpoint.h>
#endif /* sun */
#if SUN4_OS5_GNU_C
#include <signal.h>
#include <siginfo.h>
#include <ucontext.h>
#endif /* solaris */
#if IBM_LINUX_C+IBM_EMX_C
#include <signal.h>
#endif /* IBM_LINUX_C */
#if IBM_RS6000_C
#include <signal.h>
#include <fpxcp.h>
#include <fptrap.h>
#endif /* IBM_RS6000_C */
#endif

_PROTOTYPE( a_intv  i_acos,(a_intv a));
_PROTOTYPE( a_intv  i_acot,(a_intv a));
_PROTOTYPE( a_intv  i_acsh,(a_intv a));
_PROTOTYPE( a_intv  i_acth,(a_intv a));
_PROTOTYPE( a_intv  i_asin,(a_intv a));
_PROTOTYPE( a_intv  i_asnh,(a_intv a));
_PROTOTYPE( a_intv  i_atan,(a_intv a));
_PROTOTYPE( a_intv  i_atnh,(a_intv a));
_PROTOTYPE( a_intv  i_cnst,(a_char *str));
_PROTOTYPE( a_intv  i_cns2,(a_char *stl,a_char *str));
_PROTOTYPE( a_intv  i_cos,(a_intv a));
_PROTOTYPE( a_intv  i_cosh,(a_intv a));
_PROTOTYPE( a_intv  i_cot,(a_intv a));
_PROTOTYPE( a_intv  i_coth,(a_intv a));
_PROTOTYPE( a_intv  i_ep10,(a_intv a));
_PROTOTYPE( a_intv  i_exp,(a_intv a));
_PROTOTYPE( a_intv  i_exp2,(a_intv a));
_PROTOTYPE( a_intv  i_lg10,(a_intv a));
_PROTOTYPE( a_intv  i_log,(a_intv a));
_PROTOTYPE( a_intv  i_loga,(a_intv a,a_real base));
_PROTOTYPE( a_intv  i_log2,(a_intv a));
_PROTOTYPE( a_intv  i_pow,(a_intv a,a_intv x));
_PROTOTYPE( void    i_read,(FILE *device,a_intv *s,int i));
_PROTOTYPE( a_intv  i_sin,(a_intv a));
_PROTOTYPE( a_intv  i_sinh,(a_intv a));
_PROTOTYPE( a_intv  i_sqrt,(a_intv a));
_PROTOTYPE( a_intv  i_tan,(a_intv a));
_PROTOTYPE( a_intv  i_tanh,(a_intv a));
_PROTOTYPE( a_real  r_abs ,(a_real a));
_PROTOTYPE( a_real  r_acos,(a_real a));
_PROTOTYPE( a_real  r_acot,(a_real a));
_PROTOTYPE( a_real  r_acsh,(a_real a));
_PROTOTYPE( a_real  r_acth,(a_real a));
_PROTOTYPE( a_real  r_addd,(a_real a,a_real b));
_PROTOTYPE( a_real  r_addn,(a_real a,a_real b));
_PROTOTYPE( a_real  r_addu,(a_real a,a_real b));
_PROTOTYPE( a_real  r_asin,(a_real a));
_PROTOTYPE( a_real  r_asnh,(a_real a));
_PROTOTYPE( void    r_ass ,(a_real *x,a_real y));
_PROTOTYPE( a_real  r_atan,(a_real a));
_PROTOTYPE( a_real  r_atn2,(a_real a,a_real b));
_PROTOTYPE( a_real  r_atnh,(a_real a));
_PROTOTYPE( a_real  r_aval,(a_char s[],a_intg n,a_intg rnd));
_PROTOTYPE( a_real  r_ceil,(a_real r));
_PROTOTYPE( a_intg  r_clss,(a_real x));
_PROTOTYPE( a_real  r_cnsd,(char *str));
_PROTOTYPE( a_real  r_cnst,(char *str));
_PROTOTYPE( a_real  r_cnsu,(char *str));
_PROTOTYPE( a_real  r_comp,(a_real m,a_intg e));
_PROTOTYPE( void    r_conv,(a_char *str,a_real *s,a_intg rnd,a_char **next));
_PROTOTYPE( a_real  r_cos ,(a_real a));
_PROTOTYPE( a_real  r_cosh,(a_real a));
_PROTOTYPE( a_real  r_cot ,(a_real a));
_PROTOTYPE( a_real  r_coth,(a_real a));
_PROTOTYPE( a_real  r_divd,(a_real a,a_real b));
_PROTOTYPE( a_real  r_divn,(a_real a,a_real b));
_PROTOTYPE( a_real  r_divu,(a_real a,a_real b));
_PROTOTYPE( a_real  r_ep10,(a_real a));
_PROTOTYPE( a_bool  r_eq  ,(a_real a,a_real b));
_PROTOTYPE( a_real  r_exp ,(a_real a));
_PROTOTYPE( a_real  r_exp2,(a_real a));
_PROTOTYPE( a_intg  r_expo,(a_real m));
_PROTOTYPE( void    r_fini,(void));

#ifdef IEEE_HARDWARE
_PROTOTYPE( void    r_lfsr,(void));

#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C
_PROTOTYPE( void    r_ftrp,( int sig, int code, 
			     struct sigcontext *scp , char *addr));
#endif
#if SUN4_OS5_GNU_C
_PROTOTYPE( void    r_lfsx,(int fsr));
_PROTOTYPE( void    r_ftrp,(int sig, siginfo_t *sip, ucontext_t *ucp));
#endif
#if HP_9000_C
_PROTOTYPE( void r_ftrp,(int code));
#endif
#if IBM_LINUX_C+IBM_EMX_C
_PROTOTYPE( void r_ftrp,(int sig) );
#endif /* IBM_LINUX_C */
#if IBM_RS6000_C
_PROTOTYPE( void r_ftrp,(int sig, int code, struct sigcontext *scp));
#endif /* IBM_RS6000_C */
#else /* ieee_hardware */
_PROTOTYPE( void r_ftrp,(void));
#endif
_PROTOTYPE( a_real  r_flor,(a_real r));
_PROTOTYPE( a_real  r_flot,(a_intg i));
_PROTOTYPE( a_real  r_frac,(a_real r));
_PROTOTYPE( a_bool  r_ge  ,(a_real a,a_real b));
_PROTOTYPE( a_bool  r_gt  ,(a_real a,a_real b));
_PROTOTYPE( a_bool  r_le  ,(a_real a,a_real b));
_PROTOTYPE( a_bool  r_lt  ,(a_real a,a_real b));
_PROTOTYPE( a_real  r_lg10,(a_real a));
_PROTOTYPE( a_real  r_log ,(a_real a));
_PROTOTYPE( a_real  r_log2,(a_real x));
_PROTOTYPE( a_real  r_loga,(a_real x,a_real a));
_PROTOTYPE( a_real  r_mant,(a_real m));
_PROTOTYPE( a_real  r_muld,(a_real a,a_real b));
_PROTOTYPE( a_real  r_muln,(a_real a,a_real b));
_PROTOTYPE( a_real  r_mulu,(a_real a,a_real b));
_PROTOTYPE( a_bool  r_ne  ,(a_real a,a_real b));
_PROTOTYPE( void    r_outp,(char *buffer,a_real s,a_intg TotalWidth,
               a_intg FracDigits,a_intg rnd,a_intg *length));
_PROTOTYPE( a_intg  r_pcmp,(a_real a,a_real b,a_real c,a_real d));
_PROTOTYPE( a_real  r_pow ,(a_real a,a_real b));
_PROTOTYPE( a_real  r_pred,(a_real a));
_PROTOTYPE( void    r_rdcr,(FILE *device,a_real *s,a_intg rnd,int i));
#ifdef DEC_ARITH
_PROTOTYPE( void    r_read,(f_text *device,a_real *s,a_intg rnd));
#else
_PROTOTYPE( void    r_read,(FILE *device,a_real *s,a_intg rnd,int i));
#endif
_PROTOTYPE( a_intg  r_rond,(a_real a));
_PROTOTYPE( a_real  r_rval,(s_trng s,a_intg rnd));
_PROTOTYPE( a_intg  r_sign,(a_real a));
_PROTOTYPE( a_real  r_sin ,(a_real a));
_PROTOTYPE( a_real  r_sinh,(a_real a));
_PROTOTYPE( a_real  r_sqr_,(a_real a));
_PROTOTYPE( a_real  r_sqrt,(a_real a));
_PROTOTYPE( a_real  r_subd,(a_real a,a_real b));
_PROTOTYPE( a_real  r_subn,(a_real a,a_real b));
_PROTOTYPE( a_real  r_subu,(a_real a,a_real b));
_PROTOTYPE( a_real  r_succ,(a_real a));
_PROTOTYPE( a_real  r_sval,(s_trng s,a_intg rnd,s_trng *r));
_PROTOTYPE( a_real  r_tan ,(a_real a));
_PROTOTYPE( a_real  r_tanh,(a_real a));
_PROTOTYPE( a_intg  r_trun,(a_real a));
_PROTOTYPE( a_real  r_umin,(a_real a));
_PROTOTYPE( a_real  r_valu,(a_intg x));
#ifdef DEC_ARITH
_PROTOTYPE( void    r_writ,(f_text *device,a_real s,a_intg TotalWidth,
               a_intg FracDigits,a_intg rnd));
#else
_PROTOTYPE( void    r_writ,(FILE *device,a_real s,a_intg TotalWidth,
               a_intg FracDigits,a_intg rnd));
#endif





