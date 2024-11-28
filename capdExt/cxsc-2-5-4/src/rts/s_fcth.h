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

/* CVS $Id: s_fcth.h,v 1.21 2014/01/30 17:24:13 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_fcth.h                              */
/*                                                              */
/*      Description     : function prototypes                   */
/*                                                              */
/*                   new function s_whex, s_rhex  =b=           */
/*                   s_cons prototype                           */
/****************************************************************/

_PROTOTYPE( s_trng s_date, (s_trng fmt));
_PROTOTYPE( s_trng s_whex, (a_real s,a_char mode));
_PROTOTYPE( a_real s_rhex, (s_trng s));

_PROTOTYPE( a_VOID s_add, (s_etof res,s_etof s,s_etof t));
#ifdef LINT_ARGS
#if C_P_7
a_VOID s_cons(s_etof s,a_char *e_argc,...);
#else
a_VOID s_cons();
#endif
#else
a_VOID s_cons();
#endif
_PROTOTYPE( a_bool s_aaeq, (a_char s[],a_intg n,a_char t[],a_intg m));
_PROTOTYPE( a_bool s_aceq, (a_char s[],a_intg n,a_char t));
/* a_bool s_aseq, (a_char a[],a_intg n,s_trng s)); */
_PROTOTYPE( a_bool s_caeq, (a_char s,a_char t[],a_intg m));
_PROTOTYPE( a_bool s_eteq, (s_etof s,s_etof t));
_PROTOTYPE( a_bool s_etge, (s_etof s,s_etof t));
_PROTOTYPE( a_bool s_etgt, (s_etof s,s_etof t));
_PROTOTYPE( a_bool s_etin, (a_intg e,s_etof s));
_PROTOTYPE( a_bool s_etle, (s_etof s,s_etof t));
_PROTOTYPE( a_bool s_etlt, (s_etof s,s_etof t));
_PROTOTYPE( a_bool s_etne, (s_etof s,s_etof t));
_PROTOTYPE( a_bool s_aage, (a_char s[],a_intg n,a_char t[],a_intg m));
/* a_bool s_asge, (a_char a[],a_intg n,s_trng s)); */
_PROTOTYPE( a_bool s_acge, (a_char s[],a_intg n,a_char t));
_PROTOTYPE( a_bool s_cage, (a_char s,a_char t[],a_intg m));
_PROTOTYPE( a_bool s_aagt, (a_char s[],a_intg n,a_char t[],a_intg m));
_PROTOTYPE( a_bool s_acgt, (a_char s[],a_intg n,a_char t));
/* a_bool s_asgt, (a_char a[],a_intg n,s_trng s)); */
_PROTOTYPE( a_bool s_cagt, (a_char s,a_char t[],a_intg m));
_PROTOTYPE( a_VOID s_ins1, (s_etof s,a_intg e));
_PROTOTYPE( a_VOID s_ins2, (s_etof s,a_intg rs,a_intg re));
_PROTOTYPE( a_bool s_cain, (a_char s,a_char t[],a_intg m));
_PROTOTYPE( a_bool s_aain, (a_char s[],a_intg n,a_char t[],a_intg m));
_PROTOTYPE( a_bool s_aale, (a_char s[],a_intg n,a_char t[],a_intg m));
_PROTOTYPE( a_bool s_acle, (a_char s[],a_intg n,a_char t));
/* a_bool s_asle, (a_char a[],a_intg n,s_trng s)); */
_PROTOTYPE( a_bool s_cale, (a_char s,a_char t[],a_intg m));
_PROTOTYPE( a_bool s_aalt, (a_char s[],a_intg n,a_char t[],a_intg m));
_PROTOTYPE( a_bool s_aclt, (a_char s[],a_intg n,a_char t));
/* a_bool s_aslt, (a_char a[],a_intg n,s_trng s)); */
_PROTOTYPE( a_bool s_calt, (a_char s,a_char t[],a_intg m));
_PROTOTYPE( a_VOID s_mul , (s_etof res,s_etof s,s_etof t));
_PROTOTYPE( a_bool s_aane, (a_char s[],a_intg n,a_char t[],a_intg m));
_PROTOTYPE( a_bool s_acne, (a_char s[],a_intg n,a_char t));
/* a_bool s_asne, (a_char a[],a_intg n,s_trng s)); */
_PROTOTYPE( a_bool s_cane, (a_char s,a_char t[],a_intg m));
_PROTOTYPE( a_VOID s_sub , (s_etof res,s_etof s,s_etof t));
_PROTOTYPE( a_VOID s_zero, (s_etof s));
_PROTOTYPE( a_VOID s_aacc, (a_char s[],a_char t[],a_char u[]));
_PROTOTYPE( a_VOID s_accc, (a_char s[],a_char t[],a_char u));
_PROTOTYPE( a_VOID s_cacc, (a_char s[],a_char t,  a_char u[]));

_PROTOTYPE( a_char *s_aimg, (a_char str[],a_intg n,a_real s,a_intg TotalWidth,
               a_intg FracDigits,a_intg rnd));
_PROTOTYPE( void   s_asgn, (s_trng *d,s_trng s));
_PROTOTYPE( a_intg s_alen, (s_trng *s));
_PROTOTYPE( s_trng s_char, (a_char s));
_PROTOTYPE( s_trng s_conc, (s_trng d,s_trng s));
_PROTOTYPE( a_intg s_cpos, (a_char s,s_trng t));
_PROTOTYPE( void   s_free, (s_trng *d));
_PROTOTYPE( void   s_genv, (s_trng name,s_trng *value,a_bool *exists));
_PROTOTYPE( void   s_init, (s_trng *d,size_t size));
_PROTOTYPE( s_trng s_int_, (a_intg s,a_intg TotalWidth));
_PROTOTYPE( a_char *s_inxc, (s_trng s,a_intg i));
_PROTOTYPE( a_intg s_ixch, (a_intg i,size_t length));
_PROTOTYPE( a_intg s_len_, (s_trng d));
_PROTOTYPE( s_trng s_real, (a_real s,a_intg TotalWidth,a_intg FracDigits,a_intg rnd));
_PROTOTYPE( void   s_slen, (s_trng *d,a_intg n));
_PROTOTYPE( a_intg s_spos, (s_trng s,s_trng t));
_PROTOTYPE( s_trng s_stat, (a_char s[],a_intg n));
_PROTOTYPE( s_trng s_suba, (s_trng s,a_intg pos,a_intg end));
_PROTOTYPE( s_trng s_subs, (s_trng s,a_intg pos,a_intg len));
_PROTOTYPE( void   s_temp, (s_trng *d));
_PROTOTYPE( void   s_utmp, (s_trng *d));
_PROTOTYPE( void   s_vlcp, (s_trng *d));
/* a_bool s_saeq, (s_trng s,a_char a[],a_intg n));
   a_bool s_sage, (s_trng s,a_char a[],a_intg n));
   a_bool s_sagt, (s_trng s,a_char a[],a_intg n));
   a_bool s_sale, (s_trng s,a_char a[],a_intg n));
   a_bool s_salt, (s_trng s,a_char a[],a_intg n));
   a_bool s_sane, (s_trng s,a_char a[],a_intg n));
*/
_PROTOTYPE( a_bool s_sseq, (s_trng s,s_trng t));
_PROTOTYPE( a_bool s_ssge, (s_trng s,s_trng t));
_PROTOTYPE( a_bool s_ssgt, (s_trng s,s_trng t));
_PROTOTYPE( a_bool s_ssle, (s_trng s,s_trng t));
_PROTOTYPE( a_bool s_sslt, (s_trng s,s_trng t));
_PROTOTYPE( a_bool s_ssne, (s_trng s,s_trng t));
_PROTOTYPE( a_bool s_sceq, (s_trng s,a_char t));
_PROTOTYPE( a_bool s_scge, (s_trng s,a_char t));
_PROTOTYPE( a_bool s_scgt, (s_trng s,a_char t));
_PROTOTYPE( a_bool s_scle, (s_trng s,a_char t));
_PROTOTYPE( a_bool s_sclt, (s_trng s,a_char t));
_PROTOTYPE( a_bool s_scne, (s_trng s,a_char t));
_PROTOTYPE( a_bool s_cseq, (a_char s,s_trng t));
_PROTOTYPE( a_bool s_csge, (a_char s,s_trng t));
_PROTOTYPE( a_bool s_csgt, (a_char s,s_trng t));
_PROTOTYPE( a_bool s_csle, (a_char s,s_trng t));
_PROTOTYPE( a_bool s_cslt, (a_char s,s_trng t));
_PROTOTYPE( a_bool s_csne, (a_char s,s_trng t));





