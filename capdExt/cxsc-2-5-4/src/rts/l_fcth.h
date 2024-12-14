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

/* CVS $Id: l_fcth.h,v 1.21 2014/01/30 17:24:10 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : l_fcth.h                              */
/*                                                              */
/*      Description     : function prototypes                   */
/*                                                              */
/*                 : intro of _PROTO..=b=                       */
/*                       l_rprc(),l_whex()                      */
/****************************************************************/

_PROTOTYPE(multiprecision l_abs,(multiprecision a));
_PROTOTYPE(multiprecision l_acos,(multiprecision a));
_PROTOTYPE(multiprecision l_acot,(multiprecision a));
_PROTOTYPE(multiprecision l_acsh,(multiprecision a));
_PROTOTYPE(multiprecision l_acth,(multiprecision a));
_PROTOTYPE(multiprecision l_addc,(multiprecision i1,multiprecision i2));
_PROTOTYPE(multiprecision l_addd,(multiprecision i1,multiprecision i2));
_PROTOTYPE(multiprecision l_addu,(multiprecision i1,multiprecision i2));
_PROTOTYPE(multiprecision l_asin,(multiprecision a));
_PROTOTYPE(multiprecision l_asnh,(multiprecision a));
_PROTOTYPE(void           l_ass ,(multiprecision *a,multiprecision b));
_PROTOTYPE(multiprecision l_atan,(multiprecision a));
_PROTOTYPE(multiprecision l_atn2,(multiprecision a,multiprecision b));
_PROTOTYPE(multiprecision l_atnh,(multiprecision a));
_PROTOTYPE(multiprecision l_comp,(multiprecision a,a_intg n));
_PROTOTYPE(void           l_conv,(a_char *strng,multiprecision *s,a_intg rnd));
_PROTOTYPE(multiprecision l_cos ,(multiprecision a));
_PROTOTYPE(multiprecision l_cosh,(multiprecision a));
_PROTOTYPE(multiprecision l_cot ,(multiprecision a));
_PROTOTYPE(multiprecision l_coth,(multiprecision a));
_PROTOTYPE(multiprecision l_divc,(multiprecision i1,multiprecision i2));
_PROTOTYPE(multiprecision l_divd,(multiprecision i1,multiprecision i2));
_PROTOTYPE(multiprecision l_divu,(multiprecision i1,multiprecision i2));
_PROTOTYPE(a_btyp         b_dtol,(d_otpr d,multiprecision *m,a_intg rnd));
_PROTOTYPE(a_bool         l_eq,(multiprecision a,multiprecision b));
_PROTOTYPE(void           l_exct,(multiprecision *a,multiprecision b,a_intg *r,a_intg *l));
_PROTOTYPE(a_intg         l_expo,(multiprecision a));
_PROTOTYPE(multiprecision l_exp ,(multiprecision a));
_PROTOTYPE(multiprecision l_flot,(a_intg n));
_PROTOTYPE(void           l_free,(multiprecision *a));
_PROTOTYPE(a_bool         l_ge,(multiprecision a,multiprecision b));
_PROTOTYPE(a_bool         l_gt,(multiprecision a,multiprecision b));
_PROTOTYPE(void           l_init,(multiprecision *a));
_PROTOTYPE(a_bool         l_le,(multiprecision a,multiprecision b));
_PROTOTYPE(multiprecision l_log ,(multiprecision a));
_PROTOTYPE(multiprecision l_loga,(multiprecision a,multiprecision b));
_PROTOTYPE(a_bool         l_lt,(multiprecision a,multiprecision b));
_PROTOTYPE(a_btyp         b_ltod,(multiprecision m,d_otpr *d,a_intg rnd));
_PROTOTYPE(multiprecision l_mant,(multiprecision a));
_PROTOTYPE(a_intg         l_mlen,(multiprecision a));
_PROTOTYPE(multiprecision l_mulc,(multiprecision i1,multiprecision i2));
_PROTOTYPE(multiprecision l_muld,(multiprecision i1,multiprecision i2));
_PROTOTYPE(multiprecision l_mulu,(multiprecision i1,multiprecision i2));
_PROTOTYPE(a_bool         l_ne,(multiprecision a,multiprecision b));
_PROTOTYPE(multiprecision l_pow ,(multiprecision a,multiprecision b));
_PROTOTYPE(void           l_prec,(a_intg n));
_PROTOTYPE(multiprecision l_pred,(multiprecision a));
_PROTOTYPE(a_intg         l_rprc,(void));
_PROTOTYPE(void           l_read,(FILE *device,multiprecision *s,a_intg rnd,a_intg ii));
_PROTOTYPE(a_intg         l_rond,(multiprecision a));
_PROTOTYPE(void           l_rval,(s_trng strng,multiprecision *s,a_intg rnd));
_PROTOTYPE(a_intg         l_sign,(multiprecision a));
_PROTOTYPE(multiprecision l_sin ,(multiprecision a));
_PROTOTYPE(multiprecision l_sinh,(multiprecision a));
_PROTOTYPE(multiprecision l_sqrt,(multiprecision a));
_PROTOTYPE(multiprecision l_subc,(multiprecision i1,multiprecision i2));
_PROTOTYPE(multiprecision l_subd,(multiprecision i1,multiprecision i2));
_PROTOTYPE(multiprecision l_subu,(multiprecision i1,multiprecision i2));
_PROTOTYPE(multiprecision l_succ,(multiprecision a));
_PROTOTYPE(multiprecision l_tan ,(multiprecision a));
_PROTOTYPE(multiprecision l_tanh,(multiprecision a));
_PROTOTYPE(void           l_temp,(multiprecision *i1));
_PROTOTYPE(a_intg         l_trun,(multiprecision a));
_PROTOTYPE(a_bool         l_ttmp,(multiprecision *a));
_PROTOTYPE(multiprecision l_umin,(multiprecision i1));
_PROTOTYPE(void           l_utmp,(multiprecision *i1));
_PROTOTYPE(void           l_vlcp,(multiprecision *i1));
_PROTOTYPE(void           l_whex,(f_text *desc,multiprecision r,a_char mode));
_PROTOTYPE(void           l_writ,(f_text *desc,multiprecision s, a_intg 
						    TotalWidth, a_intg FracDigits,a_intg rnd));





