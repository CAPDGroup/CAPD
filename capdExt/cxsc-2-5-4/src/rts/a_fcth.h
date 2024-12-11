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

/* CVS $Id: a_fcth.h,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_fcth.h                              */
/*                                                              */
/*      Description     : function prototypes                   */
/*                                                              */
/****************************************************************/

#ifdef LINT_ARGS
void   a_2pop(void);
void   a_2psh(a_VOID new_d,a_VOID old);
a_intg a_abs_(a_intg a);
a_intg a_add_(a_intg i,a_intg j);
a_intg a_band(a_intg a,a_intg b);
a_intg a_bclr(a_intg a,a_intg m);
a_intg a_beqv(a_intg a,a_intg b);
a_intg a_bmsb(a_intg a);
a_intg a_bnot(a_intg a);
a_intg a_bons(a_intg a);
a_intg a_bor_(a_intg a,a_intg b);
a_intg a_brtt(a_intg a,a_intg m);
a_intg a_bset(a_intg a,a_intg m);
a_intg a_bshf(a_intg a,a_intg m);
a_bool a_btst(a_intg a,a_intg m);
a_intg a_bxor(a_intg a,a_intg b);
a_intg a_clck(void);
void   a_cpsh(a_VOID new_d,a_VOID old);
a_intg a_div_(a_intg i,a_intg j);
void   a_exit(a_intg p);
void   a_free(char **);
a_btyp a_gets(a_btyp *rc);
a_intg a_gtim(void);
void   a_itim(void);
a_intg a_ival(s_trng s);
a_btyp a_ixch(a_intg i,a_intg l,a_intg u);
a_VOID a_lloc(size_t i);
a_intg a_mod_(a_intg i,a_intg j);
a_intg a_mul_(a_intg i,a_intg j);
a_VOID a_nilc(a_VOID a);
void   a_popt(f_text *desc,s_trng option);
void   a_sets(int rounding);
a_intg a_sqr_(a_intg a);
a_intg a_sub_(a_intg i,a_intg j);
a_intg a_sval(s_trng s,s_trng *t);
a_intg a_syst(s_trng s);
a_intg a_tick(void);
a_intg a_umin(a_intg a);
#else
a_bool a_btst();
a_btyp a_gets(), a_ixch();
a_intg a_abs_(), a_add_(), a_band(), a_bclr(), a_beqv(), a_bmsb(), a_bnot(),
       a_bons(), a_bor_(), a_brtt(), a_bset(), a_bshf(), a_bxor(), a_clck(),
       a_div_(), a_gtim(), a_ival(), a_mod_(), a_mul_(), a_sqr_(), a_sub_(), 
	  a_sval(), a_syst(), a_tick(), a_umin();
a_VOID a_lloc(), a_nilc(), a_ppop();
void   a_exit(), a_free(), a_itim(), a_popt(), a_ppsh(), a_sets();
#endif





