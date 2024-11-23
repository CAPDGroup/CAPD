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

/* CVS $Id: b_fcth.h,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_fcth.h                              */
/*                                                              */
/*      Description     : function prototypes                   */
/*                                                              */
/*                   b_conf(),b_coni(),b_form() new arguments   */
/*                   b_adpp() new                               */
/*                   VAX_VMS_C considered / b_imag() removed    */
/*                   b_bcdi(),b_bcid() arguments changed        */
/*                   standard functions prototypes              */
/****************************************************************/

#ifdef LINT_ARGS
#if VAX_VMS_C
void     b_accu(void);
#endif
int      b_acos(multiprecision xi,multiprecision ri);
int      b_acot(multiprecision xi,multiprecision ri);
int      b_acsh(multiprecision xi,multiprecision ri);
int      b_acth(multiprecision xi,multiprecision ri);
void     b_addc(a_btyp *a);
a_bool   b_addm(a_intg k,a_btyp *a,a_btyp *b);
void     b_addu(a_btyp a,a_btyp b,a_intg cin,a_btyp *c,a_intg *cout);
a_btyp   b_adj (a_btyp *lang,a_intg *expo);
a_btyp   b_adpp(a_btyp **dc,a_intg *size,a_intg expo,a_intg dp,
                a_intg length,a_intg *p_start,a_intg *p_dp,
                a_intg *p_length);
int      b_asgn(multiprecision res);
int      b_asin(multiprecision xi,multiprecision ri);
int      b_asiv(multiprecision xi,multiprecision LPiOv2);
int      b_asnh(multiprecision xi,multiprecision ri);
int      b_atan(multiprecision xi,multiprecision ri);
int      b_atav(multiprecision xi,multiprecision LPiOv2);
int      b_atn2(multiprecision xi,multiprecision yi,multiprecision ri);
int      b_atnh(multiprecision xi,multiprecision ri);
int      b_baad(multiprecision i1,multiprecision i2,multiprecision r);
int      b_bacm(multiprecision i1,multiprecision i2);
int      b_badd(multiprecision i1,multiprecision i2,multiprecision r);
int      b_badj(a_btyp n,multiprecision i);
int      b_ball(a_btyp n,a_btyp **i);
int      b_banx(multiprecision i,multiprecision r);
int      b_basu(multiprecision i1,multiprecision i2,multiprecision r);
int      b_bcad(a_intg n,a_btyp *r);
int      b_bcat(a_intg n,a_btyp *r);
int      b_bcdi(a_real d,multiprecision *i,a_intg rnd);
int      b_bcid(multiprecision i,a_real *d,a_intg rnd);
int      b_bclr(multiprecision i);
int      b_bcmp(multiprecision i1,multiprecision i2);
int      b_bcpy(multiprecision i,multiprecision r);
int      b_bcsu(a_intg n,a_btyp *r);
int      b_bdiv(multiprecision i1,multiprecision i2,multiprecision r);
int      b_bdvn(multiprecision i1,a_btyp n,multiprecision r);
int      b_bini(multiprecision i);
int      b_biv_(a_intv a);
a_btyp   b_biv2(int (*func)(multiprecision,multiprecision,multiprecision),
                       a_real a1, a_real a2, a_real *rlb, a_real *rub);
a_btyp   b_bivp(int (*func)(multiprecision,multiprecision),
                       a_real arg, a_real *rlb, a_real *rub);
a_btyp   b_bldx(a_real a,a_real *res);
a_btyp   b_blgx(a_real a,a_real *res);
int      b_bmat(a_intg n,a_btyp *i1,a_btyp *i2,a_intg carry,a_intg *r);
int      b_bmcm(a_btyp n,a_btyp *i1,a_btyp *i2);
int      b_bmdv(a_btyp *i1,a_btyp *i2,a_btyp *r);
int      b_bms1(a_btyp i,a_btyp *r);
int      b_bms2(a_btyp i,a_btyp *r);
int      b_bmsh(a_intg l,a_btyp *s,a_intg n);
int      b_bmsp(a_intg n,a_btyp *i,a_btyp u,a_btyp *r);
int      b_bmts(a_intg n,a_btyp *i);
int      b_bmul(multiprecision i1,multiprecision i2,multiprecision r1,
                multiprecision r2);
int      b_bmun(multiprecision i,a_btyp n,multiprecision r);
int      b_bnxt(multiprecision i,multiprecision r);
int      b_bpnt( a_intv a );
int      b_brnd(multiprecision i,multiprecision r);
int      b_bshf(a_intg n,multiprecision i,multiprecision r);
int      b_bsub(multiprecision i1,multiprecision i2,multiprecision r);
int      b_busp(a_btyp i,a_btyp u,a_btyp *r);
int      b_chck(a_char *strng,char **buffer,a_intg *size,
                       a_intg *expo,a_intg *dp,a_intg *length,
                       a_bool *sign, a_char **next);
#if VAX_VMS_C
void     b_cmcp(void);
#endif
void     b_comp(a_real *x,a_intg expo,a_btyp *mant,a_bool vz);
void     b_conf(a_intg dp,a_btyp *pd,a_intg *a_b,a_intg *a_e,d_otpr d,
                a_intg *bits);
void     b_coni(a_intg dp,a_btyp *pd,a_intg *a_b,a_intg *a_e,d_otpr d,
                a_intg *bits);
int      b_cos_(multiprecision xi,multiprecision ri);
int      b_cosh(multiprecision xi,multiprecision ri);
int      b_cot_(multiprecision xi,multiprecision ri);
int      b_coth(multiprecision xi,multiprecision ri);
int      b_cscn(FILE *device,char *buffer,a_intg *expo,
                       a_intg *dp,a_intg *length,a_bool *sign,int i);
a_bool   b_deko(a_real x,a_intg *expo,a_btyp *mant,a_bool *vz);
int      b_drop(int n);
int      b_errr(a_btyp err);
int      b_exp_(multiprecision xi,multiprecision ri);
int      b_expe(multiprecision Larg);
a_btyp   b_form(a_btyp *dc,a_intg *size,a_intg expo,a_intg dp,
                a_intg length,a_bool sign,a_intg rnd,a_real *s);
void     b_freh(a_char *ptr,a_char *addr,a_char *where);
multiprecision b_get_(void);
a_bool   b_geta(dotprecision a,a_btyp *result,a_intg *expo,a_bool *vz);
void     b_geth(a_char *ptr,a_char *addr,a_char *where);
int      b_gini(void);
#if VAX_VMS_C
void     b_glbl(void);
#endif
a_btyp   b_ifrm(a_btyp *dc,a_intg expo,a_intg dp,a_intg length,
                a_bool sign,a_intv *s);
a_btyp   b_inv1(int (*func)(multiprecision,multiprecision),a_real arg,
                a_real *res,a_intg rnd);
a_btyp   b_inv2(int (*func)(multiprecision,multiprecision,multiprecision),
                a_real a1,a_real a2,a_real *res,a_intg rnd);
void     b_irnd(multiprecision value,multiprecision *lb,multiprecision *ub);
a_real   b_klog(a_real x);
a_real   b_ksqt(a_real x);
a_real   b_katn(a_real x);
a_real   b_kasn(a_real x);
int      b_lnva(multiprecision t);
int      b_lnve(multiprecision Larg);
int      b_log_(multiprecision xi,multiprecision ri);
int      b_loga(multiprecision xi,multiprecision bi,multiprecision ri);
a_btyp   b_ltor(multiprecision i,a_real *r,a_intg rnd);
void     b_mdiv(a_btyp *ma,a_btyp *mb,a_btyp *mc,a_intg *expo);
void     b_muad(a_btyp a,a_btyp b,a_btyp *lang);
a_intg   b_op88(f_text *desc,s_trng name,a_intg level);
void     b_out (a_btyp *mant,a_intg expo,a_intg digits,char *buffer,
                a_intg *bdp,a_intg *dexpo);
void     b_outf(a_intg *digits,char *buffer,a_intg *bdp,a_intg *dexpo,
                dotprecision c);
void     b_outi(a_intg *digits,char *buffer,a_intg *bdp,a_intg *dexpo,
                dotprecision c);
void     b_outm(a_btyp *mant,a_intg len,a_intg expo,a_intg digits,
                char *buffer,a_intg *bdp,a_intg *dexpo);
int      b_pi__(multiprecision pi);
int      b_pign(void);
int      b_popt(FILE *device,char *option);
int      b_pow_(multiprecision xi,multiprecision yi,multiprecision ri);
void     b_prod(a_btyp *a,a_btyp *b,a_btyp *lang);
void     b_rnd (a_intg rnd,char *buffer,a_intg digits,a_intg pos,
                a_intg *bdp,a_intg *dexpo);
a_btyp   b_rndd(a_btyp *lang,a_intg *expo,a_bool vz);
a_btyp   b_rndn(a_btyp *lang,a_intg *expo);
a_btyp   b_rndu(a_btyp *lang,a_intg *expo,a_bool vz);
a_btyp   b_rtol(a_real r,multiprecision *i,a_intg rnd);
int      b_scan(FILE *device,char **buffer,a_intg *size,a_intg *expo,
                a_intg *dp,a_intg *length,a_bool *sign,int i);
int      b_sscn(char *strng,char *buffer,a_intg *expo,a_intg *dp,
                a_intg *length,a_bool *sign,a_intg *proc);
void     b_shl1(a_btyp *lang,a_intg laenge);
void     b_shlu(a_btyp *lang,a_intg laenge,a_intg dist);
void     b_shr1(a_btyp *lang,a_intg laenge);
void     b_shru(a_btyp *lang,a_intg laenge,a_intg dist);
int      b_sico(multiprecision Larg);
int      b_sin_(multiprecision xi,multiprecision ri);
int      b_sinh(multiprecision xi,multiprecision ri);
int      b_snhv(multiprecision xi);
int      b_sqrt(multiprecision xi,multiprecision ri);
int      b_sqrv(multiprecision Larg);
void     b_subc(a_btyp *a);
a_bool   b_subm(a_intg k,a_btyp *a,a_btyp *b);
void     b_subu(a_btyp a,a_btyp b,a_intg bin,a_btyp *c,a_intg *bout);
int      b_tadd(tenbyte *a,tenbyte *b,tenbyte *res);
void     b_tadj(a_btyp *mant,a_intg *expo);
int      b_tan_(multiprecision xi,multiprecision ri);
int      b_tanh(multiprecision xi,multiprecision ri);
void     b_tcom(tenbyte *a,a_intg expo,a_btyp *mant,a_bool vz);
a_bool   b_tdek(tenbyte *a,a_intg *expo,a_btyp *mant,a_bool *vz);
int      b_tdiv(tenbyte *a,tenbyte *b,tenbyte *res);
a_bool   b_test(a_intg n,a_btyp *a);
a_bool   b_text(f_text *desc,a_bool in);
void     b_tmpf(char *name,a_intg length);
void     b_tmdv(a_btyp *ma,a_btyp *mb,a_btyp *mc,a_intg *expo);
void     b_tmph(a_char *ptr);
int      b_tmul(tenbyte *a,tenbyte *b,tenbyte *res);
void     b_trnd(a_btyp *mant,a_intg *expo,a_bool vz);
int      b_tsub(tenbyte *a,tenbyte *b,tenbyte *res);
void     b_varh(a_char *ptr,a_char *addr);
void     o_user(void);
#else
#if VAX_VMS_C
void     b_accu(), b_cmcp(), b_glbl(), o_user();
#endif
a_btyp   b_adj(),  b_form(), b_ifrm(), b_inv1(), b_inv2(), b_ltor(),
         b_rndd(), b_rndn(), b_rndu(), b_rtol();
a_bool   b_addm(), b_deko(), b_geta(), b_subm(),
         b_tdek(), b_test(), b_text();
int      b_acos(), b_acot(), b_acsh(), b_acth(), b_asgn(), b_asin(),
         b_asiv(), b_asnh(), b_atan(), b_atav(), b_atn2(), b_atnh(),
         b_baad(), b_badd(), b_bacm(), b_badj(), b_ball(), b_banx(),
         b_basu(), b_bcad(), b_bcat(), b_bcdi(), b_bcid(), b_bclr(),
         b_bcmp(), b_bcpy(), b_bcsu(), b_bdiv(), b_bdvn(), b_bini(),
         b_bmad(), b_bmap(), b_bmat(), b_brnd(), b_bmcm(), b_bmcp(),
         b_bmdv(), b_bmsh(), b_bmsp(), b_bms1(), b_bms2(), b_bmts(),
         b_bmul(), b_bmun(), b_bnxt(), b_bshf(), b_bsub(), b_busp(),
         b_chck(), b_cos_(), b_cosh(), b_cot_(), b_coth(), b_cscn(),
         b_drop(), b_errr(), b_exp_(), b_expe(), b_gini(), b_lnva(),
         b_lnve(), b_log_(), b_loga(), b_pi__(), b_pign(),
         b_popt(), b_pow_(),
         b_scan(), b_sscn(), b_sico(), b_sin_(), b_sinh(), b_snhv(),
         b_sqrt(), b_sqrv(),
         b_tadd(), b_tan_(), b_tanh(),
         b_tdiv(), b_tmul(), b_tsub();
multiprecision b_get_();
void     b_addc(), b_addu(), b_cnfu(), b_cniu(), b_comp(),
         b_conf(), b_coni(), b_freh(), b_geth(), b_irnd(),
         b_mdiv(), b_muad(),
         b_out(),  b_outf(), b_outi(), b_outm(), b_prod(),
         b_rnd(),  b_shlu(), b_shl1(),
         b_shru(), b_shr1(), b_subc(), b_subu(),
         b_tadj(), b_tcom(), b_tmdv(),
         b_tmph(), b_tmpf(), b_trnd(), b_varh();
a_btyp   b_adpp(), b_biv2(), b_bivp(), b_bldx(), b_blgx();
int      b_biv_(), b_bpnt();
a_intg   b_op88();
a_real   b_klog(), b_ksqt(), b_katn(), b_kasn();
#endif





