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

/* CVS $Id: d_fcth.h,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : d_fcth.h                              */
/*                                                              */
/*      Description     : function prototypes                   */
/*                                                              */
/*                       i_dsta, z_dadd,z_dsta,z_dsub    =b=    */
/****************************************************************/
_PROTOTYPE(void         d_assi,(d_otpi *a,d_otpi b));
_PROTOTYPE(void         d_assc,(d_otpc *a,d_otpc b));
_PROTOTYPE(void         d_assz,(d_otpz *a,d_otpz b));

#ifdef LINT_ARGS
void         c_cadd(dotprecision *cr,dotprecision *ci,a_cmpx a);
void         c_csub(dotprecision *cr,dotprecision *ci,a_cmpx a);
void         c_dadd(dotprecision *cr,dotprecision *ci,d_otpc a);
d_otpc       c_dsta(dotprecision cr,dotprecision ci);
void         c_dsub(dotprecision *cr,dotprecision *ci,d_otpc a);
void         c_padd(dotprecision *cr,dotprecision *ci,a_cmpx a,a_cmpx b);
void         c_psub(dotprecision *cr,dotprecision *ci,a_cmpx a,a_cmpx b);
void         c_rcad(dotprecision *cr,dotprecision *ci,a_real a,a_cmpx b);
void         c_rcsb(dotprecision *cr,dotprecision *ci,a_real a,a_cmpx b);
a_cmpx       c_scps(a_cmpx r[],a_cmpx s[],a_intg n,a_intg rnd);
a_cmpx       c_scpy(y_desc *r,y_desc *s,a_intg rnd);
a_cmpx       c_stad(dotprecision cr,dotprecision ci);
a_cmpx       c_stan(dotprecision cr,dotprecision ci);
a_cmpx       c_stau(dotprecision cr,dotprecision ci);
dotprecision d_add(dotprecision a,dotprecision b);
void         d_ass(dotprecision *a,dotprecision b);
dotprecision d_chs(dotprecision a);
void         d_clr(dotprecision *a);
void         d_dadd(dotprecision *a,dotprecision b);
dotprecision d_dot(a_real *a,a_real *b,a_intg n);
a_intg       d_sign(dotprecision a);
void         d_dsub(dotprecision *a,dotprecision b);
a_bool       d_eq(dotprecision a,dotprecision b);
void         d_free(dotprecision *a);
a_bool       d_ge(dotprecision a,dotprecision b);
a_bool       d_gt(dotprecision a,dotprecision b);
void         d_init(dotprecision *a);
a_bool       d_le(dotprecision a,dotprecision b);
a_bool       d_lt(dotprecision a,dotprecision b);
a_bool       d_ne(dotprecision a,dotprecision b);
void         d_padd(dotprecision *c,a_real a,a_real b);
void         d_psub(dotprecision *c,a_real a,a_real b);
void         d_radd(dotprecision *c,a_real a);
void         d_rsub(dotprecision *c,a_real a);
a_real       d_stad(dotprecision a);
a_real       d_stan(dotprecision a);
a_real       d_stau(dotprecision a);
dotprecision d_sub(dotprecision a,dotprecision b);
dotprecision d_sum(a_real *a,a_intg n);
void         d_temp(dotprecision *d);
void         d_utmp(dotprecision *d);
void         d_vlcp(dotprecision *d);
void         i_dadd(dotprecision *cl,dotprecision *cu ,d_otpi a);
d_otpi       i_dsta(dotprecision cl,dotprecision cu);
void         i_dsub(dotprecision *cl,dotprecision *cu ,d_otpi a);
void         i_iadd(dotprecision *cl,dotprecision *cu,a_intv a);
a_intv       i_ista(dotprecision cl,dotprecision cu);
void         i_isub(dotprecision *cl,dotprecision *cu,a_intv a);
void         i_padd(dotprecision *cl,dotprecision *cu,a_intv a,a_intv b);
void         i_psub(dotprecision *cl,dotprecision *cu,a_intv a,a_intv b);
void         i_riad(dotprecision *cl,dotprecision *cu,a_real a,a_intv b);
void         i_risb(dotprecision *cl,dotprecision *cu,a_real a,a_intv b);
a_intv       i_scps(a_intv r[],a_intv s[],a_intg n,a_intg rnd);
a_intv       i_scpy(y_desc *r,y_desc *s,a_intg rnd);
a_real       r_scps(a_real r[],a_real s[],a_intg n,a_intg rnd);
a_real       r_scpy(y_desc *r,y_desc *s,a_intg rnd);
void         z_ciad(dotprecision *crl,dotprecision *cil,dotprecision *cru,
                    dotprecision *ciu,a_cmpx a,a_intv b);
void         z_cisb(dotprecision *crl,dotprecision *cil,dotprecision *cru,
                    dotprecision *ciu,a_cmpx a,a_intv b);
void         z_czad(dotprecision *crl,dotprecision *cil,dotprecision *cru,
                    dotprecision *ciu,a_cmpx a,a_cinv b);
void         z_czsb(dotprecision *crl,dotprecision *cil,dotprecision *cru,
                    dotprecision *ciu,a_cmpx a,a_cinv b);
void         z_dadd(dotprecision *ZRI,dotprecision* ZII,dotprecision* ZRS,
                    dotprecision *ZIS,d_otpz A);
d_otpz       z_dsta(dotprecision ZRI,dotprecision ZII ,dotprecision ZRS,
                    dotprecision ZIS);
void         z_dsub(dotprecision *ZRI,dotprecision* ZII,dotprecision* ZRS,
                    dotprecision *ZIS,d_otpz A);
void         z_izad(dotprecision *crl,dotprecision *cil,dotprecision *cru,
                    dotprecision *ciu,a_intv a,a_cinv b);
void         z_izsb(dotprecision *crl,dotprecision *cil,dotprecision *cru,
                    dotprecision *ciu,a_intv a,a_cinv b);
void         z_padd(dotprecision *crl,dotprecision *cil,dotprecision *cru,
                    dotprecision *ciu,a_cinv a,a_cinv b);
void         z_psub(dotprecision *crl,dotprecision *cil,dotprecision *cru,
                    dotprecision *ciu,a_cinv a,a_cinv b);
void         z_rzad(dotprecision *crl,dotprecision *cil,dotprecision *cru,
                    dotprecision *ciu,a_real a,a_cinv b);
void         z_rzsb(dotprecision *crl,dotprecision *cil,dotprecision *cru,
                    dotprecision *ciu,a_real a,a_cinv b);
a_cinv       z_scps(a_cinv r[],a_cinv s[],a_intg n,a_intg rnd);
a_cinv       z_scpy(y_desc *r,y_desc *s,a_intg rnd);
void         z_zadd(dotprecision *crl,dotprecision *cil,dotprecision *cru,
                    dotprecision *ciu,a_cinv a);
a_cinv       z_zsta(dotprecision crl,dotprecision cil,dotprecision cru,
                    dotprecision ciu);
void         z_zsub(dotprecision *crl,dotprecision *cil,dotprecision *cru,
                    dotprecision *ciu,a_cinv a);
#else
a_cinv       z_scps(), z_scpy(), z_zsta();
a_cmpx       c_scps(), c_scpy(), c_stad(), c_stan(), c_stau();
a_intv       i_ista(), i_scps(), i_scpy();
a_real       d_stad(), d_stan(), d_stau(), r_scps(), r_scpy();
a_bool       d_eq(),   d_ge(),   d_gt(),   d_le(),   d_lt(),   d_ne();
void         c_cadd(), c_csub(), c_dadd(), c_dsub(), c_padd(), c_psub(), 
             c_rcad(), c_rcsb(),
             d_ass(),  d_clr(),  d_dadd(), d_dsub(), d_free(), d_init(),
             d_padd(), d_psub(), d_radd(), d_rsub(), d_temp(),
             d_utmp(), d_vlcp(), 
             i_dadd(), i_dsub(), i_iadd(), i_isub(), i_padd(), i_psub(), 
             i_riad(), i_risb(), 
             z_ciad(), z_cisb(), z_czad(), z_czsb(), z_dadd(), z_dsub(),
             z_izad(), z_izsb(), z_padd(), z_psub(), z_rzad(), z_rzsb(), 
             z_zadd(), z_zsub();
dotprecision d_add(),  d_chs(),  d_dot(),  d_sub(),  d_sum();
a_intg       d_sign();
d_otpc       c_dsta();
d_otpi       i_dsta();
d_otpz       z_dsta();
#endif





