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

/* CVS $Id: f_fcth.h,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_fcth.h                              */
/*                                                              */
/*      Description     : function prototypes                   */
/*                                                              */
/****************************************************************/
_PROTOTYPE( void   p_init,(int argc,char **argv));

_PROTOTYPE(void   f_args,(s_trng *name));
_PROTOTYPE(void   f_assg,(f_text *desc,char *name,size_t len));
_PROTOTYPE(void   f_back,(f_text *desc));
_PROTOTYPE(void   f_bhex,(f_text *desc,a_char r,a_char mode));
_PROTOTYPE(void   f_eofp,(void));
_PROTOTYPE(a_bool f_exst,(a_char name[]));
_PROTOTYPE(void   f_free,(f_text *desc));
_PROTOTYPE(void   f_get_,(f_text *desc));
_PROTOTYPE(void   f_getc,(f_text *desc));
_PROTOTYPE(a_intg f_op88,(f_text *desc,s_trng name,a_intg level));
_PROTOTYPE(void   f_popt,(s_trng *name));
_PROTOTYPE(void   f_put_,(f_text *desc));
_PROTOTYPE(void   f_putc,(a_char ch,f_text *desc));
_PROTOTYPE(void   f_quer,(f_text *desc,a_intg *status,s_trng *name));
_PROTOTYPE(void   f_rdf1,(f_text *desc,a_real *r));
_PROTOTYPE(void   f_rdh1,(f_text *desc,a_char *h));
_PROTOTYPE(void   f_rdi1,(f_text *desc,a_intg *i));
_PROTOTYPE(void   f_rdl1,(f_text *desc,multiprecision *r));
_PROTOTYPE(void   f_rdl2,(f_text *desc,multiprecision *r,a_intg rnd));
_PROTOTYPE(void   f_rdln,(f_text *desc));
_PROTOTYPE(void   f_rdr1,(f_text *desc,a_real *r));
_PROTOTYPE(void   f_rdr2,(f_text *desc,a_real *r,a_intg rnd));
_PROTOTYPE(void   f_rds1,(f_text *desc,s_trng *s));
_PROTOTYPE(void   f_rdv1,(f_text *desc,a_intv *r));
_PROTOTYPE(void   f_read,(f_text *desc,a_VOID value));
_PROTOTYPE(void   f_rhex,(f_text *desc,a_real *r,a_char mode));
_PROTOTYPE(a_intg f_rint,(FILE *device,a_intg *i));
_PROTOTYPE(void   f_rset,(f_text *desc,a_char *name,a_char *device));
_PROTOTYPE(void   f_rstn,(f_text *desc,a_intg spec));
_PROTOTYPE(void   f_rwri,(f_text *desc,a_char *name,a_char *device));
_PROTOTYPE(void   f_rwrn,(f_text *desc,a_intg spec));
_PROTOTYPE(a_bool f_sexs,(s_trng name));
_PROTOTYPE(void   f_srse,(f_text *desc,s_trng device));
_PROTOTYPE(void   f_srwi,(f_text *desc,s_trng device));
_PROTOTYPE(void   f_whex,(f_text *desc,a_real r,a_char mode));
_PROTOTYPE(void   f_wint,(f_text *device,a_intg i,a_intg TotalWidth));
_PROTOTYPE(void   f_wrb1,(f_text *desc,a_bool b));
_PROTOTYPE(void   f_wrb2,(f_text *desc,a_bool b,a_intg w));
_PROTOTYPE(void   f_wrc1,(f_text *desc,a_char s[],a_intg n));
_PROTOTYPE(void   f_wrc2,(f_text *desc,a_char s[],a_intg n,a_intg w));
_PROTOTYPE(void   f_wrd1,(f_text *desc,d_otpr r,a_char mode));
_PROTOTYPE(void   f_wrf1,(f_text *desc,a_real r));
_PROTOTYPE(void   f_wrf2,(f_text *desc,a_real r,a_intg w));
_PROTOTYPE(void   f_wrf3,(f_text *desc,a_real r,a_intg w,a_intg f));
_PROTOTYPE(void   f_wrf4,(f_text *desc,a_real r,a_intg w,a_intg f,a_intg rnd));
_PROTOTYPE(void   f_wrh1,(f_text *desc,a_char h));
_PROTOTYPE(void   f_wrh2,(f_text *desc,a_char h,a_intg w));
_PROTOTYPE(void   f_wri1,(f_text *desc,a_intg i));
_PROTOTYPE(void   f_wri2,(f_text *desc,a_intg i,a_intg w));
_PROTOTYPE(void   f_wrid,(FILE   *desc,d_otpr r,a_char mode));
_PROTOTYPE(void   f_writ,(f_text *desc,a_VOID value));
_PROTOTYPE(void   f_wrl1,(f_text *desc,multiprecision r));
_PROTOTYPE(void   f_wrl2,(f_text *desc,multiprecision r,a_intg w));
_PROTOTYPE(void   f_wrl3,(f_text *desc,multiprecision r,a_intg w,a_intg f));
_PROTOTYPE(void   f_wrl4,(f_text *desc,multiprecision r,a_intg w,a_intg f,a_intg rnd));
_PROTOTYPE(void   f_wrln,(f_text *desc));
_PROTOTYPE(void   f_wrr1,(f_text *desc,a_real r));
_PROTOTYPE(void   f_wrr2,(f_text *desc,a_real r,a_intg w));
_PROTOTYPE(void   f_wrr3,(f_text *desc,a_real r,a_intg w,a_intg f));
_PROTOTYPE(void   f_wrr4,(f_text *desc,a_real r,a_intg w,a_intg f,a_intg rnd));
_PROTOTYPE(void   f_wrs1,(f_text *desc,s_trng s));
_PROTOTYPE(void   f_wrs2,(f_text *desc,s_trng s,a_intg w));





