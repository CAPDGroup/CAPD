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

/* CVS $Id: e_fcth.h,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_fcth.h                              */
/*                                                              */
/*      Description     : function prototypes                   */
/*                                                              */
/****************************************************************/

#ifdef LINT_ARGS
void    e_actn(a_btyp action,a_btyp code,
               void (*function)(a_btyp,int,va_list));
void    e_back(FILE *device);
void    e_bmsg(FILE *device);
#if VAX_VMS_C
void    e_data(void);
#endif
a_bool  e_dz_e(void);
a_bool  e_dz_o(void);
a_bool  e_ie_e(void);
a_bool  e_ie_o(void);
void    e_ienv(a_intg action,a_intg code,a_bool mode);
a_bool  e_io_e(void);
a_bool  e_io_o(void);
a_bool  e_of_e(void);
a_bool  e_of_o(void);
void    e_popp(void);
void    e_push(char *function,char *filename);
void    e_rall(void);
void    e_rdze(void);
void    e_rdzo(void);
void    e_rest(a_intg d);
void    e_riee(void);
void    e_rieo(void);
void    e_rioe(void);
void    e_rioo(void);
void    e_rofe(void);
void    e_rofo(void);
void    e_rufe(void);
void    e_rufo(void);
void    e_sdze(void);
void    e_sall(void);
void    e_save(a_intg *p);
void    e_sdzo(void);
void    e_siee(void);
void    e_sieo(void);
void    e_sioe(void);
void    e_sioo(void);
void    e_sofe(void);
void    e_sofo(void);
void    e_sufe(void);
void    e_sufo(void);
void    e_tmsg(int msgid);
void    e_tmrt(int e_argc,va_list e_argv,a_bool print);
void    e_tprt(int e_argc,va_list e_argv);
#if C_P_7
void    e_trap(a_btyp code,int e_argc,...);
#else
void    e_trap();
#endif
a_bool  e_uf_e(void);
a_bool  e_uf_o(void);
void    e_xall(a_btyp,int,va_list);
void    e_xarg(a_btyp,int,va_list);
void    e_xdbz(a_btyp,int,va_list);
void    e_xine(a_btyp,int,va_list);
void    e_xiob(a_btyp,int,va_list);
void    e_xios(a_btyp,int,va_list);
void    e_xiop(a_btyp,int,va_list);
void    e_xnor(a_btyp,int,va_list);
void    e_xofl(a_btyp,int,va_list);
void    e_xset(a_btyp,int,va_list);
void    e_xufl(a_btyp,int,va_list);
#else
#if VAX_VMS_C
void    e_data();
#endif
void    e_actn(), e_back(), e_ienv(), e_popp(), e_push(), e_trap(),
        e_trpt(), e_sioo(), e_sdzo(), e_sofo(), e_sufo(), e_sieo(),
        e_rioo(), e_tmrt(), e_bmsg(),
        e_rdzo(), e_rofo(), e_rufo(), e_rieo(), e_sioe(), e_sdze(),
        e_sofe(), e_sufe(), e_siee(), e_rioe(), e_rdze(), e_rofe(),
        e_rufe(), e_riee(), e_sall(), e_rall(), e_save(), e_rest(),
        e_tmsg(), e_xall(), e_xarg(), e_xdbz(), e_xine(), e_xiob(),
        e_xiop(), e_xios(), e_xnor(), e_xofl(), e_xset(), e_xufl();
a_bool  e_io_e(), e_dz_e(), e_of_e(), e_uf_e(), e_ie_e(), e_io_o(),
        e_dz_o(), e_of_o(), e_uf_o(), e_ie_o();
#endif





