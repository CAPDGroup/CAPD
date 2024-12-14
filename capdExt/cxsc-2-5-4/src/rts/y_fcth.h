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

/* CVS $Id: y_fcth.h,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : y_fcth.h                              */
/*                                                              */
/*      Description     : function prototypes                   */
/*                                                              */
/*                   prototypes for variable length arg. lists  */
/****************************************************************/

_PROTOTYPE(void     y_asgn, (y_dscp d,y_dscp s));
_PROTOTYPE(void     y_free, (y_dscp d));
_PROTOTYPE(void     y_init, (y_dscp d,a_byte dim,size_t elsize));
_PROTOTYPE(void     y_inid, (y_dscp d,a_byte dim,size_t elsize));
_PROTOTYPE(y_dscp   y_alck, (y_dscp d));

_PROTOTYPE(a_btyp   y_ixch, (a_intg index,y_bnds range));
_PROTOTYPE(void     y_temp, (y_dscp d));
_PROTOTYPE(void     y_utmp, (y_dscp d));
_PROTOTYPE(void     y_vlcp, (y_dscp d));
_PROTOTYPE(a_btyp   y_yxch, (a_intg index,y_bnds range));

#ifdef LINT_ARGS
#if C_P_7
a_VOID   y_inxc(y_dscp d,...);
a_VOID   y_inxn(y_dscp d,...);
a_VOID   y_suba(a_VOID m,a_VOID s,a_char *mode,...);
a_VOID   y_stat(y_dscp d,a_VOID s,size_t z,a_byte dim,...);
a_VOID   y_yxcn(y_dscp d,...);
a_VOID   y_ynxn(y_dscp d,...);
void     y_new (y_dscp d,...);
#else
a_VOID   y_inxc();
a_VOID   y_inxn();
a_VOID   y_suba();
a_VOID   y_stat();
a_VOID   y_yxcn();
a_VOID   y_ynxn();
void     y_new ();
#endif
#else
a_VOID   y_inxc(), y_inxn(), y_stat(), y_suba(), y_ynxn(), y_yxcn();
void     y_new();
#endif





