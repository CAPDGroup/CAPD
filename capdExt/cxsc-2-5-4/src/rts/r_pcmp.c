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

/* CVS $Id: r_pcmp.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_pcmp.c                              */
/*                                                              */
/*      Entry           : a_intg r_pcmp(a,b,c,d)                */
/*                        a_real a,b,c,d;                       */
/*                                                              */
/*      Arguments       : a = operand of first product          */
/*                        b = operand of first product          */
/*                        c = operand of second product         */
/*                        d = operand of second product         */
/*                                                              */
/*      Function value  : -1  if  a*b<c*d                       */
/*                         0  if  a*b=c*d                       */
/*                         1  if  a*b>c*d                       */
/*                                                              */
/*      Description     : Compare products of real numbers      */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

#ifdef LINT_ARGS
local a_intg r_pcmp(a_real a,a_real b,a_real c,a_real d)
#else
local a_intg r_pcmp(a,b,c,d)

a_real a;
a_real b;
a_real c;
a_real d;
#endif
        {
        a_bool za,zb,zc,zd,va,vb,vc,vd;
        a_btyp ma[D_U_RATIO],mb[D_U_RATIO],
               mc[D_U_RATIO],md[D_U_RATIO];
        a_intg ea,eb,ec,ed;
        a_btyp p[4*D_U_RATIO]; 
        int i;   

        za = b_deko(a,&ea,ma,&va);
        zb = b_deko(b,&eb,mb,&vb);
        zc = b_deko(c,&ec,mc,&vc);
        zd = b_deko(d,&ed,md,&vd);

        va ^= vb;
        vc ^= vd;
         
        if (ea>EXPO_MAX || eb>EXPO_MAX || ec>EXPO_MAX || ed>EXPO_MAX)
           {
           if ((ea>EXPO_MAX && SIGNALING(ma[0])) ||
               (eb>EXPO_MAX && SIGNALING(mb[0])) ||
               (ec>EXPO_MAX && SIGNALING(mc[0])) ||
               (ed>EXPO_MAX && SIGNALING(md[0])))
              {
              e_trap(INV_OP+E_IEEE,10,E_TMSG,5,
                     E_TDBL,&a,E_TDBL,&b,E_TDBL,&c,E_TDBL,&d);
              }
           else if ((ea>EXPO_MAX && !MANT_INFINITY(ma)) ||
               (eb>EXPO_MAX && !MANT_INFINITY(mb)) ||
               (ec>EXPO_MAX && !MANT_INFINITY(mc)) ||
               (ed>EXPO_MAX && !MANT_INFINITY(md)))
              {
              e_trap(INV_ARG,10,E_TMSG,14,
                     E_TDBL,&a,E_TDBL,&b,E_TDBL,&c,E_TDBL,&d);
              }
           else if ((ea>EXPO_MAX && zb) ||
               (eb>EXPO_MAX && za) ||
               (ec>EXPO_MAX && zd) ||
               (ed>EXPO_MAX && zc))
              {
              e_trap(INV_OP+E_IEEE,10,E_TMSG,10,
                     E_TDBL,&a,E_TDBL,&b,E_TDBL,&c,E_TDBL,&d);
              }
           else if (ea>EXPO_MAX || eb>EXPO_MAX)
              {
              if (va==vc && (ec>EXPO_MAX || ed>EXPO_MAX)) 
                 return(0);
              return( (va) ? -1 : 1 );
              }
           else
              {
              return( (vc) ? 1 : -1);
              }
           return(0);   
           }

        if (za || zb) return( (zc || zd) ? 0 : (vc) ? 1 : -1);
        if (zc || zd) return( (va) ? -1 : 1);
        if (va!=vc) return( (va) ? -1 : 1);

        if (ea==EXPO_MIN && !(ma[0] & HIDDEN_BIT))
           {
           do {
              b_shl1(ma,D_U_RATIO);
              ea--;
              }
           while (!(ma[0] & HIDDEN_BIT));
           }
        if (eb==EXPO_MIN && !(mb[0] & HIDDEN_BIT))
           {
           do {
              b_shl1(mb,D_U_RATIO);
              eb--;
              }
           while (!(mb[0] & HIDDEN_BIT));
           }
        if (ec==EXPO_MIN && !(mc[0] & HIDDEN_BIT))
           {
           do {
              b_shl1(mc,D_U_RATIO);
              ec--;
              }
           while (!(mc[0] & HIDDEN_BIT));
           }
        if (ed==EXPO_MIN && !(md[0] & HIDDEN_BIT))
           {
           do {
              b_shl1(md,D_U_RATIO);
              ed--;
              }
           while (!(md[0] & HIDDEN_BIT));
           }

        ea += eb;
        ec += ed;

        if (ea>ec+1) return( (va) ? -1 : 1);
        if (ea+1<ec) return( (va) ? 1 : -1);

        if (ea<ec) b_shl1(mc,D_U_RATIO);
        else if (ea>ec) b_shl1(ma,D_U_RATIO);

        B_CLEAR(p,4*D_U_RATIO)

        b_prod(ma,mb,p);
        b_prod(mc,md,p+4);

        for (i=0;i<2*D_U_RATIO;i++)
           if (p[i]<p[i+2*D_U_RATIO]) 
              return( (va) ? 1 : -1);
           else if (p[i]>p[i+2*D_U_RATIO]) 
              return( (va) ? -1 : 1);

        return(0);
        }





