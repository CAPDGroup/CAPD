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

/* CVS $Id: t_ilg2.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_ilg2.c                              */
/*                                                              */
/*      Interval version of 2**x                                */
/*                                                              */
/****************************************************************/
 
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif
 
#define log2e      /* log2(e) */           \
           EXTREAL(0x3f, 0xff,              \
                   0xB8, 0xAA, 0x3B, 0x29,  \
                   0x5C, 0x17, 0xF0, 0xBB)  /* +1.44269504088896340735*/\
 
/* ->29 3B AA B8   BB F0 17 5C   D0 FE 87 BE INTEL Darst>!!*/
 
extern a_real *r_pinf;
extern a_real *r_zero;
extern a_real *r_eps_;
 
 
/*-----------------------------------------------------------------*
 | log2(x) = ln(x)*log2(e) gerichted rerundet                      |
 *-----------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int log2rd(const ExtReal *arg, ExtReal *res, int rndg)
#else
static int log2rd(arg, res, rndg)
const ExtReal   *arg;
      ExtReal   *res;
      int rnding;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int log2rd(ExtReal *arg, ExtReal *res, int rndg)
#else
static int log2rd(arg, res, rndg)
ExtReal   *arg;
ExtReal   *res;
int      rndg;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int         rnd;        /* RundungsMode                       */
   int         rr;         /* Rueckgabe RundungsArt, hier dummy  */
   int         ret;        /* Rueckgabe                          */
   int         check;      /* Rueckgabe von Makro ArgCheck       */
#ifdef ANSI_C
#if SUN4_CPP_C
static ExtReal logcnst;
static char    logvalue[sizeof(ExtReal)] = log2e;
 
memcpy(&logcnst,logvalue,sizeof(logcnst));
#else
static const ExtReal logcnst = log2e;
#endif
#else
static ExtReal logcnst;
static char    logvalue[sizeof(ExtReal)] = log2e;
 
memcpy(&logcnst,logvalue,sizeof(logcnst));
#endif
   /* --- pruefe Argument --- */
   ArgCheck1(Ln, arg, res);
 
   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
 
   setrndmode(rndg);
 
   /* --- Log2 gerichtet gerundet --- */
   ret = _s_ln(arg, res, &rr);
   mulee(res, &logcnst, res);
   if ((rndg == DOWN))
      {
      round_ln(DOWN, rr, res, res);
      }
   else
      {
      round_ln(UP, rr, res, res);
      }
 
   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);
 
   return ret;
} /* log2rd() */


#ifdef LINT_ARGS
a_intv t_ilg2(a_intv ai)
#else
a_intv t_ilg2(ai)
 
a_intv ai;
#endif
   {
   a_intv   res;
   int      rnd, rc;
   IExtReal  a, r;
   a_real half, maxx, min;
   a_intg expo;
 
   E_SPUSH("t_ilg2")

   ((a_btyp *)&min)[B_HPART] = 0xC090C800;   /* -1074   */
   ((a_btyp *)&min)[B_LPART] = 0;
   ((a_btyp *)&maxx)[B_HPART] = 0x40900000;   /*  1024   */
   ((a_btyp *)&maxx)[B_LPART] = 0;
   ((a_btyp *)&half)[B_HPART] = 0x3FE00000;  /* 0.5     */
   ((a_btyp *)&half)[B_LPART] = 0;
 
   expo = r_expo(ai.INF);
   if (r_eq(r_comp(half, expo), ai.INF)==TRUE) /* exact result */
      {
      R_ASSIGN(res.INF, r_flot(expo-1)); /* ai.INF = 2**(expo-1) */
      }
   else
      {
      rnd = getrndmode();
      longreal_to_extreal((LongReal *)&ai.INF, &a.l);
 
      if ((rc = log2rd(&a.l, &r.l, DOWN))!=0)
         ieee_aborti1(rc, &ai);
 
      setrndmode(DOWN);
      if ((rc = extreal_to_longreal(&r.l, (LongReal *)&res.INF))!=0)
        ieee_aborti1(rc, &ai);
      setrndmode(rnd);  /* restore rounding mode */
      }    /* at this point: lower bound of result has been computed */
 
   /* consider upper bound */
   expo = r_expo(ai.SUP);
   if (r_eq(r_comp(half, expo), ai.SUP)==TRUE) /* exact result */
      {
      R_ASSIGN(res.SUP, r_flot(expo-1)); /* ai.SUP = 2**(expo-1) */
      }
   else
      {
      rnd = getrndmode();
      longreal_to_extreal((LongReal *)&ai.SUP, &a.u);
 
      if ((rc = log2rd(&a.u, &r.u, UP))!=0)
         ieee_aborti1(rc, &ai);
 
      setrndmode(UP);
      if ((rc = extreal_to_longreal(&r.u, (LongReal *)&res.SUP))!=0)
        ieee_aborti1(rc, &ai);
      setrndmode(rnd);  /* restore rounding mode */
      }  /* upper bound */
 
 
      E_SPOPP("t_ilg2")
   return res;
   }
 
/* ------------------------------------------------------------ */
 





