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

/* CVS $Id: t_il10.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_il10.c                              */
/*                                                              */
/*      Interval version of log10                               */
/*                                                              */
/****************************************************************/
 
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif
 
#define log10e      /* log10(e) */          \
           EXTREAL(0x3F, 0xFD,              \
                   0xDE, 0x5B, 0xD8, 0xA9,  \
                   0x37, 0x28, 0x71, 0x95)  /* +0.43429448190325182765 */ \
/* ->A9 D8 5B DE  95 71 28 37  AF AA 5B 35   INTEL Darst.!!      */
 
extern a_real *r_pinf;
extern a_real *r_zero;
extern a_real *r_one_;
extern a_real *r_eps_;
 
 
/*-----------------------------------------------------------------*
 | log10(x) = ln(x)*log10(e) gerichted rerundet                    |
 *-----------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int log10r(const ExtReal *arg, ExtReal *res, int rndg)
#else
static int log10r(arg, res, rndg)
const ExtReal   *arg;
      ExtReal   *res;
      int rnding;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int log10r(ExtReal *arg, ExtReal *res, int rndg)
#else
static int log10r(arg, res, rndg)
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
static char    logvalue[sizeof(ExtReal)] = log10e;
 
memcpy(&logcnst,logvalue,sizeof(logcnst));
#else
static const ExtReal logcnst = log10e;
#endif
#else
static ExtReal logcnst;
static char    logvalue[sizeof(ExtReal)] = log10e;
 
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
} /* log10r() */


#ifdef LINT_ARGS
a_intv t_il10(a_intv ai)
#else
a_intv t_il10(ai)
 
a_intv ai;
#endif
   {
   a_intv   res;
   int      rnd, rc, i;
   IExtReal  a, r;
   a_real half, maxx, min;
   a_intg expo, expo_max;
   a_real ten, p;

   E_SPUSH("t_il10")
 
   ((a_btyp *)&ten)[B_HPART] = 0x40240000;   /*    10   */
   ((a_btyp *)&ten)[B_LPART] = 0;
   ((a_btyp *)&min)[B_HPART] = 0xC090C800;   /* -1074   */
   ((a_btyp *)&min)[B_LPART] = 0;
   ((a_btyp *)&maxx)[B_HPART] = 0x40900000;   /*  1024   */
   ((a_btyp *)&maxx)[B_LPART] = 0;
   ((a_btyp *)&half)[B_HPART] = 0x3FE00000;  /* 0.5     */
   ((a_btyp *)&half)[B_LPART] = 0;
 
   if (r_eq(ai.INF, *r_one_) == TRUE)
      {
      R_ASSIGN(res.INF, *r_zero); /* ai.INF = 1 = 10**0 */
      }
   else
      {
      if (r_eq(r_frac(ai.INF), *r_zero) == TRUE)
        {
        expo_max = 23;
        R_ASSIGN(p, *r_one_);
        for (expo = 0; expo < 23; expo++)
           {
           R_ASSIGN(p, r_muld(p, ten));
           if ( r_eq(p, ai.INF)==TRUE )
              {
              R_ASSIGN(res.INF, r_flot(expo)); /* ai.INF = 10**expo */
              goto ubound;
              }
           } /* for i, ... */
        }
      rnd = getrndmode();
      longreal_to_extreal((LongReal *)&ai.INF, &a.l);
 
      if ((rc = log10r(&a.l, &r.l, DOWN))!=0)
         ieee_aborti1(rc, &ai);
 
      setrndmode(DOWN);
      if ((rc = extreal_to_longreal(&r.l, (LongReal *)&res.INF))!=0)
        ieee_aborti1(rc, &ai);
      setrndmode(rnd);  /* restore rounding mode */
   }    /* at this point: lower bound of result has been computed */
 
ubound:
   if (r_eq(ai.SUP, *r_one_) == TRUE)
      {
      R_ASSIGN(res.SUP, *r_zero); /* ai.SUP = 1 = 10**0 */
      }
   else
      {
      if (r_eq(r_frac(ai.SUP), *r_zero) == TRUE)
        {
        expo = 0;
        R_ASSIGN(p, *r_one_);
        for (i=0; i < expo_max; i++)
           {
           expo++;
           R_ASSIGN(p, r_muld(p, ten));
           if ( r_eq(p, ai.SUP)==TRUE )
              {
              R_ASSIGN(res.SUP, r_flot(expo)); /* ai.SUP = 10**expo */
              goto end;
              }
           } /* for i, ... */
        }
      rnd = getrndmode();
      longreal_to_extreal((LongReal *)&ai.SUP, &a.l);
 
      if ((rc = log10r(&a.l, &r.l, UP))!=0)
         ieee_aborti1(rc, &ai);
 
      setrndmode(UP);
      if ((rc = extreal_to_longreal(&r.l, (LongReal *)&res.SUP))!=0)
        ieee_aborti1(rc, &ai);
      setrndmode(rnd);  /* restore rounding mode */
   }    /* both bounds are computed */
 
end:
      E_SPOPP("t_il10")
   return res;
   }
 
/* ------------------------------------------------------------ */
 





