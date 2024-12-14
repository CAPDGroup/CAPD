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

/* CVS $Id: t_iex2.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_iex2.c                              */
/*                                                              */
/*      Interval version of 2**x                                */
/*                                                              */
/****************************************************************/
 
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif
extern a_real *r_pinf;
extern a_real *r_zero;
extern a_real *r_eps_;
 
/*--------------------------------------------------------------*
 | exp2rd   exp2 with rounding mode                             |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int exp2rd(const ExtReal *arg, ExtReal *res, int rnd)
#else
static int exp2rd(arg, res, rnd)
const ExtReal *arg;
      ExtReal *res;
      int rnd;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int exp2rd(ExtReal *arg, ExtReal *res, int rnd)
#else
static int exp2rd(arg, res, rnd)
ExtReal *arg;
ExtReal *res;
int rnd;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal     r;          /* Ergebnis vor Rundung              */
   int         ret;        /* Rueckgabe                         */
   int         check;      /* Rueckgabe von Makro ArgCheck      */
 
   /* --- pruefe Argument, dann Pruefung aus --- */
   ArgCheck1(Exp, arg, res);
   /* arg_check = Off; */
 
   /* --- exp --- */
   ret = t_2exp(arg, &r);
 
   /* --- Abbruch bei Fehler --- */
   if(ret!=NoErr) {
      copyee(&r, res);
      arg_check = On;
      return (ret);
   }
 
   /* nicht runden, falls Argument Null ist */
   if (!cmpee(arg, &Zero)) ret |= LB_Exact;
 
   /* --- Rundungs-Fehler --- */
   if ( rnd == UP )
      {
      if(!(ret&LB_Exact)) round_rel(UP,  &r, &EpsExp, res);
      else copyee(&r, res);
      }
   else
      {
      if(!(ret&LB_Exact)) round_rel(DOWN,  &r, &EpsExp, res);
      else copyee(&r, res);
      }
   /* Verbesserung des Ergebnis' in Spezialfaellen */
#if INT_HPREC
   { int vz;
     vz = SGNE(arg);
     if (vz==POS&&-1==cmpee(res, &One)) /* exp2 x >= 1 */
        copyee(&One, res);
   }
#endif /* INT_HPREC */
 
   /* Unterlaufbehandlung */
   if(0==cmpee(res,&Zero))
      copyee(&LongRealDenormMin,res);
   if (rnd == UP)
      {
      /* --- ArgPruefung wieder an --- */
      arg_check = On;
      }
 
   /* --- kein Fehler mehr moeglich --- */
   return NoErr;
} /* exp2rd() */


#ifdef LINT_ARGS
a_intv t_iex2(a_intv ai)
#else
a_intv t_iex2(ai)
 
a_intv ai;
#endif
   {
   a_intv   res;
   int      rnd, rc;
   IExtReal  a, r;
   a_real half, maxx, min;
   a_intg expo;

   E_SPUSH("t_iex2")
 
   ((a_btyp *)&min)[B_HPART] = 0xC090C800;   /* -1074   */
   ((a_btyp *)&min)[B_LPART] = 0;
   ((a_btyp *)&maxx)[B_HPART] = 0x40900000;   /*  1024   */
   ((a_btyp *)&maxx)[B_LPART] = 0;
   ((a_btyp *)&half)[B_HPART] = 0x3FE00000;  /* 0.5     */
   ((a_btyp *)&half)[B_LPART] = 0;
 
   if (r_ge(ai.INF, maxx) == TRUE)
      {
      e_trap(INV_ARG, 6, E_TMSG, 48,
             E_TDBL+E_TEXT(5), &ai.INF, E_TDBL+E_TEXT(6), &ai.SUP);
      R_ASSIGN(res.INF, *r_pinf);
      R_ASSIGN(res.SUP, *r_pinf);
      }
   else if (r_le(ai.INF, min) == TRUE)
      {
      R_ASSIGN(res.INF, *r_zero); /* underflow */
      }
   else if (r_sign(r_frac(ai.INF)) == 0)
      {
      expo = r_trun(ai.INF);
      R_ASSIGN(res.INF,r_comp(half, expo+1)); /* exact result */
      }
   else
      {
      rnd = getrndmode();
      longreal_to_extreal((LongReal *)&ai.INF, &a.l);
 
      if ((rc = exp2rd(&a.l, &r.l, DOWN))!=0)
         ieee_aborti1(rc, &ai);
 
      setrndmode(DOWN);
      if ((rc = extreal_to_longreal(&r.l, (LongReal *)&res.INF))!=0)
        ieee_aborti1(rc, &ai);
      setrndmode(rnd);  /* restore rounding mode */
      }    /* at this point: lower bound of result has been computed */
 
   if (r_ge(ai.SUP, maxx) == TRUE)   /* consider upper bound */
      {
      e_trap(INV_ARG, 6, E_TMSG, 48,
             E_TDBL+E_TEXT(5), &ai.INF, E_TDBL+E_TEXT(6), &ai.SUP);
      R_ASSIGN(res.INF, *r_pinf);
      R_ASSIGN(res.SUP, *r_pinf);
      }
   else if (r_le(ai.SUP, min) == TRUE)
      {
      R_ASSIGN(res.INF, *r_eps_); /* underflow */
      }
   else if (r_sign(r_frac(ai.SUP)) == 0)
      {
      expo = r_trun(ai.SUP);
      R_ASSIGN(res.SUP,r_comp(half, expo+1));
      }
   else
      {
      rnd = getrndmode();
      longreal_to_extreal((LongReal *)&ai.SUP, &a.u);
 
      if ((rc = exp2rd(&a.u, &r.u, UP))!=0)
         ieee_aborti1(rc, &ai);
 
      setrndmode(UP);
      if ((rc = extreal_to_longreal(&r.u, (LongReal *)&res.SUP))!=0)
         ieee_aborti1(rc, &ai);
      setrndmode(rnd);
      }  /* upper bound */
 
      E_SPOPP("t_iex2")
   return res;
   }
 
/* ------------------------------------------------------------ */
 





