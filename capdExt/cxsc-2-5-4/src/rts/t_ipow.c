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

/* CVS $Id: t_ipow.c,v 1.22 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_ipow.c                              */
/*                                                              */
/*      Synopsis        : t_ipow(ai,bi  )                       */
/*                        a_intv ai;                            */
/*                        a_intv bi;                            */
/*                                                              */
/*      Description     : ai^bi                                 */
/*                                                              */
/*                                                              */
/*                   protoype for minmax()                      */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* StdFctInterval2(t_ipow,ipowee) */
#ifdef LINT_ARGS
a_intv t_ipow(a_intv ai, a_intv bi)
#else
a_intv t_ipow(ai,bi)
a_intv ai;
a_intv bi;
#endif
        {
        a_intv   res;
        int      rnd, rc;
        IExtReal  a, b, r;

        E_SPUSH("t_ipow")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&ai.INF, &a.l);
        longreal_to_extreal((LongReal *)&bi.INF, &b.l);
        longreal_to_extreal((LongReal *)&ai.SUP, &a.u);
        longreal_to_extreal((LongReal *)&bi.SUP, &b.u);

        if ((rc = ipowee(&a, &b, &r))!=0)
             ieee_aborti2(rc, &ai, &bi);

        setrndmode(DOWN);
        if ((rc = extreal_to_longreal(&r.l, (LongReal *)&res.INF))!=0)
             ieee_aborti2(rc, &ai, &bi);
        setrndmode(UP);
        if ((rc = extreal_to_longreal(&r.u, (LongReal *)&res.SUP))!=0)
             ieee_aborti2(rc, &ai, &bi);
        setrndmode(rnd);

        E_SPOPP("t_ipow")
        return res;
        }

/* ------------------------------------------------------------ */

#ifdef LINT_ARGS
#ifdef ANSI_C
static int typ1(const IExtReal *bas, const IExtReal *exp,
                IExtReal *res, IExtReal *ri);
static int typ2(const IExtReal *bas, const  ExtReal *ei,
                IExtReal *res, IExtReal *ri);
static int minmax(int rnd, const IExtReal *bas, const IExtReal *exp,
                ExtReal *res, ExtReal *ri);
#else  /* NOT ANSI_C */
static int typ1(IExtReal *bas, IExtReal *exp, IExtReal *res, IExtReal *ri);
static int typ2(IExtReal *bas,  ExtReal *ei,  IExtReal *res, IExtReal *ri);
static int minmax(int rnd, IExtReal *bas, IExtReal *exp,
                  ExtReal *res, ExtReal *ri);
#endif /* ANSI_C */
#else  /* NOT LINT_ARGS */
static int typ1(), typ2(), minmax();
#endif /* LINT_ARGS */

/*--------------------------------------------------------------*
 | ipowee                                                       |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int ipowee (const IExtReal *bas, const IExtReal *exp, IExtReal *res)
#else
int ipowee (bas, exp, res)
const IExtReal *bas;
const IExtReal *exp;
      IExtReal *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int ipowee (IExtReal *bas, IExtReal *exp, IExtReal *res)
#else
int ipowee (bas, exp, res)
IExtReal *bas;
IExtReal *exp;
IExtReal *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   IExtReal    r;          /* Ergebnis vor Rundung              */
   IExtReal    ri;         /* RundungsInfo                      */
   ExtReal     eps;        /* RundungsFehler                    */
   ExtReal     ei;         /* int von exp->u                    */
   int         ret;        /* Rueckgabe                         */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

   /* --- pruefe Argument, dann Pruefung aus --- */
   ArgCheckI2(IPow, bas, exp, res);
   arg_check = Off;

   /* --- ganzzahlige Anteile --- */
   rndintee(&(exp->u), &ei);

   /* --- TypenUnterscheidung --- */
   /* e.sup ==e.inf == flor (e.sup) */
   if(0==cmpee(&(exp->l), &ei) && 0==cmpee(&(exp->u), &ei))    
         ret = typ2(bas, &ei, &r, &ri); /* ret=bas^integer */
    else ret = typ1(bas, exp, &r, &ri); /* ret=bas^interval*/

   /* --- Abbruch bei Fehler --- */
   if(ret!=NoErr) {
      icopyee(&r, res);
      arg_check = On;
      return ret;
   }


   /* --- Rundungs-Fehler --- */
   mulee(&EpsPow, &ri.u, &eps);
   round_rel(UP,   &r.u, &eps, &(res->u));
   mulee(&EpsPow, &ri.l, &eps);
   round_rel(DOWN, &r.l, &eps, &(res->l));

   /* --- ArgPruefung wieder an --- */
   arg_check = On;

   return ret;
} /* ipowee() */

/*--------------------------------------------------------------*
 | typ1                                                         |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int typ1(const IExtReal *bas, const IExtReal *exp,
                IExtReal *res, IExtReal *ri)
#else
static int typ1(bas, exp, res, ri)
const IExtReal *bas;
const IExtReal *exp;
      IExtReal *res;
      IExtReal *ri;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int typ1(IExtReal *bas, IExtReal *exp, IExtReal *res, IExtReal *ri)
#else
static int typ1(bas, exp, res, ri)
IExtReal *bas;
IExtReal *exp;
IExtReal *res;
IExtReal *ri;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
   const ExtReal *basu;    /* Argumente PunktPower              */
   const ExtReal *basl;    /*                                   */
   const ExtReal *expu;    /*                                   */
   const ExtReal *expl;    /*                                   */
#else
   ExtReal *basu;
   ExtReal *basl;
   ExtReal *expu;
   ExtReal *expl;
#endif
   int         cmpexp;     /* Ergebnis Vergleich Basis          */
   int         cmpbas;     /* Ergebnis Vergleich Basis          */
   int         retu;       /* Rueckgabe                         */
   int         retl;       /* Rueckgabe                         */

   /* ******************************************** 1. Runde *** */
   /* --- bas->u <= 1 --- */
   if (1!=cmpee(&(bas->u), &One)) {
      cmpbas = -1; 
      expu = &(exp->l);
      expl = &(exp->u);
   }
   else
   /* --- bas->l >= 1 --- */
   if (-1!=cmpee(&(bas->l), &One)) {
      cmpbas = 1; /* bas->l >= 1 */
      expu = &(exp->u);
      expl = &(exp->l);
   }
   else
   /* --- bas->l < 1 < bas->u --- */
      cmpbas = 0;

   /* ******************************************* 2. Runde *** */
   /* --- exp->u <= 0 --- */
   if (1!=cmpee(&(exp->u), &Zero)) {
      cmpexp = -1; 
      if (1!=cmpabsee(&(bas->l), &(bas->u))) 
	 { /* |l|<=|u| */
         basu = &(bas->l);
         basl = &(bas->u);
	 }
	 else {
         basu = &(bas->u);
         basl = &(bas->l);
	 }
   }
   else
   /* --- exp->l >= 0 --- */
   if (-1!=cmpee(&(exp->l), &Zero)) {
      cmpexp = 1; 
      basu = &(bas->u);
      basl = &(bas->l);
   }
   /* --- exp->l < 0 < exp->u --- */
   else {
      cmpexp = 0;
      switch(cmpbas) {
      case -1:/* bas <= 1 */
         basu = &(bas->l);
         basl = &(bas->l);
         break;
      case  0: /* 1 in bas */
         break;
      case  1: /* bas >=1 */
         basu = &(bas->u);
         basl = &(bas->u);
         break;
      } /* switch */
   }

   /* ********************************************** 3. Runde *** */
   if (cmpbas == 0) /* 1 in bas */
      switch(cmpexp) {
      case -1: /* exp <=0 */
         expu = &(exp->l);
         expl = &(exp->l);
         break;
      case  0: /* 0 in exp */
         break;
      case  1: /* exp >=0 */
         expu = &(exp->u);
         expl = &(exp->u);
         break;
      } /* switch */

/* --- Aufruf PunktFunktion Power -----------------------------*/
   /* --- Faelle mit DoppelAuswertung --- */
   if (cmpexp == 0 && cmpbas == 0) {
      retu = minmax(UP,   bas, exp, &(res->u), &(ri->u));
      retl = minmax(DOWN, bas, exp, &(res->l), &(ri->l));
   }
   /* --- alle anderen Faelle --- */
   else {
      retu = powsub(basu, expu, &(res->u), &(ri->u));
      retl = powsub(basl, expl, &(res->l), &(ri->l));
   }

   return max(retu, retl);
} /* typ1() */

/*--------------------------------------------------------------*
 | typ2                                                         |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int typ2(const IExtReal *bas, const ExtReal *ei,
                IExtReal *res, IExtReal *ri)
#else
static int typ2(bas, ei, res, ri)
const IExtReal *bas;
const  ExtReal *ei;  /* Integer Exponent */
      IExtReal *res;
      IExtReal *ri;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int typ2(IExtReal *bas, ExtReal *ei, IExtReal *res, IExtReal *ri)
#else
static int typ2(bas, ei, res, ri)
IExtReal *bas;
 ExtReal *ei;  /* Integer Exponent */
IExtReal *res;
IExtReal *ri;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
   const ExtReal *basu;    /* Argumente PunktPower              */
   const ExtReal *basl;    /*                                   */
#else
   ExtReal *basu;
   ExtReal *basl;
#endif
   ExtReal     ir, r;      /* Ganzzahlige Argumente             */
   ExtReal     abu;        /* |bas->u|                          */
   ExtReal     abl;        /* |bas->l|                          */
   int         retu;       /* Rueckgabe                         */
   int         retl;       /* Rueckgabe                         */
   int         flag  ;     /* bas contains zero                 */
					  /* 0 = no zero, 1=zero, 2=neg. exp   */
   mulee(ei, &Half, &r);
   rndintee(&r, &ir);

   flag   = (( 1!=cmpee(&(bas->l), &Zero)) && /* 0 in [l,u]? */
             (-1!=cmpee(&(bas->u), &Zero)))? 1:0;

   /* --- exp < 0 --- */
   if (-1==cmpee(ei, &Zero))
   {
      if (flag ==1) return DOMAIN; /* error in argument bas */
	 else flag=2;
   }

   /* --- exp gerade --- */
   if ( 0==cmpee(&r, &ir)) {
      absee(&(bas->u), &abu);
      absee(&(bas->l), &abl);
      if (1==cmpee(&abu, &abl)) /* abu > abl */
	 {
	    if (flag==2)
	    {/* exp<0, 0 not in [bas] */
	       basl = &abu;
	       basu = &abl;
	    }
	    else
	    {
	       basl = (flag==1) ? &Zero:&abl;
	       basu = &abu;
	    }
	 }
      else
	 {
	    if (flag==2)
	    {/* exp<0, 0 not in [bas] */
	       basl = &abl;
	       basu = &abu;
	    }
	    else
	    {
	       basl = (flag==1) ? &Zero:&abu;
            basu = &abl;
	    }
	 }
   }
   /* --- exp ungerade --- */
   else 
   {
	    if (flag==2)
	    {/* exp<0, 0 not in [bas] */
            basu = &(bas->l);
            basl = &(bas->u);
	    }
	    else
	    {
            basu = &(bas->u);
            basl = &(bas->l);
	    }
   }

   /* --- PunktFunktion --- */
   retu = powsub(basu, ei, &(res->u), &(ri->u));
   retl = powsub(basl, ei, &(res->l), &(ri->l));

   return max(retu, retl);
} /* typ2() */

/*--------------------------------------------------------------*
 | DoppelAuswertung                                             |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int minmax(int rnd, const IExtReal *bas, const IExtReal *exp,
                  ExtReal *res, ExtReal *ri)
#else
static int minmax(rnd, bas, exp, res, ri)
      int rnd;
const IExtReal *bas;
const IExtReal *exp;
      ExtReal  *res;
      ExtReal  *ri;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int minmax(int rnd, IExtReal *bas, IExtReal *exp,
                  ExtReal *res, ExtReal *ri)
#else
static int minmax(rnd, bas, exp, res, ri)
int rnd;
IExtReal *bas;
IExtReal *exp;
ExtReal  *res;
ExtReal  *ri;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal     rlnl;       /* ZwischenErgebnis ln               */
   ExtReal     rlnu;       /* ZwischenErgebnis ln               */
   ExtReal     r1;         /* ZwischenErgebnis                  */
   ExtReal     r2;         /* ZwischenErgebnis                  */
   int         ret;        /* Rueckgabe                         */

   /* --- doppelte Auswertung --- */
   if(NoErr!=(ret=lnee(&(bas->l), &rlnl))) return ret;
   if(NoErr!=(ret=lnee(&(bas->u), &rlnu))) return ret;
   if(rnd==UP) {
      mulee(&rlnu, &(exp->u), &r1);
      mulee(&rlnl, &(exp->l), &r2);
   } else{ /* DOWN */
      mulee(&rlnl, &(exp->u), &r1);
      mulee(&rlnu, &(exp->l), &r2);
   }

   /* --- min bzw max --- */
   if (rnd==cmpee(&r1, &r2)) {
      if(1==cmpee(&r1, &MaxArgExp)) return OVER_FLOW;
      ret = expee(&r1, res);
      copyee(&r1, ri);
   }
   else {
      if(1==cmpee(&r2, &MaxArgExp)) return OVER_FLOW;
      ret = expee(&r2, res);
      copyee(&r2, ri);
   }

   return ret;
} /* minmax() */





