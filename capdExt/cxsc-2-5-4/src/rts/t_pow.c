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

/* CVS $Id: t_pow.c,v 1.22 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_pow.c                               */
/*                                                              */
/*      Entries         : a_real t_pow(bas,exp)                 */
/*                        a_real bas;                           */
/*                        a_real exp;                           */
/*                                                              */
/*      Arguments       : bas  = base of power                  */
/*                        exp  = exponent of power              */
/*                                                              */
/*      Description     : Power function.                       */
/*                                                              */
/*                   problems on i.d. -32^0.2 !!                */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

#define OUT(exc) return exc_handle_2(Pow, exc, bas, exp, res)

/* StdFctReal2(t_pow,powee) */
#ifdef LINT_ARGS
a_real t_pow(a_real arg1, a_real arg2)
#else
a_real t_pow(arg1,arg2)

a_real arg1;
a_real arg2;
#endif
        {
        a_real   res;
        int      rnd, rc;
        ExtReal  a1, a2, r;

        E_SPUSH("t_pow")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg1, &a1);
        longreal_to_extreal((LongReal *)&arg2, &a2);

        if ((rc = powee(&a1, &a2, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
             ieee_abortr2(rc, &arg1, &arg2);

        setrndmode(rnd);

        E_SPOPP("t_pow")
        return res;
        }
/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | powee                                                        |
 | res = expee(exp*lnee(bas))                                   |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int powee(const ExtReal *bas, const ExtReal *exp, ExtReal *res)
#else
int powee(bas, exp, res)
const ExtReal  *bas;
const ExtReal  *exp;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int powee(ExtReal *bas, ExtReal *exp, ExtReal *res)
#else
int powee(bas, exp, res)
ExtReal  *bas;
ExtReal  *exp;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal     rr;         /* RundungsInfo hier dummy           */
   int         rnd;        /* RundungsMode                      */
   int         ret;        /* Rueckgabe                         */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

   rnd = getrndmode();

   /* --- pruefe Argument --- */
   ArgCheck2(Pow, bas, exp, res);

   ret = powsub(bas, exp, res, &rr);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   return ret;
} /* powee() */
/*--------------------------------------------------------------*
 | powee                                                        |
 | res = expee(exp*lnee(bas))                                   |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int powsub(const ExtReal *bas, const ExtReal *exp, ExtReal *res, ExtReal *rr)
#else
int powsub(bas, exp, res, rr)
const ExtReal  *bas;
const ExtReal  *exp;
      ExtReal  *res;
      ExtReal  *rr;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int powsub(ExtReal *bas, ExtReal *exp, ExtReal *res, ExtReal *rr)
#else
int powsub(bas, exp, res, rr)
ExtReal  *bas;
ExtReal  *exp;
ExtReal  *res;
ExtReal  *rr;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal     r;          /* Ergebnis ln                       */
   ExtReal     b;          /* |bas|                             */
   ExtReal     e;          /* Ergebnis Exponent                 */
   ExtReal     i;          /* int (exp)                         */
   int         esgn;       /* Vorzeichen des Exponenten         */
   int         eint;       /* == 0 falls Exponent ganzzahlig    */
   int         e_even;     /* == 0 falls Exponent Geradzahlig   */
   int         bsgn;       /* Vorzeichen der Basis              */
   int         rnd;        /* RundungsMode                      */
   int         ret;        /* Rueckgabe                         */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

   /* --- pruefe Argument --- */
   ArgCheck2(Pow, bas, exp, res);

   rndintee(exp, &i);           /* i = int(exp)              */
   eint=cmpee(exp, &i);         /* eint = 0: exp Ganzzahlig  */
   esgn=cmpee(exp, &Zero);      /* esgn = Vorzeichen Exp     */
   e_even=mod2e(exp);           /* e_even = 0 falls gerade   */
   
   if(1!=(bsgn=cmpee(bas, &Zero)))      /* bas <=0                   */
   {     
       if(0!=eint || POS!=esgn)         /* exp<=0 OR exp not in Z    */
       {
         if(0==bsgn) 
	    {                              /* bas = 0                   */
            if(POS!=esgn)
               OUT(DOMAIN);             /* raus, falls exp <=0       */
         }
         else                          /* bas < 0                   */
           if(0!=eint)                 /* exp in R-Z                 */
		 {                           /* if exp is odd root, no error*/
		    divee(&One,exp,&b);      /* note:exp <>0 , numerical   */
							    /*      problems awaited !    */
              rndintee(&b,&e);         /* e = int(b)                 */
              if ((mod2e(&e)!=0) && (0==cmpee(&b,&e)));       
		    else OUT(DOMAIN);
           }
       }
   }

   /* --- bei Null sofort fertig --- */
   if (0==cmpee(bas, &Zero))
   {
	 copyee (&One, rr);
      return copyee(&Zero, res);
   }

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* --- b := |bas| --- */
   bsgn = SGNE(bas);
   absee(bas, &b);

   /* --- pow --- */
   if(NoErr==(ret=lnee(&b, &r))) {
      mulee(&r, exp, &e);
      if(1==cmpee(&e, &MaxArgExp)) return OVER_FLOW;
      ret = expee(&e, res);

      if(/*esgn==NEG &&*/ eint==0)      /* exp < 0 und ganzzahlig */
         if (!(bsgn==POS || e_even==0)) /* falls nicht bas>0      */
                                        /* oder exp gerade        */
            chsee(res, res);            /* Vorzeichenwechsel      */
   }

   /* --- Rundungsinfo --- */
   copyee(&e, rr);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   return ret;
} /* powsub() */





