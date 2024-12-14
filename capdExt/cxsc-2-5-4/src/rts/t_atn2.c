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

/* CVS $Id: t_atn2.c,v 1.22 2014/01/30 17:24:15 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_atn2.c                              */
/*                                                              */
/*      Entries         : a_real t_atn2(x,y)                    */
/*                        a_real x,y;                           */
/*                                                              */
/*      Arguments       : x/y  = argument of arcus tangent      */
/*                                                              */
/*      Description     : Arcus tangent function with quotient  */
/*                        argument                              */
/*                                                              */
/*                   break after return deleted                 */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

#define OUT(exc) return exc_handle_2(Atan, exc, x, y, res)

/* StdFctReal(t_atn2,atanee2) */
#ifdef LINT_ARGS
a_real t_atn2(a_real x, a_real y)
#else
a_real t_atn2(x, y)

a_real x,y;
#endif
        {
        a_real   res;
        int      rnd, rc;
        ExtReal  a1, a2, r;

        E_SPUSH("t_atn2")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&x, &a1);
        longreal_to_extreal((LongReal *)&y, &a2);

        if ((rc = atanee2(&a1, &a2, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
             ieee_abortr2(rc, &x, &y);

        setrndmode(rnd);

        E_SPOPP("t_atn2")
        return res;
        }

/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | Punkt Arcus Tangens                                          |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int atanee2 (const ExtReal *x, const ExtReal *y, ExtReal *res)
#else
int atanee2 (x, y, res)
const ExtReal   *x;
const ExtReal   *y;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int atanee2 (ExtReal *x, ExtReal *y, ExtReal *res)
#else
int atanee2 (x,y, res)
ExtReal   *x;
ExtReal   *y;
ExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int         rnd;        /* RundungsMode                      */
   int         ret=0;      /* Rueckgabe                         */
   /* Initialisation to avoid warning 'may be used uninitialized in this function' */

   int         check;      /* Rueckgabe von Makro ArgCheck      */

   /* --- pruefe Argument --- */
   ArgCheck2(Atan, x, y, res);

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

  /**************************************************************
   * Check and Handle Range Error and Special Cases             *
   **************************************************************/


  /*----------------------------*
   | Check and Handle Case x=0  |
   *----------------------------*/
   if(0==cmpee(x, &Zero))         /* x   = 0                   */
   {
      switch (cmpee(y, &Zero))
	 {
	  case 0:                     /* y = 0 */
		OUT(DOMAIN);          /* 0/0 undefined        */
		/* break; */

       case 1:                      /* y>0 */
          ret=copyee(&Zero, res);   /* ri = arctan(0/|y|) = 0 */
		break;

       case -1:                    /* y<0 */
          ret=copyee(&Pi, res);   /* ri = arctan(0/-|y|) = PI */
		break;
	 }
	 return ret;
   }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/
  /*---------------------------------*
   | Check and Handle Case y=0, x<>0 |
   *---------------------------------*/

   if(0==cmpee(y, &Zero))         /* y   = 0                   */
   {
	 if (SGNE(x)==NEG)
         ret=copyee(&MinusPiHalf, res);  /* ri = sgn(x)*arcctg(0/|x|) = -Pi/2 */
	 else
         ret=copyee(&PiHalf, res);   /* ri = sgn(x)*arcctg(0/|x|) = Pi/2 */
      return ret;
   }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/
  /*-----------------------------*
   | Check and Handle Case x =y  |
   *-----------------------------*/

   {
   ExtReal     arg;        /* quotient of x/y, y/x resp.        */

   switch (cmpabsee(x, y)) 
   {     
    case 0: /* x = y */
      if(SGNE(y)==POS)
	 {
         ret = copyee(&PiQuart, res);      /* ri=sgn(x)*arctan(1) =+-pi/4  */
	 }
      else
	 {
         ret = copyee(&Pi, &arg);         /* ri=sgn(x)*arctan(-1)=+-3*pi/4*/
         ret = mulee(&arg,&ThreeQuart, res);
      }
	 if (SGNE(x)==NEG) ret = chsee(res,res);

      return ret;

    case 1:
      /* --- x>y -> argument = y/x --- */
      ret = divee(y,x, &arg);
	 break;

    case -1:
      /* --- x<y -> argument = x/y --- */
      ret = divee(x,y,&arg);
	 break;

   }

   /* --- Arcus Tangens --- */
   ret = _s_atan(&arg, res);
   }

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   return ret;
} /* atanee2() */





