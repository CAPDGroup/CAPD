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

/* CVS $Id: t_scee.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_scee.c                              */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* #include <stdlib.h> */
/* already included by o_defs.h
#ifdef LINT_ARGS
extern int abs(int);
#else
extern int abs();
#endif
*/

/*--------------------------------------------------------------*
 | scaliee                                                      |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int scaliee (const ExtReal *arg, int intarg, ExtReal *res)
#else
int scaliee (arg, intarg, res)
const ExtReal  *arg;
      int      intarg;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int scaliee (ExtReal *arg, int intarg, ExtReal *res)
#else
int scaliee (arg, intarg, res)
ExtReal  *arg;
int      intarg;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int         e;          /* Exponent des Argumentes           */

   /* --- bei 0 fertig --- */
   if(0==cmpee(arg, &Zero))
      return copyee(arg, res);

   /* --- extrahiere Exponenten --- */
   e = (int)((arg->s.exp & (ExtExpMask)) - ExtExpBias);

   /* --- Pruefe auf DOMAIN --- */
   if(abs(e+intarg) > ExtExpBias) return DOMAIN;

   /* --- Setzt res --- */
   copyee(arg, res);
   res->s.exp += intarg;

   /* --- kein Fehler mehr moeglich --- */
   return NoErr;
} /* scaliee() */

/*--------------------------------------------------------------*
 | scalee                                                       |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int scalee (const ExtReal *arg, const ExtReal *intarg, ExtReal *res)
#else
int scalee (arg, intarg, res)
const ExtReal  *arg;
const ExtReal  *intarg;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int scalee (ExtReal *arg, ExtReal *intarg, ExtReal *res)
#else
int scalee (arg, intarg, res)
ExtReal  *arg;
ExtReal  *intarg;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int         a;          /* (int) intarg                      */
   int         rnd;        /* RundungsMode                      */
   int         ret;        /* Rueckgabe                         */

   /* --- RundungsMode sichern, CHOP setzen --- */
   rnd  = getrndmode();
   setrndmode(CHOP);

   /* --- Konvertierung nach int, falls fehlerfrei scaliee --- */
   if(NoErr==(ret=extreal_to_int(intarg, &a)))
      ret = scaliee(arg, a, res);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   return ret;
} /* scalee() */





