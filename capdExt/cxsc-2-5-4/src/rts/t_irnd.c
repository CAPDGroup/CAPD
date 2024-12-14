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

/* CVS $Id: t_irnd.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_irnd.c                              */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/*--------------------------------------------------------------*
 | iround_rel, arg*(1+eps_rel)) bzw arg*(1-eps_rel))            |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int iround_rel(const IExtReal *arg, const ExtReal *eps_rel, IExtReal *res)
#else
int iround_rel(arg, eps_rel, res)
const IExtReal  *arg;
const  ExtReal  *eps_rel;
      IExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int iround_rel(IExtReal *arg, ExtReal *eps_rel, IExtReal *res)
#else
int iround_rel(arg, eps_rel, res)
IExtReal  *arg;
 ExtReal  *eps_rel;
IExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   round_rel(UP,   &(arg->u), eps_rel, &(res->u));
   round_rel(DOWN, &(arg->l), eps_rel, &(res->l));

   return NoErr;
} /* iround_rel() */

/*--------------------------------------------------------------*
 | iround_abs, arg+eps_abs bzw arg+eps_abs                      |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int iround_abs(const IExtReal *arg, const ExtReal *eps_abs, IExtReal *res)
#else
int iround_abs(arg, eps_abs, res)
const IExtReal  *arg;
const  ExtReal  *eps_abs;
      IExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int iround_abs(IExtReal *arg, ExtReal *eps_abs, IExtReal *res)
#else
int iround_abs(arg, eps_abs, res)
IExtReal  *arg;
 ExtReal  *eps_abs;
IExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   round_abs(UP,   &(arg->u), eps_abs, &(res->u));
   round_abs(DOWN, &(arg->l), eps_abs, &(res->l));

   return NoErr;
} /* iround_abs() */

/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int round_rel(int rndmod, const ExtReal *arg, const ExtReal *eps_rel,
              ExtReal *res)
#else
int round_rel(rndmod, arg, eps_rel, res)
      int      rndmod;
const ExtReal  *arg;
const ExtReal  *eps_rel;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int round_rel(int rndmod, ExtReal *arg, ExtReal *eps_rel, ExtReal *res)
#else
int round_rel(rndmod, arg, eps_rel, res)
int      rndmod;
ExtReal  *arg;
ExtReal  *eps_rel;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal  e;
   int      rnd;

   rnd  = getrndmode();     /* Rundungs-Mode sichern */

   setrndmode(NEAR);
   mulee(arg, eps_rel, &e);
   absee(&e, &e);
   setrndmode(rndmod);
   if (rndmod==UP)
      addee(arg, &e, res);
   else if (rndmod==DOWN)
      subee(arg, &e, res);

   setrndmode(rnd); /* auf urspruenglichen Mode zuruecksetzten */

   return NoErr;
} /* round_rel() */

/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int round_abs(int rndmod, const ExtReal *arg, const ExtReal *eps_abs,
              ExtReal *res)
#else
int round_abs(rndmod, arg, eps_abs, res)
      int      rndmod;
const ExtReal  *arg;
const ExtReal  *eps_abs;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int round_abs(int rndmod, ExtReal *arg, ExtReal *eps_abs, ExtReal *res)
#else
int round_abs(rndmod, arg, eps_abs, res)
int      rndmod;
ExtReal  *arg;
ExtReal  *eps_abs;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int      rnd;

   rnd  = getrndmode();     /* Rundungs-Mode sichern */

   setrndmode(rndmod);
   if(rndmod==UP)
      addee(arg, eps_abs, res);
   else if(rndmod==DOWN)
      subee(arg, eps_abs, res);
   setrndmode(rnd); /* auf urspruenglichen Mode zuruecksetzten */

   return NoErr;
} /* round_abs() */

/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int sround(int rndmod, const ExtReal *arg, const ExtReal *eps_rel,
           const ExtReal *eps_abs, ExtReal *res)
#else
int sround(rndmod, arg, eps_rel, eps_abs, res)
      int      rndmod;
const ExtReal  *arg;
const ExtReal  *eps_rel;
const ExtReal  *eps_abs;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int sround(int rndmod, ExtReal *arg, ExtReal *eps_rel, ExtReal *eps_abs,
           ExtReal *res)
#else
int sround(rndmod, arg, eps_rel, eps_abs, res)
int      rndmod;
ExtReal  *arg;
ExtReal  *eps_rel;
ExtReal  *eps_abs;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal  e;
   int      rnd;

   rnd  = getrndmode();     /* Rundungs-Mode sichern */

   setrndmode(NEAR);
   mulee(arg, eps_rel, &e);
   absee(&e, &e);
   setrndmode(rndmod);
   if(rndmod==UP) {
      addee(arg, &e, res);
      addee(res, eps_abs, res);
   }
   else if(rndmod==DOWN) {
      subee(arg, &e, res);
      subee(res, eps_abs, res);
   }
   setrndmode(rnd); /* auf urspruenglichen Mode zuruecksetzten */

   return NoErr;
} /* round_abs() */





