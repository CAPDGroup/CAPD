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

/* CVS $Id: t_xtre.c,v 1.22 2014/01/30 17:24:17 cxsc Exp $ */

/*****************************************************************
**                                                              **
**  extract.c   17.01.91 Baumhof                                **
**                                                              **
**  liefert Mantisse und Exponent als ExtReal                   **
**                                                              **
**  Format:                                                     **
**      int xtracte(const ExtReal *arg, ExtReal *rmant,         **
**                                      ExtReal *rexp);         **
**                                                              **
*****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C
/* #define __CB_SUN4_OS4_HARDWARE */
#endif
#endif

/* #include <stdlib.h> */
/* already included by o_defs.h
#ifdef LINT_ARGS
extern int abs(int);
#else
extern int abs();
#endif
*/

#ifdef ANSI_C
#ifdef LINT_ARGS
int xtracte(const ExtReal *arg, ExtReal *rmant, ExtReal *rexp)
#else
int xtracte(arg, rmant, rexp)
const ExtReal *arg;
      ExtReal *rmant;
      ExtReal *rexp;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int xtracte(ExtReal *arg, ExtReal *rmant, ExtReal *rexp)
#else
int xtracte(arg, rmant, rexp)
ExtReal *arg;
ExtReal *rmant;
ExtReal *rexp;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    unsigned char vz;           /* Vorzeichen von arg bzw. exp */
    short int exp;              /* wahrer Exponent von arg */
    unsigned short int e;       /* Betrag von exp, um 15-s bit nach links */
    int s;                      /* Exponent von exp */

    /* Null ?  ==>  Fehler */
    if (0==cmpee(&Zero, arg))
    {
        copyee(&MInfty, rexp);
        copyee(&Zero, rmant);
        return(ExcDBZ);
    }

    /* Dekomposition */
    vz = (arg->s.exp&ExtSignMask) ? 1 : 0;
    exp = (arg->s.exp&ExtExpMask) - ExtExpBias;

    /* Mantisse */
    copyee(arg, rmant);
    rmant->s.exp = ExtExpBias + (vz ? ExtSignMask : 0);

    /* Exponent, konvertieren nach ExtReal */
    copyee(&Zero, rexp);
    if (exp)
    {
        vz = (exp<0 ? 1 : 0);
        e = abs(exp);
        s = 15;
        while (!(e&0x8000)) { s--; e <<= 1; }
#ifndef __CB_SUN4_OS4_HARDWARE
        rexp->s.exp = ((ExtExp)s + ExtExpBias) + (vz ? ExtSignMask : 0);
        rexp->s.DIGIT(6) = e & 0x0ff;
        rexp->s.DIGIT(7) = (e & 0x0ff00)>>8;
#else
	*(a_btyp *)rexp =
		(((a_btyp)
		  (((ExtExp)s + ExtExpBias) + (vz ? ExtSignMask : 0)))<<16)
		+ (a_btyp)((e<<1)&0x0ffff);
#endif

    }
    return(0);
}





