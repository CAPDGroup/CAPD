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

/* CVS $Id: t_exp.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_exp.c                               */
/*                                                              */
/*      Entries         : a_real t_exp(arg)                     */
/*                        a_real arg;                           */
/*                                                              */
/*      Arguments       : arg  = argument of exponential        */
/*                                                              */
/*      Description     : Exponential function.                 */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* #include "exp.h"        */
/* ------------------------------------------------------------- */
/* Header Exp und Expm1 Konstanten                               */
/* ------------------------------------------------------------- */
/* Konstanten mit INTERN verglichen */
#define ExpKonstanten  {                 /* c[i]=2**(i/8)           */ \
           EXTREAL(0x3F, 0xFF,              \
                   0x80, 0x00, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00), /* +1.000000000000000E+000 */ \
           EXTREAL(0x3F, 0xFF,              \
                   0x8B, 0x95, 0xC1, 0xE3,  \
                   0xEA, 0x8B, 0xD6, 0xE7), /* +1.090507732665258E+000 */ \
           EXTREAL(0x3F, 0xFF,              \
                   0x98, 0x37, 0xF0, 0x51,  \
                   0x8D, 0xB8, 0xA9, 0x6F), /* +1.189207115002721E+000 */ \
           EXTREAL(0x3F, 0xFF,              \
                   0xA5, 0xFE, 0xD6, 0xA9,  \
                   0xB1, 0x51, 0x38, 0xEA), /* +1.296839554651010E+000 */ \
           EXTREAL(0x3F, 0xFF,              \
                   0xB5, 0x04, 0xF3, 0x33,  \
                   0xF9, 0xDE, 0x64, 0x84), /* +1.414213562373095E+000 */ \
           EXTREAL(0x3F, 0xFF,              \
                   0xC5, 0x67, 0x2A, 0x11,  \
                   0x55, 0x06, 0xDA, 0xDD), /* +1.542210825407941E+000 */ \
           EXTREAL(0x3F, 0xFF,              \
                   0xD7, 0x44, 0xFC, 0xCA,  \
                   0xD6, 0x9D, 0x6A, 0xF4), /* +1.681792830507429E+000 */ \
           EXTREAL(0x3F, 0xFF,              \
                   0xEA, 0xC0, 0xC6, 0xE7,  \
                   0xDD, 0x24, 0x39, 0x2F)  /* +1.834008086409342E+000 */ \
} /* --- Ende ExpKonstanten --- */

/* Konstanten mit INTERN verglichen */
#define Expm1Konstanten {                /* c[i]=2**(i/8) - 1       */ \
           EXTREAL(0x00, 0x00,              \
                   0x00, 0x00, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00), /*  0.000000000000000E+000 */ \
           EXTREAL(0x3F, 0xFB,              \
                   0xB9, 0x5C, 0x1E, 0x3E,  \
                   0xA8, 0xBD, 0x6E, 0x70), /* +9.050773266525766E-002 */ \
           EXTREAL(0x3F, 0xFC,              \
                   0xC1, 0xBF, 0x82, 0x8C,  \
                   0x6D, 0xC5, 0x4B, 0x7A), /* +1.892071150027211E-001 */ \
           EXTREAL(0x3F, 0xFD,              \
                   0x97, 0xFB, 0x5A, 0xA6,  \
                   0xC5, 0x44, 0xE3, 0xA8), /* +2.968395546510096E-001 */ \
           EXTREAL(0x3F, 0xFD,              \
                   0xD4, 0x13, 0xCC, 0xCF,  \
                   0xE7, 0x79, 0x92, 0x11), /* +4.142135623730950E-001 */ \
           EXTREAL(0x3F, 0xFE,              \
                   0x8A, 0xCE, 0x54, 0x22,  \
                   0xAA, 0x0D, 0xB5, 0xBA), /* +5.422108254079409E-001 */ \
           EXTREAL(0x3F, 0xFE,              \
                   0xAE, 0x89, 0xF9, 0x95,  \
                   0xAD, 0x3A, 0xD5, 0xE8), /* +6.817928305074291E-001 */ \
           EXTREAL(0x3F, 0xFE,              \
                   0xD5, 0x81, 0x8D, 0xCF,  \
                   0xBA, 0x48, 0x72, 0x5E)  /* +8.340080864093424E-001 */ \
} /* --- Ende Expm1Konstanten --- */

/* StdFctReal(t_exp,expee) */
#ifdef LINT_ARGS
a_real t_exp(a_real arg)
#else
a_real t_exp(arg)

a_real arg;
#endif
        {
        int      rnd, rc;
        a_real   res;
        ExtReal  a, r;

        E_SPUSH("t_exp")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg, &a);

        if ((rc = expee(&a, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
           ieee_abortr1(rc, &arg);

        setrndmode(rnd);

        E_SPOPP("t_exp")
        return res;
        }
/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | expee                                                        |
 | res = 2**j * c[i] * 2**t                                     |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int expee (const ExtReal *arg, ExtReal *res)
#else
int expee (arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int expee (ExtReal *arg, ExtReal *res)
#else
int expee (arg, res)
ExtReal   *arg;
ExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
#if SUN4_CPP_C
   static ExtReal c[8];
   static char cc[8][sizeof(ExtReal)] = ExpKonstanten;
#else
   static const ExtReal c[8]=ExpKonstanten;     /* c[i]=2**(i/8) */
#endif
#else
   static ExtReal c[8];
   static char cc[8][sizeof(ExtReal)] = ExpKonstanten;
#endif
   ExtReal     absa;       /* absolutes Argument                */
   ExtReal     t;          /* reduziertes Argument fuer Potenzreihe 2**t  */
   DReal       td;         /* wie t                             */
   ExtReal     j;          /* ganzzahliger Exponent 2**j        */
   int         i;          /* Index fuer Konstanten             */
   ExtReal     ie;         /* wie ie                            */
   DReal       v;          /* absa*8/ln2                        */
   ExtReal     ve;         /* wie v                             */
   DReal       gza_v;      /* ganzzahliger Anteil von v         */
   ExtReal     gza_ve;     /* wie gza_v                         */
   ExtReal     gza8tel;    /* gza_v / 8                         */
   ExtReal     r;          /* Ergebnis 2**t                     */
   int         sgn;        /* Vorzeichen des Argumentes arg     */
   int         rnd;        /* RundungsMode                      */
   int         ret;        /* Rueckgabe                         */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

#ifndef ANSI_C
   {
      int i;
      for (i=0; i<8; i++)
         memcpy(&(c[i]), cc[i], sizeof(c[0]));
   }
#else
#if SUN4_CPP_C
   {
      int i;
      for (i=0; i<8; i++)
         memcpy(&(c[i]), cc[i], sizeof(c[0]));
   }
#endif
#endif

   /* --- pruefe Argument --- */
   ArgCheck1(Exp, arg, res);

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* --- absa := |a| --- */
   sgn = SGNE(arg);
   absee(arg, &absa);

   /* --- v := arg*8/ln(2) --- */
   mulLdE(&absa, &v);            /* v:=arg*(LdE=log2(e)=1/ln(2)) */
   v.e += 3;                     /* v*=8                         */

   /* --- gza_v = int(v) --- */
   setrndmode(DOWN);
   dreal_to_extreal(&v, &ve);
   rndintee(&ve, &gza_ve);
   extreal_to_dreal(&gza_ve, &gza_v);

   /* --- t := (v - gza_v)/8 --- */
   setrndmode(NEAR);
   subdd(&v, &gza_v, &td);
   dreal_to_extreal(&td, &t);
   scaliee(&t, -3, &t);          /* t/=8                        */

   /* --- r := 2**t = (2**t-1)+1 --- */
   ret = _s_2xm1(&t, &r);
   addee(&r, &One, &r);

   /* --- j := int(gza_v/8) --- */
   scaliee(&gza_ve, -3, &gza8tel); /* gza8tel = gza_ve / 8      */
   setrndmode(DOWN);
   rndintee(&gza8tel, &j);

   /* --- i := index(gza_v, j) = int(gza_v/8-j)*8 --- */
   subee(&gza8tel, &j, &ie);
   scaliee(&ie, 3, &ie);         /* i*=8                        */
   extreal_to_int(&ie, &i);

   /* --- res := 2**j * c[i] * r --- */
   setrndmode(NEAR);
   scalee(&r, &j, &r);
   mulee(&r, &c[i], res);

   /* --- if arg < 0 ==> res := 1/res --- */
   if(NEG==sgn)
      divee(&One, res, res);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   return ret;
} /* expee() */





