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

/* CVS $Id: t_sinv.c,v 1.22 2014/01/30 17:24:17 cxsc Exp $ */

/*****************************************************************
**                                                              **
**  exam.c   16.01.91 Baumhof                                   **
**                                                              **
**  Zahlklassifizierung (fxam im 80387)                         **
**                                                              **
**  Format:                                                     **
**      int _s_xam(const ExtReal *arg);                         **
**      int _s_chk_invalid(const ExtReal *arg);                 **
**                                                              **
*****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

#if SHORTABTYP
#define UMSB (a_btyp)0x80000000
#else
#define UMSB 0x80000000L
#endif

#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C
/* #define __CB_SUN4_OS4_HARDWARE */
#endif
#endif


#ifndef __CB_SUN4_OS4_HARDWARE


#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_xam(const ExtReal *arg)
#else
int _s_xam(arg)
const ExtReal *arg;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_xam(ExtReal *arg)
#else
int _s_xam(arg)
ExtReal *arg;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    unsigned char vz;           /* Vorzeichen von arg */
    unsigned short int exp;     /* Exponent von arg */
    a_btyp mant[2];  /* Mantisse von arg, 64=32+32 bit */

    /* Dekomposition */
    vz = (arg->s.exp&ExtSignMask) ? 1 : 0;
    exp = arg->s.exp&ExtExpMask;
    mant[0] = ((((((a_btyp)arg->s.DIGIT(3)<<8)
                + (a_btyp)arg->s.DIGIT(2))<<8)
                + (a_btyp)arg->s.DIGIT(1))<<8)
                + (a_btyp)arg->s.DIGIT(0);
    mant[1] = ((((((a_btyp)arg->s.DIGIT(7)<<8)
                + (a_btyp)arg->s.DIGIT(6))<<8)
                + (a_btyp)arg->s.DIGIT(5))<<8)
                + (a_btyp)arg->s.DIGIT(4);

    /* Unendlich */
    if (exp==32767 && !mant[0] &&
        mant[1]==(a_btyp)UMSB)
       return(vz ? MInf : PInf);

    /* NAN */
    else if (exp==32767)
        return(vz ? MNAN : PNAN);

    /* Null */
    else if (!(mant[0]|mant[1]))
        return(vz ? MZero : PZero);

    /* denormal oder pseudodenormal */
    else if (exp==0)
        return(vz ? MDenorm : PDenorm);

    /* 8087 unnormal (unsupported) */
    else if (!(mant[1] & UMSB))
        return(vz ? MUnorm : PUnorm);

    /* normal */
    else
        return(vz ? MNorm : PNorm);
}



#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_chk_invalid(const ExtReal *arg)
#else
int _s_chk_invalid(arg)
const ExtReal *arg;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_chk_invalid(ExtReal *arg)
#else
int _s_chk_invalid(arg)
ExtReal *arg;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    return(_s_xam(arg)&(PUnorm|MUnorm|PNAN|MNAN|PInf|MInf));
}


#else	/* __CB_SUN4_OS4_HARDWARE */


#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_xam(const ExtReal *arg)
#else
int _s_xam(arg)
const ExtReal *arg;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_xam(ExtReal *arg)
#else
int _s_xam(arg)
ExtReal *arg;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    unsigned char vz;           /* Vorzeichen von arg */
    unsigned short int exp;     /* Exponent von arg */
    a_btyp mant[4];  /* Mantisse von arg, 17+32+32+32 bit */

    /* Dekomposition */
    vz = (arg->s.exp&ExtSignMask) ? 1 : 0;
    exp = arg->s.exp&ExtExpMask;
    mant[0] = ((*(a_btyp *)arg)&0x0ffff) /*+0x010000*/ ;	/* 17 bit */
    mant[1] = *(((a_btyp *)arg)+1);			/* 3*32 bit */
    mant[2] = *(((a_btyp *)arg)+2);
    mant[3] = *(((a_btyp *)arg)+3);

    /* Unendlich */
    if (exp==32767 && !(mant[0]|mant[1]|mant[2]|mant[3]))
       return(vz ? MInf : PInf);

    /* NAN */
    else if (exp==32767)
        return(vz ? MNAN : PNAN);

    /* Null */
    else if (!exp && !(mant[0]|mant[1]|mant[2]|mant[3]))
        return(vz ? MZero : PZero);

    /* denormal oder pseudodenormal */
    else if (exp==0)
        return(vz ? MDenorm : PDenorm);

    /* normal */
    else
        return(vz ? MNorm : PNorm);
}



#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_chk_invalid(const ExtReal *arg)
#else
int _s_chk_invalid(arg)
const ExtReal *arg;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_chk_invalid(ExtReal *arg)
#else
int _s_chk_invalid(arg)
ExtReal *arg;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    return(_s_xam(arg)&(PUnorm|MUnorm|PNAN|MNAN|PInf|MInf));
}


#endif	/* __CB_SUN4_OS4_HARDWARE */
#undef UMSB





