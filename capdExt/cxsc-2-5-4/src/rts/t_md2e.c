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

/* CVS $Id: t_md2e.c,v 1.23 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_md2e.c                              */
/*                                                              */
/****************************************************************/
/*                                                              */
/*      Routinen:                                               */
/*      int t_md2e(ExtReal *a)  [mod2e]                         */
/*              [-1/2,1/2]   -->  0                             */
/*              +-(1/2,3/2)  -->  1                             */
/*              +-[3/2,5/2]  -->  0 usw.                        */
/*      int t_md4e(ExtReal *a)  [mod4e]                         */
/*              [-1/2,1/2]   -->  0                             */
/*              (1/2,3/2)    -->  1,   -(1/2,3/2)   -->  3      */
/*              +-[3/2,5/2]  -->  2                             */
/*              (5/2,7/2)    -->  3,   -(5/2,7/2)   -->  1      */
/*              +-[7/2,9/2]  -->  0 usw.                        */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

#if SHORTABTYP
#define UONE (a_btyp)1
#else
#define UONE (unsigned)1L
#endif

#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C
/* #define __CB_SUN4_OS4_HARDWARE */
#endif
#endif

/*--------------------------------------------------------------*
 | Modulo                                                       |
 *--------------------------------------------------------------*/

#ifndef __CB_SUN4_OS4_HARDWARE
#ifdef ANSI_C
#ifdef LINT_ARGS
int mod2e(const ExtReal *arg)
#else
int mod2e(arg)
const ExtReal *arg;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int mod2e(ExtReal *arg)
#else
int mod2e(arg)
ExtReal *arg;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    unsigned short int exp;     /* Exponent von arg */
    a_btyp  m[2];     /* Mantisse von arg */
    int bits;                   /* bits d(0) d(-1) von arg */
    short int trueexp;          /* wahrer Exponent von arg */
    a_btyp  rest;     /* Rest d(-2)... von arg */

    /* Dekomposition */
    exp = arg->s.exp&ExtExpMask;
    m[0] = ((((((a_btyp)arg->s.DIGIT(3)<<8)
             + (a_btyp)arg->s.DIGIT(2))<<8)
             + (a_btyp)arg->s.DIGIT(1))<<8)
             + (a_btyp)arg->s.DIGIT(0);
    m[1] = ((((((a_btyp)arg->s.DIGIT(7)<<8)
             + (a_btyp)arg->s.DIGIT(6))<<8)
             + (a_btyp)arg->s.DIGIT(5))<<8)
             + (a_btyp)arg->s.DIGIT(4);

    trueexp = exp-16383;
    if (trueexp >= 62) 
        rest = 0;
    else if (trueexp >= 31)
        rest = m[0] & ((UONE<<(62-trueexp))-1);
    else if (trueexp == 30)
        rest = m[0];
    else if (trueexp >= -1)
        rest = m[0] | (m[1] & ((UONE<<(30-trueexp))-1));
    else
        rest = m[0] | m[1];

    bits = 0;
    if ((trueexp>=0 && trueexp<=31 && (m[1]&(UONE<<(31-trueexp))))
        || (trueexp>=32 && trueexp<=63 && (m[0]&(UONE<<(63-trueexp)))))
        bits |= 2;
    if ((trueexp>=-1 && trueexp<=30 && (m[1]&(UONE<<(30-trueexp))))
        || (trueexp>=31 && trueexp<=62 && (m[0]&(UONE<<(62-trueexp)))))
        bits |= 1;

    switch (bits)
    {
        case 0:
        case 3: return(0);
        case 1: return(rest ? 1 : 0);
        case 2: return(1);
    }
    /*NOTREACHED*/
    return 0;
}


#ifdef ANSI_C
#ifdef LINT_ARGS
int mod4e(const ExtReal *arg)
#else
int mod4e(arg)
const ExtReal *arg;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int mod4e(ExtReal *arg)
#else
int mod4e(arg)
ExtReal *arg;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    unsigned char vz;           /* Vorzeichen von arg */
    unsigned short int exp;     /* Exponent von arg */
    a_btyp  m[2];     /* Mantisse von arg */
    int bits;                   /* bits d(1) d(0) d(-1) von arg */
    short int trueexp;          /* wahrer Exponent von arg */
    a_btyp  rest;     /* Rest d(-2)... von arg */
    int res = 0;                /* Rueckgabe */

    /* Dekomposition */
    vz = arg->s.exp&ExtSignMask ? 1 : 0;
    exp = arg->s.exp&ExtExpMask;
    m[0] = ((((((a_btyp)arg->s.DIGIT(3)<<8)
             + (a_btyp)arg->s.DIGIT(2))<<8)
             + (a_btyp)arg->s.DIGIT(1))<<8)
             + (a_btyp)arg->s.DIGIT(0);
    m[1] = ((((((a_btyp)arg->s.DIGIT(7)<<8)
             + (a_btyp)arg->s.DIGIT(6))<<8)
             + (a_btyp)arg->s.DIGIT(5))<<8)
             + (a_btyp)arg->s.DIGIT(4);

    trueexp = exp-16383;
    if (trueexp >= 62) 
        rest = 0;
    else if (trueexp >= 31)
        rest = m[0] & ((UONE<<(62-trueexp))-1);
    else if (trueexp == 30)
        rest = m[0];
    else if (trueexp >= -1)
        rest = m[0] | (m[1] & ((UONE<<(30-trueexp))-1));
    else
        rest = m[0] | m[1];

    bits = 0;
    if ((trueexp>=1 && trueexp<=32 && (m[1]&(UONE<<(32-trueexp))))
        || (trueexp>=33 && trueexp<=64 && (m[0]&(UONE<<(64-trueexp)))))
        bits |= 4;
    if ((trueexp>=0 && trueexp<=31 && (m[1]&(UONE<<(31-trueexp))))
        || (trueexp>=32 && trueexp<=63 && (m[0]&(UONE<<(63-trueexp)))))
        bits |= 2;
    if ((trueexp>=-1 && trueexp<=30 && (m[1]&(UONE<<(30-trueexp))))
        || (trueexp>=31 && trueexp<=62 && (m[0]&(UONE<<(62-trueexp)))))
        bits |= 1;

    switch (bits)
    {
        case 0:
        case 7: res = 0; break;
        case 1: res = rest ? 1 : 0; break;
        case 2: res = 1; break;
        case 3:
        case 4: res = 2; break;
        case 5: res = rest ? 3 : 2; break;
        case 6: res = 3; break;
    }
    if (res && vz)
        return(4-res);
    else
        return(res);
}


#else   /* __CB_SUN4_OS4_HARDWARE */


#ifdef ANSI_C
#ifdef LINT_ARGS
int mod2e(const ExtReal *arg)
#else
int mod2e(arg)
const ExtReal *arg;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int mod2e(ExtReal *arg)
#else
int mod2e(arg)
ExtReal *arg;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    unsigned short int exp;     /* Exponent von arg */
    a_btyp  m[4];     /* Mantisse von arg, 17+32+32+32 bit */
    int bits;                   /* bits d(0) d(-1) von arg */
    short int trueexp;          /* wahrer Exponent von arg */
    a_btyp  rest;     /* Rest d(-2)... von arg */

    /* Dekomposition */
    exp = arg->s.exp&ExtExpMask;
    m[0] = ((*(a_btyp *)arg)&0x0ffff) + 0x010000;        /* 17 bit */
    m[1] = *(((a_btyp *)arg)+1);                         /* 3*32 bit */
    m[2] = *(((a_btyp *)arg)+2);
    m[3] = *(((a_btyp *)arg)+3);

    trueexp = exp-16383;
    if (trueexp >= 111) 
        rest = 0;
    else if (trueexp >= 80)
        rest = m[3] & ((UONE<<(111-trueexp))-1);
    else if (trueexp == 79)
        rest = m[3];
    else if (trueexp >= 48)
        rest = m[3] | (m[2] & ((UONE<<(79-trueexp))-1));
    else if (trueexp == 47)
        rest = m[3] | m[2];
    else if (trueexp >= 16)
        rest = m[3] | m[2] | (m[1] & ((UONE<<(47-trueexp))-1));
    else if (trueexp == 15)
        rest = m[3] | m[2] | m[1];
    else if (trueexp >= -1)
        rest = m[3] | m[2] | m[1] | (m[0] & ((UONE<<(15-trueexp))-1));
    else
        rest = m[3] | m[2] | m[1] | m[0];

    bits = 0;
    if (trueexp>=81&&trueexp<=112&&(m[3]&(UONE<<(112-trueexp)))
        || trueexp>=49&&trueexp<=80&&(m[2]&(UONE<<(80-trueexp)))
        || trueexp>=17&&trueexp<=48&&(m[1]&(UONE<<(48-trueexp)))
        || trueexp>=0&&trueexp<=16&&(m[0]&(UONE<<(16-trueexp))))
        bits |= 2;
    if (trueexp>=80&&trueexp<=111&&(m[3]&(UONE<<(111-trueexp)))
        || trueexp>=48&&trueexp<=79&&(m[2]&(UONE<<(79-trueexp)))
        || trueexp>=16&&trueexp<=47&&(m[1]&(UONE<<(47-trueexp)))
        || trueexp>=-1&&trueexp<=15&&(m[0]&(UONE<<(15-trueexp))))
        bits |= 1;

    switch (bits)
    {
        case 0:
        case 3: return(0);
        case 1: return(rest ? 1 : 0);
        case 2: return(1);
    }
    /*NOTREACHED*/
    return 0;
}


#ifdef ANSI_C
#ifdef LINT_ARGS
int mod4e(const ExtReal *arg)
#else
int mod4e(arg)
const ExtReal *arg;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int mod4e(ExtReal *arg)
#else
int mod4e(arg)
ExtReal *arg;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    unsigned char vz;           /* Vorzeichen von arg */
    unsigned short int exp;     /* Exponent von arg */
    a_btyp  m[4];     /* Mantisse von arg, 17+32+32+32 bit */
    int bits;                   /* bits d(1) d(0) d(-1) von arg */
    short int trueexp;          /* wahrer Exponent von arg */
    a_btyp  rest;     /* Rest d(-2)... von arg */
    int res = 0;                /* Rueckgabe */

    /* Dekomposition */
    vz = arg->s.exp&ExtSignMask ? 1 : 0;
    exp = arg->s.exp&ExtExpMask;
    m[0] = ((*(a_btyp *)arg)&0x0ffff) + 0x010000;        /* 17 bit */
    m[1] = *(((a_btyp *)arg)+1);                         /* 3*32 bit */
    m[2] = *(((a_btyp *)arg)+2);
    m[3] = *(((a_btyp *)arg)+3);

    trueexp = exp-16383;
    if (trueexp >= 111) 
        rest = 0;
    else if (trueexp >= 80)
        rest = m[3] & ((UONE<<(111-trueexp))-1);
    else if (trueexp == 79)
        rest = m[3];
    else if (trueexp >= 48)
        rest = m[3] | (m[2] & ((UONE<<(79-trueexp))-1));
    else if (trueexp == 47)
        rest = m[3] | m[2];
    else if (trueexp >= 16)
        rest = m[3] | m[2] | (m[1] & ((UONE<<(47-trueexp))-1));
    else if (trueexp == 15)
        rest = m[3] | m[2] | m[1];
    else if (trueexp >= -1)
        rest = m[3] | m[2] | m[1] | (m[0] & ((UONE<<(15-trueexp))-1));
    else
        rest = m[3] | m[2] | m[1] | m[0];

    bits = 0;
    if (trueexp>=82&&trueexp<=113&&(m[3]&(UONE<<(113-trueexp)))
        || trueexp>=50&&trueexp<=81&&(m[2]&(UONE<<(81-trueexp)))
        || trueexp>=18&&trueexp<=49&&(m[1]&(UONE<<(49-trueexp)))
        || trueexp>=1&&trueexp<=17&&(m[0]&(UONE<<(17-trueexp))))
        bits |= 4;
    if (trueexp>=81&&trueexp<=112&&(m[3]&(UONE<<(112-trueexp)))
        || trueexp>=49&&trueexp<=80&&(m[2]&(UONE<<(80-trueexp)))
        || trueexp>=17&&trueexp<=48&&(m[1]&(UONE<<(48-trueexp)))
        || trueexp>=0&&trueexp<=16&&(m[0]&(UONE<<(16-trueexp))))
        bits |= 2;
    if (trueexp>=80&&trueexp<=111&&(m[3]&(UONE<<(111-trueexp)))
        || trueexp>=48&&trueexp<=79&&(m[2]&(UONE<<(79-trueexp)))
        || trueexp>=16&&trueexp<=47&&(m[1]&(UONE<<(47-trueexp)))
        || trueexp>=-1&&trueexp<=15&&(m[0]&(UONE<<(15-trueexp))))
        bits |= 1;

    switch (bits)
    {
        case 0:
        case 7: res = 0; break;
        case 1: res = rest ? 1 : 0; break;
        case 2: res = 1; break;
        case 3:
        case 4: res = 2; break;
        case 5: res = rest ? 3 : 2; break;
        case 6: res = 3; break;
    }
    if (res && vz)
        return(4-res);
    else
        return(res);
}


#endif  /* __CB_SUN4_OS4_HARDWARE */
#undef UONE





