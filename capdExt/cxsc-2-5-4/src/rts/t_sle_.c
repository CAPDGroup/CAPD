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

/* CVS $Id: t_sle_.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_sle_.c                              */
/*                                                              */
/****************************************************************/

/*****************************************************************
**                                                              **
**  etoltoe.c   11.01.91 Baumhof                                **
**                                                              **
**  Konvertierung LongReal (8 Byte) <--> ExtReal (10 Byte)      **
**                                                              **
**  Format:                                                     **
**      int _s_ltoe(const LongReal *arg, ExtReal *res);         **
**      int _s_etol(const ExtReal *arg, LongReal *res);         **
**                                                              **
*****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

#if SHORTABTYP
#define UONE 		(a_btyp)1
#define UMSB 		(a_btyp)0x80000000
#define U1FFFFF	(a_btyp)0x1fffff
#define U0FFFFF	(a_btyp)0x0fffff
#define U7FFFFF	(a_btyp)0x7fffffff
#else
#define UONE 		(unsigned)1L
#define UMSB 		0x80000000L
#define U1FFFFF	0x1fffffL
#define U0FFFFF	0x0fffffL
#define U7FFFFF	0x7fffffffL
#endif

/* aus t.h Tenbyte-Arithmetik D. Cordes */
extern int b_rflg;              /* Rundung: Down, Near, Up, Chop */

#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_ltoe(const LongReal *arg, ExtReal *res)
#else
int _s_ltoe(arg, res)
const LongReal *arg;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_ltoe(LongReal *arg, ExtReal *res)
#else
int _s_ltoe(arg, res)
LongReal *arg;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    Digit lr[8];                /* arg byteweise */
    unsigned char vz;           /* Vorzeichen von arg */
    unsigned short int exp;     /* Exponent von arg */
    a_btyp mant[2];  /* Mantisse von arg, 52=20+32 bit */
    unsigned short int e;       /* Exponent von res */
    a_btyp m[2];     /* Mantisse von res, 64=32+32 bit */

    Digit *p;
    int i, s;
    a_btyp k;

    /* arg nach lr umspeichern */
    for (i=0, p=(Digit *)arg; i<8; )
#if LSBFIRST
        lr[i++] = *(p++);
#else
        lr[7-i++] = *(p++);
#endif
    /* Dekomposition */
    vz = (lr[7]&0x80) ? 1 : 0;
    exp = ((lr[6]&0xf0)>>4) + ((lr[7]&0x7f)<<4); 
    mant[0] = ((((((a_btyp)lr[3]<<8)
                + (a_btyp)lr[2])<<8)
                + (a_btyp)lr[1])<<8)
                + (a_btyp)lr[0];
    mant[1] = (((((a_btyp)lr[6]&0x0f)<<8)
                + (a_btyp)lr[5])<<8)
                + (a_btyp)lr[4];

    /* Unendlich oder NaN */
    if (exp==2047)
    {
        e = 32767;
        m[1] = (UONE<<31) + (mant[1]<<11) + (mant[0]>>21);
        m[0] = mant[0]<<11;
    }
    /* Null */
    else if (!exp && !(mant[0]|mant[1]))
    {
        e = 0;
        m[0] = m[1] = 0;
    }
    /* normale oder denormalisierte Zahl */
    else
    {
        /* denormalisierte Zahl */
        if (!exp)
        {
            exp = 1;
            /* vorderstes Bit in mant[1],mant[0] suchen */
            s=0; k=UONE<<20; while (k>mant[1]) { k>>=1; s++; }
            if (!k) { k=UONE<<31; while (k>mant[0]) { k>>=1; s++; } }
            /* Mantisse um s+11 Bit schieben */
            exp -= s;
            s += 11;
            if (s>=32) { mant[1] = mant[0]; mant[0] = 0; s-=32; }
            mant[1] = (mant[1]<<s) + (mant[0]>>(32-s));
            mant[0] <<= s;
            /* jetzt steht eine 1 im MSB von mant[1] */
        }
        /* normale Zahl */
        else
        {
            /* Mantisse um 11 Bit schieben, MSB dazu */
            mant[1] = (mant[1]<<11) + (mant[0]>>21);
            mant[0] <<= 11;
            mant[1] |= (UONE<<31);
        }
        /* Exponentenanpassung */
        e = exp + 15360;        /* -1023+16383 */
        m[0] = mant[0];
        m[1] = mant[1];
    }
    /* Komposition */
    res->s.exp = e + (vz ? ExtSignMask : 0);
    res->s.DIGIT(0) = (Digit) (m[0] & 0x0ff); m[0] >>= 8;
    res->s.DIGIT(1) = (Digit) (m[0] & 0x0ff); m[0] >>= 8;
    res->s.DIGIT(2) = (Digit) (m[0] & 0x0ff); m[0] >>= 8;
    res->s.DIGIT(3) = (Digit) (m[0]);
    res->s.DIGIT(4) = (Digit) (m[1] & 0x0ff); m[1] >>= 8;
    res->s.DIGIT(5) = (Digit) (m[1] & 0x0ff); m[1] >>= 8;
    res->s.DIGIT(6) = (Digit) (m[1] & 0x0ff); m[1] >>= 8;
    res->s.DIGIT(7) = (Digit) (m[1]);

    return(0);
}



#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_etol(const ExtReal *arg, LongReal *res)
#else
int _s_etol(arg, res)
const ExtReal *arg;
      LongReal *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_etol(ExtReal *arg, LongReal *res)
#else
int _s_etol(arg, res)
ExtReal *arg;
LongReal *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    unsigned char vz;           /* Vorzeichen von arg */
    unsigned short int exp;     /* Exponent von arg */
    a_btyp mant[2];  /* Mantisse von arg, 64=32+32 bit */
    unsigned short int e;       /* Exponent von res */
    a_btyp m[2];     /* Mantisse von res, 52=20+32 bit */
    a_btyp rest;     /* != 0, falls nicht exakt */
    a_btyp nbit;     /* Bit, das Near-Rundung entscheidet */
    short int trueexp;          /* Wahrer Exponent von arg */

    Digit *p;
    int s;

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

    /* Unendlich oder NaN */
    if (exp==32767)
    {
        e = 2047;
        m[1] = (mant[1]>>11)& U0FFFFF;
        m[0] = ((mant[1]&0x07ff)<<21) + (mant[0]>>11);
    }
    /* Null */
    else if (!(mant[0]|mant[1]))
    {
        e = 0;
        m[0] = m[1] = 0;
    }
    /* pseudodenormal */
    else if (!exp && mant[1]& UMSB)
    {
        return(DOMAIN);
    }
    /* denormalisierte Zahl */
    else if (!exp)
    {
        /* return(UNDER_FLOW); */
        if ((vz&&b_rflg==DOWN)||(!vz&&b_rflg==UP))
        {
            m[0] = UONE;
            m[1] = 0;
            e = 0;
        }
        else
        {
            e = 0;
            m[0] = m[1] = 0;
        }
    }
    /* 8087 unnormal (unsupported) */
    else if (!(mant[1]& UMSB))
    {
        return(DOMAIN);
    }
    /* normale Zahl */
    else
    {
        trueexp = exp-16383;
        /* Ueberlauf */
        if (trueexp>=1024)
        {
            return(OVER_FLOW);
        }
        /* Unterlauf */
        else if (trueexp<=-1075)
        {
            /* return(UNDER_FLOW); */
            if ((vz&&b_rflg==DOWN)||(!vz&&b_rflg==UP)
                ||(trueexp==-1075 && b_rflg==NEAR))
            {
                m[0] = UONE;
                m[1] = 0;
                e = 0;
            }
            else
            {
                e = 0;
                m[0] = m[1] = 0;
            }
            s = 0;
        }
        /* muss denormalisiert werden ? */
        else if (trueexp<=-1023)
        {
            s = -1022-trueexp+11;
            e = 0;
        }
        /* Zahl passt */
        else
        {
            s = 11;
            e = trueexp+1023;
        }
        /* Mantisse um s Bit nach rechts schieben,
           vordere Eins noch dran */
        /* nbit = hoechstwertiges abgeschnittenes Bit */
        /* rest zeigt an, ob danach auch noch Bits weggefallen sind */
        if (s>0)
        {
           if (s>=33)
               nbit = mant[1]&(UONE<<(s-33));
           else
               nbit = mant[0]&(UONE<<(s-1));
           if (s>=32)
           {
               if (s==32) rest = mant[0] & U7FFFFF ;
               else rest = mant[0];
               mant[0] = mant[1];
               mant[1] = 0;
               s -= 32;
           }
           else
               rest = 0;
           if (s>=2) rest |= (mant[0]&((UONE<<(s-1))-1));
           m[0] = (mant[0]>>s) + (mant[1]<<(32-s));
           m[1] = (mant[1]>>s) & U1FFFFF;

           /* Rundung: Addition von 1 notwendig ? */
           if ((rest|nbit) &&      /* es wurde etwas abgeschnitten */
               ((b_rflg==NEAR&&nbit&&(rest||m[0]&1))
                   /* NEAR immer auf gerade runden, falls in der Mitte */
                || (b_rflg==UP&&!vz) || (b_rflg==DOWN&&vz)))
           {
               m[0]++;
               /* Uebertrag nach mant[1] ? */
               if (!m[0])
               {
                   m[1]++;
                   /* Ueberlauf ? */
                   if (!e && m[1]==(UONE<<20))
                   {
                       e = 1;      /* diese Zahl ist wieder normalisiert */
                   }
                   else if (e && m[1]==(UONE<<21))
                   {
                       m[1] >>= 1;
                       e++;
                   }
               }
           }
        }
    }
    /* Komposition, vordere Eins weg, nach res speichern */
#if LSBFIRST
    p = (Digit *)res;
    *(p++) = (Digit) (m[0]&0x0ff); m[0] >>= 8;
    *(p++) = (Digit) (m[0]&0x0ff); m[0] >>= 8;
    *(p++) = (Digit) (m[0]&0x0ff); m[0] >>= 8;
    *(p++) = (Digit) (m[0]);
    *(p++) = (Digit) (m[1]&0x0ff); m[1] >>= 8;
    *(p++) = (Digit) (m[1]&0x0ff); m[1] >>= 8;
    *(p++) = (Digit) ((m[1]&0x00f) + ((e&0x00f)<<4));
    *(p++) = (Digit) (((e>>4)&0x07f) + (vz ? 0x080 : 0));
#else
    p = ((Digit *)res) + 8;
    *(--p) = (Digit) (m[0]&0x0ff); m[0] >>= 8;
    *(--p) = (Digit) (m[0]&0x0ff); m[0] >>= 8;
    *(--p) = (Digit) (m[0]&0x0ff); m[0] >>= 8;
    *(--p) = (Digit) (m[0]);
    *(--p) = (Digit) (m[1]&0x0ff); m[1] >>= 8;
    *(--p) = (Digit) (m[1]&0x0ff); m[1] >>= 8;
    *(--p) = (Digit) ((m[1]&0x00f) + ((e&0x00f)<<4));
    *(--p) = (Digit) (((e>>4)&0x07f) + (vz ? 0x080 : 0));
#endif

    return(0);
}
#undef UONE
#undef UMSB
#undef U1FFFFF	
#undef U0FFFFF
#undef U7FFFFF





