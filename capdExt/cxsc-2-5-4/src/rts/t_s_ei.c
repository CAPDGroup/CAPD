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

/* CVS $Id: t_s_ei.c,v 1.22 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_s_ei.c                              */
/*                                                              */
/*                   ambigous assignment removed                */
/****************************************************************/

/*****************************************************************
**                                                              **
**  eetoi.c   15.01.91 Baumhof                                  **
**                                                              **
**  Rundung ExtReal (10 Byte) <--> Integer                      **
**                                                              **
**  Format:                                                     **
**      int _s_etoi(const ExtReal *arg, int *res);              **
**        Rundung: Chop                                         **
**      int _s_etoie(const ExtReal *arg, ExtReal *res);         **
**        Rundung gemaess b_rflg                                **
**                                                              **
*****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* aus t.h Tenbyte-Arithmetik D. Cordes */
extern int b_rflg;              /* Rundung: Down, Near, Up, Chop */

#if SHORTABTYP
#define UONE (a_btyp)1
#define UMSB (a_btyp) 0x80000000
#else
#define UONE (unsigned)1L
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
int _s_etoi(const ExtReal *arg, int *res)
#else
int _s_etoi(arg, res)
const ExtReal *arg;
      int     *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_etoi(ExtReal *arg, int *res)
#else
int _s_etoi(arg, res)
ExtReal *arg;
int     *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    unsigned char vz;           /* Vorzeichen von arg */
    unsigned short int exp;     /* Exponent von arg */
    a_btyp m[2];     /* Mantisse von arg, 64=32+32 bit */
    short int trueexp;          /* Wahrer Exponent von arg */

    /* Dekomposition */
    vz = (arg->s.exp&ExtSignMask) ? 1 : 0;
    exp = arg->s.exp&ExtExpMask;
    m[0] = ((((((a_btyp)arg->s.DIGIT(3)<<8)
             + (a_btyp)arg->s.DIGIT(2))<<8)
             + (a_btyp)arg->s.DIGIT(1))<<8)
             + (a_btyp)arg->s.DIGIT(0);
    m[1] = ((((((a_btyp)arg->s.DIGIT(7)<<8)
             + (a_btyp)arg->s.DIGIT(6))<<8)
             + (a_btyp)arg->s.DIGIT(5))<<8)
             + (a_btyp)arg->s.DIGIT(4);
    *res = 0;

    /* Unendlich oder NaN */
    if (exp==32767)
    {
        return(DOMAIN);
    }
    /* Null */
    else if (!(m[0]|m[1]))
    {
        return(0);
    }
    /* pseudodenormal */
    else if (!exp && m[1] & UMSB)
    {
        return(DOMAIN);
    }
    /* denormalisierte Zahl */
    else if (!exp)
    {
        return(0);
    }
    /* 8087 unnormal (unsupported) */
    else if (!(m[1]  & UMSB))
    {
        return(DOMAIN);
    }
    /* normale Zahl */
    else
    {
        trueexp = exp-16383;
        /* Ueberlauf */
        if (trueexp>=31)
        {
            return(OVER_FLOW);
        }
        /* Unterlauf */
        else if (trueexp<=-1)
        {
            return(0);
        }
        /* Zahl passt */
        *res = (int)(m[1]>>(31-trueexp));
        if (vz)
            *res = -(*res);
        return(0);
    }
}



#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_etoie(const ExtReal *arg, ExtReal *res)
#else
int _s_etoie(arg, res)
const ExtReal *arg;
      ExtReal *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_etoie(ExtReal *arg, ExtReal *res)
#else
int _s_etoie(arg, res)
ExtReal *arg;
ExtReal *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    unsigned char vz;           /* Vorzeichen von arg */
    unsigned short int exp;     /* Exponent von arg */
    a_btyp mant[2];  /* Mantisse von arg, 64=32+32 bit */
    unsigned short int e;       /* Exponent von res */
    a_btyp m[2];     /* Mantisse von res, 64=32+32 bit */
    a_btyp rest;     /* != 0, falls nicht exakt */
    a_btyp nbit;     /* Bit, das Near-Rundung entscheidet */
    short int trueexp;          /* Wahrer Exponent von arg */

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
        return(DOMAIN);
    }
    /* Null */
    else if (!(mant[0]|mant[1]))
    {
        e = 0;
        m[0] = m[1] = 0;
    }
    /* pseudodenormal */
    else if (!exp && mant[1] & UMSB)
    {
        return(DOMAIN);
    }
    /* denormalisierte Zahl */
    else if (!exp)
    {
        e=0;m[0]=m[1]=0;
        return(UNDER_FLOW);
    }
    /* 8087 unnormal (unsupported) */
    else if (!(mant[1] & UMSB))
    {
        return(DOMAIN);
    }
    /* normale Zahl */
    else
    {
        trueexp = exp-16383;
        /* Ueberlauf  -->  Zahl ist schon ganz */
        if (trueexp>=63)
        {
            e = exp;
            m[0] = mant[0];
            m[1] = mant[1];
        }
        /* Unterlauf */
        else if (trueexp<=-1)
        {
            /* Ergebnis = 0, auch bei NEAR und Argument=0.5 */
            if (b_rflg==CHOP || (b_rflg==UP&&vz) || (b_rflg==DOWN && !vz)
                || (b_rflg==NEAR &&
                   (trueexp<=-2 ||
                    (trueexp==-1 &&
                    mant[1]==(0x80000000l) &&
                    !mant[0]))))
                {
                e = 0;
                m[0] = m[1] = 0;
            }
            /* Ergebnis = +-1 */
            else
            {
                e = 16383;
                m[0] = 0;
                m[1] = UONE<<31;
            }
        }
        /* Zahl passt */
        else
        {
            s = 63-trueexp;
            if (s>=33)
            {
                nbit = mant[1]&(UONE<<(s-33));
                rest = mant[0] | (mant[1]&((UONE<<(s-33))-1));
                m[0] = 0;
                m[1] = mant[1] & (-(UONE<<(s-32)));
            }
            else
            {
                nbit = mant[0]&(UONE<<(s-1));
                rest = mant[0]&((UONE<<(s-1))-1);
                if (s==32) m[0] = 0;
                else m[0] = mant[0] & (-(UONE<<s));
                m[1] = mant[1];
            }
            e = exp;
        
            /* Rundung notwendig ? */
            if ((rest|nbit) &&  /* es wurde etwas abgeschnitten */
                ((b_rflg==NEAR && nbit &&
                 (rest ||
                  (s>=32 ? m[1]&(UONE<<(s-32)) : m[0]&(UONE<<s))))
                 || (b_rflg==UP&&!vz) || (b_rflg==DOWN&&vz)))
            {
                if (s>=32)
                    m[1] += (UONE<<(s-32));
                else
                {
                    m[0] += (UONE<<s);
                    if (!m[0])
                        m[1]++;
                }
                /* Ueberlauf ? */
                if (!m[1])
                {
                    m[1] = (UONE<<31);
                    e++;
                }
            }
        }
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


#else	/* __CB_SUN4_OS4_HARDWARE */


#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_etoi(const ExtReal *arg, int *res)
#else
int _s_etoi(arg, res)
const ExtReal *arg;
      int     *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_etoi(ExtReal *arg, int *res)
#else
int _s_etoi(arg, res)
ExtReal *arg;
int     *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    unsigned char vz;           /* Vorzeichen von arg */
    unsigned short int exp;     /* Exponent von arg */
    a_btyp m[4];     /* Mantisse von arg, 17+32+32+32 bit */
    short int trueexp;          /* Wahrer Exponent von arg */

    /* Dekomposition */
    vz = (arg->s.exp&ExtSignMask) ? 1 : 0;
    exp = arg->s.exp&ExtExpMask;
    m[0] = ((*(a_btyp *)arg)&0x0ffff) + 0x010000;	/* 17 bit */
    m[1] = *(((a_btyp *)arg)+1);				/* 3*32 bit */
    m[2] = *(((a_btyp *)arg)+2);
    m[3] = *(((a_btyp *)arg)+3);
    *res = 0;

    /* Unendlich oder NaN */
    if (exp==32767)
    {
        return(DOMAIN);
    }
    /* Null */
    else if (!exp && !((m[0]&0x0ffff)|m[1]|m[2]|m[3]))
    {
        return(0);
    }
    /* denormalisierte Zahl */
    else if (!exp)
    {
        return(0);
    }
    /* normale Zahl */
    else
    {
        trueexp = exp-16383;
        /* Ueberlauf */
        if (trueexp>=31)
        {
            return(OVER_FLOW);
        }
        /* Unterlauf */
        else if (trueexp<=-1)
        {
            return(0);
        }
        /* Zahl passt */
        *res = (int)(((m[0]<<15)+((m[1]>>17)&0x7fff))>>(31-trueexp));
        if (vz)
            *res = -(*res);
        return(0);
    }
}



#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_etoie(const ExtReal *arg, ExtReal *res)
#else
int _s_etoie(arg, res)
const ExtReal *arg;
      ExtReal *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_etoie(ExtReal *arg, ExtReal *res)
#else
int _s_etoie(arg, res)
ExtReal *arg;
ExtReal *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
    unsigned char vz;           /* Vorzeichen von arg */
    unsigned short int exp;     /* Exponent von arg */
    a_btyp mant[4];  /* Mantisse von arg, 17+32+32+32 bit */
    unsigned short int e;       /* Exponent von res */
    a_btyp m[4];     /* Mantisse von erg, 17+32+32+32 bit */
    a_btyp rest;     /* != 0, falls nicht exakt */
    a_btyp nbit;     /* Bit, das Near-Rundung entscheidet */
    short int trueexp;          /* Wahrer Exponent von arg */

    int s;

    /* Dekomposition */
    vz = (arg->s.exp&ExtSignMask) ? 1 : 0;
    exp = arg->s.exp&ExtExpMask;
    mant[0] = ((*(a_btyp *)arg)&0x0ffff) + 0x010000;	/* 17 bit */
    mant[1] = *(((a_btyp *)arg)+1);			/* 3*32 bit */
    mant[2] = *(((a_btyp *)arg)+2);
    mant[3] = *(((a_btyp *)arg)+3);

    /* Unendlich oder NaN */
    if (exp==32767)
    {
        return(DOMAIN);
    }
    /* Null */
    else if (!exp && !((mant[0]&0x0ffff)|mant[1]|mant[2]|mant[3]))
    {
        e = 0;
        m[0] = m[1] = m[2] = m[3] = 0;
    }
    /* denormalisierte Zahl */
    else if (!exp)
    {
        /* e=0; m[0]=m[1]=m[2]=m[3]=0; */
        return(UNDER_FLOW);
    }
    /* normale Zahl */
    else
    {
        trueexp = exp-16383;
        /* Ueberlauf  -->  Zahl ist schon ganz */
        if (trueexp>=112)
        {
            e = exp;
            m[0] = mant[0];
            m[1] = mant[1];
            m[2] = mant[2];
            m[3] = mant[3];
        }
        /* Unterlauf */
        else if (trueexp<=-1)
        {
            /* Ergebnis = 0, auch bei NEAR und Argument=0.5 */
            if (b_rflg==CHOP || (b_rflg==UP&&vz) || (b_rflg==DOWN && !vz)
                || (b_rflg==NEAR &&
                   (trueexp<=-2 ||
                    trueexp==-1 &&
                    mant[0]==0x010000 && !mant[1] && !mant[2] && !mant[3])))
	    {
                e = 0;
                m[0] = m[1] = m[2] = m[3] = 0;
            }
            /* Ergebnis = +-1 */
            else
            {
                e = 16383;
		m[0] = 0x010000;
                m[1] = m[2] = m[3] = 0;
            }
        }
        /* Zahl passt */
        else
        {
            s = 112-trueexp;	/* so viele bits fallen weg */
	    if (s <= 32)
	    {
		nbit = mant[3]&(UONE<<(s-1));
		rest = mant[3]&((UONE<<(s-1))-1);
		if (s==32) m[3] = 0; else m[3] = mant[3]&(-(UONE<<s));
		m[0] = mant[0];
		m[1] = mant[1];
		m[2] = mant[2];
	    }
	    else if (s <= 64)
	    {
		nbit = mant[2]&(UONE<<(s-33));
		rest = mant[3] | (mant[2]&((UONE<<(s-33))-1));
		m[3] = 0;
		if (s==64) m[2] = 0; else m[2] = mant[2]&(-(UONE<<(s-32)));
		m[0] = mant[0];
		m[1] = mant[1];
	    }
	    else if (s <= 96)
	    {
		nbit = mant[1]&(UONE<<(s-65));
		rest = mant[3] | mant[2] | (mant[1]&((UONE<<(s-65))-1));
		m[3] = m[2] = 0;
		if (s==96) m[1] = 0; else m[1] = mant[1]&(-(UONE<<(s-64)));
		m[0] = mant[0];
	    }
	    else
	    {
		nbit = mant[0]&(UONE<<(s-97));
		rest = mant[3]|mant[2]|mant[1] | (mant[0]&((UONE<<(s-97))-1));
		m[3] = m[2] = m[1] = 0;
		m[0] = mant[0]&(-(UONE<<(s-96)));
	    }
            e = exp;
        
            /* Rundung notwendig ? */
            if ((rest|nbit) &&  /* es wurde etwas abgeschnitten */
                ((b_rflg==NEAR && nbit &&
                 (rest ||
		  (s<32 ? m[3]&(UONE<<s) :
		  (s<64 ? m[2]&(UONE<<(s-32)) :
		  (s<96 ? m[1]&(UONE<<(s-64)) : m[0]&(UONE<<(s-96)))))))
                 || (b_rflg==UP&&!vz) || (b_rflg==DOWN&&vz)))
            {
		int i = 3;
		while (s>=32) { i--; s-=32; }
		m[i] += (UONE<<s);
		while (i>0 && !m[i])
		    m[--i]++;
		/* Ueberlauf ? */
		if (i==0 && m[0]==0x020000)
		{
		    e++;
		    m[0] = 0x010000;
		}
            }
        }
    }
    /* Komposition */
    *(a_btyp *)res = ((e + (vz ? ExtSignMask : 0)) << 16)
			    + (m[0] & 0x0ffff);
    *(((a_btyp *)res)+1) = m[1];
    *(((a_btyp *)res)+2) = m[2];
    *(((a_btyp *)res)+3) = m[3];

    return(0);
}


#endif	/* __CB_SUN4_OS4_HARDWARE */
#undef UONE
#undef UMSB





