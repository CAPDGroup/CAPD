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

/* CVS $Id: t_ieee.h,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_ieee.h                              */
/*                                                              */
/*                   o_defs.h to be included                    */
/*                   do not #undef extern                       */
/****************************************************************/
/*#include "c:\cmp\c\include\stdio.h"  */
/*#include "c:\cmp\c\include\math.h"   */
/*#include "c:\cmp\c\include\limits.h" */
/* #include <math.h> */ /* 121090 Baumhof, siehe auch exc.h */
#include <stdio.h>      /* FILE in fprinte()                */
/* entfernt 080591 wegen Kollision von size_t in sys/types.h und o_defs.h */
/* #include <limits.h>  */ /* CHAR_BIT fuer BitsPerDigit       */
#ifndef CHAR_BIT
#define CHAR_BIT 8         /* aus limits.h 080591 cb           */
#endif
#ifdef AIX
/* #ifndef CLOCKS_PER_SEC */
#include "/u/p88c/runtime/o_defs.h"
/* #endif*/
#include "/u/p88c/runtime/tbyte/t_exc_.h"
#include "/u/p88c/runtime/tbyte/t_mach.h"
#else
/* #ifndef CLOCKS_PER_SEC */
#include "o_defs.h"
/* #endif */
#include "t_exc_.h"        /* ieee exception codes             */
#include "t_mach.h"       /* Maschinenabhaengigkeiten         */
#endif
/* ---------------------------------------------------------*/

#ifdef ANSI_C
#define const
#else
#define const
#endif

/* ------------------------------------------------------------- */
/* Typ Declarationen                                             */
/* ------------------------------------------------------------- */
#define LongReal EByte  /* IEEE basic format                     */
#if SUN4_CPP_C
#define ExtReal  TByte  /* IEEE extended format                  */
/*typedef char ExtReal[t_size];*/
#else
#define ExtReal  TByte  /* IEEE extended format                  */
#endif

/* --- Digit --- */
/* moved to o_type.h
typedef unsigned char Digit;
*/
                                             /* eine MantissenStelle  */

/* Anzahl Bits pro Digit */
#define BitsPerDigit         ((int)(sizeof(Digit)*CHAR_BIT))

/* --- ExtReal Format --- */   /* 030990 Baumhof "short" eingefuegt */
/* moved to o_type.h
typedef unsigned short int  ExtExp;
*/
/* Typ Exponent        */
#define ExtExpBias (ExtExp) (0x3FFF)         /* Exponent Bias       */
#define ExtExpMask (ExtExp) (0x7FFF)         /* Exponent Maske      */
/* moved to o_type.h
#define ExtMantLen 8
*/
                                             /* Laenge Mantisse     */
#define ExtSignMask (ExtExp) (0x8000)        /* Vorzeichen Maske    */

/* --- EightByte: 64-bit-IEEE-Gleitpunktzahl                        */

#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C
/* #define __CB_SUN4_OS4_HARDWARE */
#endif
#endif

#ifdef __CB_SUN4_OS4_HARDWARE
typedef double EByte;
#else
typedef struct { char c[8]; } EByte;
#endif

/* --- TenByte: Typ wird im Intel 80x86 Assembler zur Darstellung von
   Extended Gleitpunktzahlen fuer die 80x87 Coprozessoren verwendet ---  */

/* moved to o_type.h
typedef union {
   char c[10];
   struct s{
#if LSBFIRST
      Digit  digit[ExtMantLen];
      ExtExp exp;
#else
      ExtExp exp;
      Digit  digit[ExtMantLen];
#endif
   } s;
} TByte; 
*/

/* Makro zum Zugriff auf die Digits im Tenbyte-Format  080391 cb */
/* Verwendung: tenbyte.s.DIGIT(0..7) */
#if LSBFIRST
#  define DIGIT(a) digit[a]
#else
#  define DIGIT(a) digit[7-(a)]
#endif

/* --- Intervalle --- */
typedef struct IExtReal{
   ExtReal u;
   ExtReal l;
} IExtReal;

typedef struct ILongReal{
   LongReal u;
   LongReal l;
} ILongReal;
/* --- end of type declarations --- */

/* ------------------------------------------------------------- */
/* Grundoperationen umleiten auf Cordes-80-bit-Routinen          */
/* ------------------------------------------------------------- */

#define addee b_tadd
#define subee b_tsub
#define mulee b_tmul
#define divee b_tdiv

#define setrndmode t_srnd
#define getrndmode t_grnd

/* -------------------------------------------------------------- */
/* Spezialfallabfrage bei Intervallfunktionen                     */
/* -------------------------------------------------------------- */
#define INT_HPREC       1

/* -------------------------------------------------------------- */
/* Funktions Declarationen                                        */
/* -------------------------------------------------------------- */

#ifdef LINT_ARGS

/* --- Arithmetik --- */
/* already declared in b_fcth.h
extern int  addee(const  ExtReal *arg1, const  ExtReal *arg2,  ExtReal *res);
extern int  subee(const  ExtReal *arg1, const  ExtReal *arg2,  ExtReal *res);
extern int  mulee(const  ExtReal *arg1, const  ExtReal *arg2,  ExtReal *res);
extern int  divee(const  ExtReal *arg1, const  ExtReal *arg2,  ExtReal *res);
*/

extern int  sqrtee  (const  ExtReal *arg,  ExtReal *res);
extern int  sqrtm1ee(const  ExtReal *arg,  ExtReal *res);
extern int isqrtee  (const IExtReal *arg, IExtReal *res);

/* --- Punktstandardfunktionen-- */
extern int sinee(const ExtReal *arg,  ExtReal *res);
extern int cosee(const ExtReal *arg,  ExtReal *res);
extern int tanee(const ExtReal *arg,  ExtReal *res);
extern int cotee(const ExtReal *arg,  ExtReal *res);

extern int asinee(const ExtReal *arg,  ExtReal *res);
extern int acosee(const ExtReal *arg,  ExtReal *res);
extern int atanee(const ExtReal *arg,  ExtReal *res);
extern int atanee2(const ExtReal *arg,  const ExtReal *y, ExtReal *res);
extern int acotee(const ExtReal *arg,  ExtReal *res);

extern int expm1ee(const ExtReal *arg,  ExtReal *res);
extern int lnp1ee (const ExtReal *arg,  ExtReal *res);

extern int expee (const ExtReal *arg,  ExtReal *res);
extern int lnee  (const ExtReal *arg,  ExtReal *res);
extern int sinhee(const ExtReal *arg,  ExtReal *res);
extern int coshee(const ExtReal *arg,  ExtReal *res);
extern int tanhee(const ExtReal *arg,  ExtReal *res);
extern int cothee(const ExtReal *arg,  ExtReal *res);

extern int asinhee(const ExtReal *arg,  ExtReal *res);
extern int acoshee(const ExtReal *arg,  ExtReal *res);
extern int atanhee(const ExtReal *arg,  ExtReal *res);
extern int acothee(const ExtReal *arg,  ExtReal *res);

extern int powee  (const ExtReal *bas,  const ExtReal *exp, ExtReal *res);
extern int logaee (const ExtReal *arg,  const ExtReal *bas, ExtReal *res);

/* --- Intervallstandardfunktionen-- */
extern int isinee   (const IExtReal *arg,  IExtReal *res);
extern int icosee   (const IExtReal *arg,  IExtReal *res);
extern int itanee   (const IExtReal *arg,  IExtReal *res);
extern int icotee   (const IExtReal *arg,  IExtReal *res);

extern int iasinee  (const IExtReal *arg,  IExtReal *res);
extern int iacosee  (const IExtReal *arg,  IExtReal *res);
extern int iatanee  (const IExtReal *arg,  IExtReal *res);
extern int iacotee  (const IExtReal *arg,  IExtReal *res);

extern int ilnee    (const IExtReal *arg,  IExtReal *res);
extern int iexpee   (const IExtReal *arg,  IExtReal *res);
extern int isinhee  (const IExtReal *arg,  IExtReal *res);
extern int icoshee  (const IExtReal *arg,  IExtReal *res);
extern int itanhee  (const IExtReal *arg,  IExtReal *res);
extern int icothee  (const IExtReal *arg,  IExtReal *res);

extern int iasinhee (const IExtReal *arg,  IExtReal *res);
extern int iacoshee (const IExtReal *arg,  IExtReal *res);
extern int iatanhee (const IExtReal *arg,  IExtReal *res);
extern int iacothee (const IExtReal *arg,  IExtReal *res);

extern int ipowee   (const IExtReal *bas,const IExtReal *exp,IExtReal *res);

#else

/* already declared in b_fcth.h
extern int addee(), divee(), subee(), mulee();
*/
extern int sqrtee(), sqrtm1ee(), isqrtee();

extern int sinee(), cosee(), tanee(), cotee(), asinee(), acosee();
extern int atanee(), atanee2(), acotee(), expm1ee(), lnp1ee(), expee(), lnee();
extern int sinhee(), coshee(), tanhee(), cothee(), asinhee(), acoshee();
extern int atanhee(), acothee(), powee(), logaee();

extern int isinee(), icosee(), itanee(), icotee(), iasinee(), iacosee();
extern int iatanee(), iacotee(), ilnee(), iexpee(), isinhee(), icoshee();
extern int itanhee(), icothee(), iasinhee();
extern int iacoshee(), iatanhee(), iacothee();
extern int ipowee();

#endif /* LINT_ARGS */

/* ------------------------------------------------------------ */
/* weitere Funktionen                                           */
/* ------------------------------------------------------------ */
/* --- Vorzeichen --- */
#define POS 1
#define NEG (-POS)
/* --- ExtReal Vorzeichen, result = NEG or POS --- */
#define SGNE(arg) ((int)(((arg)->s.exp&ExtSignMask)==0? POS: NEG))

#ifdef LINT_ARGS

/*--- convert --- */
extern int  extreal_to_longreal(const  ExtReal  *arg,  LongReal *res);
/* --- unused ---
extern int iextreal_to_longreal(const IExtReal  *arg, ILongReal *res);
*/
extern int  longreal_to_extreal(const  LongReal *arg,  ExtReal  *res);
/* --- unused ---
extern int ilongreal_to_extreal(const ILongReal *arg, IExtReal  *res);
*/

extern int  extreal_to_int(const  ExtReal  *arg,  int *res);

/*--- Vergleich --- */
/* cmp == -1 falls arg1<arg2; 0 falls arg1==arg2; 1 falls arg1>arg2 */
extern int  cmpee    (const  ExtReal  *arg1, const  ExtReal  *arg2);
extern int  cmpabsee (const  ExtReal  *arg1, const  ExtReal  *arg2);

/*--- print ExtReal --- */
extern int  printe(const ExtReal *arg);
extern int  fprinte(FILE *stream, const ExtReal *arg);

/*--- Vorzeichenwechsel --- */
extern int  chsee(const  ExtReal *arg, ExtReal *res);
extern int ichsee(const IExtReal *arg,IExtReal *res);

/*--- Extract: Teile in Mantisse und Exponent --- */
extern int xtracte(const ExtReal *arg, ExtReal *rmant, ExtReal *rexp);
extern int xtrexpe(const ExtReal *arg, int *exponent);

/*--- absolute value --- */
extern int  absee(const  ExtReal *arg, ExtReal *res);
extern int iabsee(const IExtReal *arg,IExtReal *res);

/*--- copy --- */
extern int  copyee(const  ExtReal *arg,  ExtReal *res);
extern int icopyee(const IExtReal *arg, IExtReal *res);

#else

extern int extreal_to_longreal() /* , iextreal_to_longreal()  */ ;
extern int longreal_to_extreal() /* , ilongreal_to_extreal()  */ ;
extern int extreal_to_int();
extern int cmpee(), cmpabsee();
extern int printe(), fprinte();
extern int chsee(), ichsee();
extern int xtracte(), xtrexpe();
extern int absee(), iabsee();
extern int copyee(), icopyee();

#endif /* LINT_ARGS */

/* --------------------------------------------------------------------*/
/* RundungsModi setzen und abfragen                                    */
/*    z.B.: setrndmode(DOWN);             Rundung nach oben setzen     */
/*          if(NEAR==getrndmode()) ...    falls Rundung nach unten ... */
/* --------------------------------------------------------------------*/
#ifdef LINT_ARGS
#ifdef ANSI_C
/* extern int setrndmode(int rndmod);     set rounding control   */
extern int getrndmode(void);           /* get rounding control   */
#else
/* extern int setrndmode(int rndmod); */
extern int getrndmode(void);
#endif
#else  /* NOT LINT_ARGS */
extern int /* setrndmode(), */ getrndmode();
#endif /* LINT_ARGS */

/* --- RoundModes --- */
/* maybe this macros have been already defined in a_defs(Ibm_370 sensibility */
#ifndef DOWN
#define DOWN (-UP)
#endif
#ifndef NEAR
#define NEAR 0
#endif
#ifndef UP  
#define UP   1
#endif
#ifndef CHOP
#define CHOP 2
#endif

#ifdef LINT_ARGS

/* --- RundungsFunktionen --- */
extern int iround_rel(const IExtReal *arg,
                      const ExtReal *eps_rel,IExtReal *res);
extern int iround_abs(const IExtReal *arg,
                      const ExtReal *eps_abs,IExtReal *res);

extern int round_rel(int rnd,const ExtReal *arg,
                     const ExtReal *eps_rel,ExtReal *res);
extern int round_abs(int rnd,const ExtReal *arg,
                     const ExtReal *eps_abs,ExtReal *res);

extern int sround(int rnd, const ExtReal *arg, const ExtReal *eps_rel,
                   const ExtReal *eps_abs, ExtReal *res);

/* --- Rundung auf ganze Zahl --- */
extern int  rndintee (const  ExtReal *arg,  ExtReal *res);
extern int irndintee (const IExtReal *arg, IExtReal *res);

/* --- end of function declarations --- */

#else

extern int iround_rel(), iround_abs(), round_rel(), round_abs(), sround();
extern int rndintee(), irndintee();

#endif /* LINT_ARGS */

/* --- macros ---------------------------------------------------*/
#define ADDULP(arg) addee(arg, &MinExtReal, arg)      /* + ulp   */
#define SUBULP(arg) subee(arg, &MinExtReal, arg)      /* - ulp   */

/* ------------------------------------------------------------------- */
/* FehlerBehandlung                                                    */
/* ieee_matherr liegt als Souce vor und ist vom Anwender veraenderbar  */
/* Rueckgabe                                                           */
/*    TakeOver:                                                        */
/*       Uebernahme type als RueckgabeWert der Abgebrochenen IeeeFunktion */
/*       Uebernahme res als Ergebnis der Abgebrochenen IeeeFunktion    */
/*       keine weitere Standard Fehlerbehandlung                       */
/*    NoTakeOver:                                                      */
/*       Standard Fehlerbehandlung                                     */
/* ------------------------------------------------------------------- */
struct ieee_exception{
         int      type;  /* FehlerArt, kann veraendert werden (s.o.)    */
   const char     *name; /* Name der Abgebrochenen Funktion (string)    */
   const IExtReal *arg1; /* Argument 1                                  */
   const IExtReal *arg2; /* Argument 2                                  */
         IExtReal *res;  /* Ergebnis, kann veraendert werden (TakeOver) */
};
#ifdef LINT_ARGS
extern int ieee_matherr(struct ieee_exception *exception);
#else
extern int ieee_matherr();
#endif /* LINT_ARGS */
#define TakeOver   1
#define NoTakeOver 0

/* -------------------------------------------------------------------- */
/* Konstanten                                                           */
/* -------------------------------------------------------------------- */
extern const ExtReal Zero;        /* 0                                  */
extern const ExtReal Half;        /* 1/2                                */
extern const ExtReal Quarter;     /* 1/4                                */
extern const ExtReal ThreeQuart;  /* 3/4                                */
extern const ExtReal One;         /* 1                                  */
extern const ExtReal Two;         /* 2                                  */
extern const ExtReal Four;        /* 4                                  */
extern const ExtReal Eight;       /* 8                                  */
extern const ExtReal MinusOne;    /* -1                                 */
extern const ExtReal MinExtReal;  /* kleinste denormale ExtReal         */
extern const ExtReal MinNormExtReal;/* kleinste normale ExtReal (positiv) */
extern const ExtReal MaxNormExtReal;/* groesste normale ExtReal (positiv) */
extern const ExtReal SqrtHalf;    /* 1/sqrt(2)                (near)    */
extern const ExtReal TwoPow63;    /* 2**63                              */
extern const ExtReal TwoPow127;   /* 2**127                             */
extern const ExtReal Pi;          /* pi                       (near)    */
extern const ExtReal MinusPi;     /* -pi                      (near)    */
extern const ExtReal PiHalf;      /* pi/2                     (near)    */
extern const ExtReal MinusPiHalf; /* -pi/2                    (near)    */
extern const ExtReal PiQuart;     /* pi/4                     (near)    */
extern const ExtReal PiDiv8;      /* pi/8                     (near)    */
extern const ExtReal PiDiv16;     /* pi/16                    (near)    */
extern const ExtReal Pi3Div16;    /* pi*3/16                  (near)    */
extern const ExtReal SqrtTwo;     /* sqrt(2)                  (near)    */
extern const ExtReal Ln2;         /* ln(2)                    (near)    */
/* log2(e)=1/ln(2)          (near)    */
/* --- eigene Konstante in dexp ---
extern const ExtReal LdE;
*/
extern const ExtReal PInfty;      /* +Infinity                          */
extern const ExtReal MInfty;      /* -Infinity                          */

extern const IExtReal IOne;
extern const IExtReal IHalf;
extern const IExtReal IPi;
extern const IExtReal IPiHalf;
extern const IExtReal ITwoDivPi;
extern const IExtReal ISqrtTwo;

/* ------------------------------------------------------------------ */

#ifdef ANSI_C
#else
#undef const
/*
#undef extern
*/
#endif

/* ------------------------------------------------------------------- */





