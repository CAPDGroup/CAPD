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

/* CVS $Id: t_defs.h,v 1.22 2014/01/30 17:24:15 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_defs.h                              */
/*                                                              */
/*      Description     : Definitions for tenbyte routines      */
/*                                                              */
/****************************************************************/


/* ---------------------------------------------------------------------- */
/* Header NUR zur Entwicklung der ieee-Funktionen                         */
/* wird vom Anwender der Funktionen NICHT benoetigt !                     */
/* ---                                                                    */
/* dem Anwender steht nur ieee.h und die dort definierten Funktionen zur  */
/* Verfuegung !                                                           */
/* ---------------------------------------------------------------------- */

/* #include <stdlib.h>     fuer Macro max()  (Nur Microsoft!!)            */
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_ieee.h"
#include "/u/p88c/runtime/tbyte/t_drea.h"
#include "/u/p88c/runtime/tbyte/t_cond.h"
#else
#include "t_ieee.h"       /*                                           */
#include "t_drea.h"       /* Doppeltes Zahlenformat                    */
#include "t_cond.h"       /* condition codes fuer ArgumentPruefung     */
#endif
/* ---------------------------------------------------------------------- */

#ifdef ANSI_C
#define const 
#else
#define const
#endif

/* ---------------------------------------------------------------------- */
/* FILE Handle fuer FehlerMeldungen (vom type FILE)                       */
/* ---------------------------------------------------------------------- */
#define MsgOut stderr

/* ---------------------------------------------------------------------- */
/* ExceptionHandle                                                        */
/* ---------------------------------------------------------------------- */
#ifdef LINT_ARGS
extern void msg_exctyp(int typ, const char *name);
/*
extern void msg_exc(const char *name, const ExtReal *arg);
*/
extern void exc_to_a(int exc, char **str);
extern void excfct_to_a(int fct, char **str);
extern int exc_handle_1(int fct, int exc, const  ExtReal *arg,ExtReal *ret);
extern int exc_handle_i1(int fct,int exc, const IExtReal *arg,IExtReal *ret);
extern int default_handle_1(int fct,int exc,const ExtReal *arg,ExtReal *ret);
extern int default_handle_i1(int fct,int exc,const IExtReal *arg,
                             IExtReal *ret);
extern int exc_handle_2(int fct, int exc, const  ExtReal *arg1,
                    const  ExtReal *arg2,  ExtReal *ret);
extern int exc_handle_i2(int fct,int exc, const IExtReal *arg1,
                    const IExtReal *arg2, IExtReal *ret);
extern int default_handle_2(int fct, int exc, const ExtReal *arg1,
                    const ExtReal *arg2,ExtReal *ret);
extern int default_handle_i2(int fct,int exc,const IExtReal *arg1,
                    const IExtReal *arg2,IExtReal *ret);
#else
extern void msg_exctyp(), /* msg_exc(), */ ext_to_a(), excfct_to_a();
extern int exc_handle_1(), exc_handle_i1();
extern int default_handle_1(), default_handle_i1();
extern int exc_handle_2(), exc_handle_i2();
extern int default_handle_2(), default_handle_i2();
#endif /* LINT_ARGS */
/* ---------------------------------------------------------------- */
/* Globale Groesse ArgumentPruefung On oder Off                     */
/* ---------------------------------------------------------------- */
extern char arg_check;
#define On  1
#define Off 0

/* ---------------------------------------------------------------- */
/* Funktionen zur ArgumentPruefung                                  */
/* ---------------------------------------------------------------- */
#ifdef LINT_ARGS
/* --- Ein IntervallArgument --- */
extern int chk_arg_i1(int fct, const IExtReal *arg, IExtReal *res);
/* --- Zwei IntervallArgumente --- */
extern int chk_arg_i2(int fct,
                    const IExtReal *arg1,const IExtReal *arg2,IExtReal *res);

/* --- Ein PunktArgument --- */
extern int chk_arg_1 (int fct, const  ExtReal *arg,  ExtReal *res);
/* --- Zwei PunktArgumente --- */
extern int chk_arg_2 (int fct,
                     const  ExtReal *arg1,const ExtReal *arg2,ExtReal *res);

extern int chk_1 (int fct, const  ExtReal *arg);

#else
extern int chk_arg_i1(), chk_arg_i2(), chk_arg_1(), chk_arg_2(), chk_1();
#endif /* LINT_ARGS */

/* -------------------------------------------------------------- */
/* Macros zur ArgumentPruefung                                    */
/* -------------------------------------------------------------- */
#define ArgCheck1(F,A,R)                                         \
   if(arg_check==On)                                             \
      { if(-1==(check=chk_arg_1(F,A,R))) return 0;               \
        if(NoErr!=check) return check;                           \
      }
#define ArgCheckI1(F,A,R)                                        \
   { if(-1==(check=chk_arg_i1(F,A,R))) return 0;                 \
        if(NoErr!=check) return check;                           \
   }
#define ArgCheck2(F,A1,A2,R)                                     \
   { if(-1==(check=chk_arg_2(F,A1,A2,R))) return 0;              \
        if(NoErr!=check) return check;                           \
   }
#define ArgCheckI2(F,A1,A2,R)                                    \
   { if(-1==(check=chk_arg_i2(F,A1,A2,R))) return 0;             \
        if(NoErr!=check) return check;                           \
   }

/* ----------------------------------------------------------------- */
/* Macro zur Definition von ExtReal Konstanten   ACHTUNG C5.1 Bug!!  */
/* ----------------------------------------------------------------- */
/*  AIX PS/2 C unterstuetzt ## nicht  240890 Baumhof
#define EXTREAL(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9) \
     {{0x##y9,0x##y8,0x##y7,0x##y6,0x##y5,0x##y4,0x##y3,0x##y2,0x##y1,0x##y0}}
*/

#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C
/* #define __CB_SUN4_OS4_HARDWARE */
#endif
#endif

#ifdef __CB_SUN4_OS4_HARDWARE	/* fuer SUN4 ist LSBFIRST=false */
#define EXTREAL(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9) \
	{{  /* Exponent+VZ */	y0, y1,	\
	   /* 64bit Mant  */    (((y2<<1)&0x0fe)+((y3>>7)&0x001)),	\
				(((y3<<1)&0x0fe)+((y4>>7)&0x001)),	\
				(((y4<<1)&0x0fe)+((y5>>7)&0x001)),	\
				(((y5<<1)&0x0fe)+((y6>>7)&0x001)),	\
				(((y6<<1)&0x0fe)+((y7>>7)&0x001)),	\
				(((y7<<1)&0x0fe)+((y8>>7)&0x001)),	\
				(((y8<<1)&0x0fe)+((y9>>7)&0x001)),	\
				((y9<<1)&0x0fe),			\
	    /* Rest 0     */	0, 0, 0, 0, 0, 0  }}
#else
#if LSBFIRST                    /* Intel-Format, 070291 cb */
#define EXTREAL(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9) \
               {{y9,y8,y7,y6,y5,y4,y3,y2,y1,y0}}
#else
#define EXTREAL(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9) \
               {{y0,y1,y2,y3,y4,y5,y6,y7,y8,y9}}
#endif /* LSBFIRST */
#endif

/* ------------------------------------------------------------------ */
/*                                                                    */
/* ------------------------------------------------------------------ */
/* Rueckgabe einer Funktion -(int)ret- enthaelt Except- und RndCtrl-Code */
#define Except (int)0x0fff /* ret&Exception>0 ==> Fehler (s. ieeeexc.h) */
#define RndCtrl ~Except /* ret&RndCtrl Unterscheidung rel Rundungs-Fehler */
/* --- RndCtrl: --- */
#define  UB_Exact   0x2000 /* Upper Bound exact ==> keine Rundung */
#define  LB_Exact   0x1000 /* Lower Bound exact ==> keine Rundung */
/* --- Ende ret --- */

/* --------------------------------------------------------------- */
/* DReal                                                           */
/* --------------------------------------------------------------- */
#ifdef LINT_ARGS
/* --- mul mit doppelt langem Ergebnis --- */
extern int dmulee (const ExtReal *arg1, const ExtReal *arg2,
                         ExtReal *resh, ExtReal *resl);
extern int dmul24dpie(const ExtReal *arg, int two_or_four,
                      ExtReal *resh, ExtReal *resl); /* mul 2/Pi bzw 4/Pi */
#else
extern int dmulee(), dmul24dpie();
#endif /* LINT_ARGS */

/* -------------------------------------------------------------- */
/*                                                                */
/* -------------------------------------------------------------- */
#ifdef LINT_ARGS
/* --- check ExtReal arg --- */
extern int _s_xam(const ExtReal *arg);     /* fxam                 */
extern int _s_chk_invalid(const ExtReal *arg); /* check invalid operation */
extern int chk_extreal(const ExtReal *arg, const int code);

extern int printec (const ExtReal *arg);  /* print ExtReal fuer c prog  */
extern int printecm(const ExtReal *arg);  /* print ExtReal fuer c prog  */

extern int iround_rel(const IExtReal *arg,const ExtReal *eps_rel,
                      IExtReal *res);
extern int iround_abs(const IExtReal *arg,const ExtReal *eps_abs,
                      IExtReal *res);
extern int round_rel(int rnd, const ExtReal *arg,
                     const ExtReal *eps_rel, ExtReal *res);
extern int round_abs(int rnd, const ExtReal *arg,
                     const ExtReal *eps_abs, ExtReal *res);
extern int sround(int rnd, const ExtReal *arg, const ExtReal *eps_rel,
                          const ExtReal *eps_abs, ExtReal *res);

extern int chk_ival(const IExtReal *arg);
#else
extern int _s_xam(), _s_chk_invalid(), chk_extreal();
extern int printec(), printecm();
extern int iround_rel(), iround_abs(), round_rel(), round_abs(), sround();
extern int chk_ival();
#endif /* LINT_ARGS */

/* -------------------------------------------------------------- */
/* Quadrat Wurzel Sqrt                                            */
/* -------------------------------------------------------------- */
#define Sqrt   161
#define Sqrtm1 162
#define ISqrt  261
#ifdef LINT_ARGS
extern int _s_sqrt(const ExtReal *arg, ExtReal *res);
#else
extern int _s_sqrt();
#endif /* LINT_ARGS */
extern const ExtReal MaxArgSqrtm1;
extern const ExtReal EpsSqrt;

/* -------------------------------------------------------------- */
/* trigonometrische Funktionen                                    */
/* -------------------------------------------------------------- */
#define Sin  111
#define Cos  112
#define Tan  113
#define Cot  114
#define ISin 211
#define ICos 212
#define ITan 213
#define ICot 214

#define Asin  121
#define Acos  122
#define Atan  123
#define Acot  124
#define IAsin 221
#define IAcos 222
#define IAtan 223
#define IAcot 224

#define  Period_PiHalf   2
#define  Period_PiQuart  4
#define  J_Init_Sin      0
#define  J_Init_Cos      1
#define  J_Init_Tan      0
#define  J_Init_Cot      2

/* Laenge von 2/Pi= 2.5*ExtMantLen*BitsPerDigit */
#define  LenOfRedConst   160

#ifdef LINT_ARGS
extern int _s_sin(const ExtReal *arg, ExtReal *res);
extern int _s_tan(const ExtReal *arg, ExtReal *resn, ExtReal *resd);
extern int gza_trg(const ExtReal *arg, int jinit, int periode,
                   DReal *v, ExtReal *j, int *jmod4);
extern int red_trg(const DReal *v, const ExtReal *j, int jmod4, ExtReal *x);
extern int sincos(const ExtReal *t, ExtReal *res);
extern int tancot(const ExtReal *t, int jmod4, ExtReal *res);
extern int isincos(const IDReal *v, const IExtReal *j,
                   int jmod4u, int jmod4l, IExtReal *res);
#else
extern int _s_sin(), _s_tan(), gza_trg(), red_trg();
extern int sincos(), tancot(), isincos();
#endif /* LINT_ARGS */
extern const ExtReal EpsSin;
extern const ExtReal EpsCos;
extern const ExtReal EpsTan;
extern const ExtReal EpsCot;

/* --------------------------------------------------------------- */
/* --- trigonometrische Umkehr-Funktionen --- */
#ifdef LINT_ARGS
extern int _s_atan(const ExtReal *arg, ExtReal *res);

extern int asinviatan(const ExtReal *arg, ExtReal *res);
extern int acosviatan(const ExtReal *arg, ExtReal *res);
#else
extern int _s_atan(), asinviatan(), acosviatan();
#endif /* LINT_ARGS */
extern const ExtReal EpsASin;
extern const ExtReal EpsACos;
extern const ExtReal EpsATan;
extern const ExtReal EpsACot;

/* ---------------------------------------------------------------- */
/* --- Exp --- */
#define  Exp   150
#define  Expm1 151
#define IExp   250
extern const ExtReal EpsExp;
#ifdef LINT_ARGS
extern int _s_2xm1 (const ExtReal *arg, ExtReal *res);
extern int scalee (const ExtReal *arg, const ExtReal *intscalfakt,
                   ExtReal *res);
extern int scaliee(const ExtReal *arg, int scalfakt, ExtReal *res);
#else
extern int _s_2xm1(), scalee(), scaliee();
#endif /* LINT_ARGS */

extern const ExtReal MaxArgExp; /* Maximales Argument expee vor OVERFLOW  */
extern const ExtReal MinArgExp; /* Minimales Argument expee vor UNDERFLOW */

/* --- Hyperbolische --- */
#define  Sinh 131
#define ISinh 231
#define  Cosh 132
#define ICosh 232
#define  Tanh 133
#define ITanh 233
#define  Coth 134
#define ICoth 234
extern const ExtReal EpsSinh;
extern const ExtReal EpsCosh;
extern const ExtReal EpsTanh;
extern const ExtReal EpsCoth;
extern const ExtReal MaxArgHyp; /* Maximales Argument */
extern const ExtReal MinArgHyp; /* Minimales Argument */

/* --- Area Hyperbolische --- */
#define  Asinh 141
#define IAsinh 241
#define  Acosh 142
#define IAcosh 242
#define  Atanh 143
#define IAtanh 243
#define  Acoth 144
#define IAcoth 244
extern const ExtReal EpsASinh;
extern const ExtReal EpsACosh;
extern const ExtReal EpsATanh;
extern const ExtReal EpsACoth;

/* --- Log --- */
#define  Ln    152
#define ILn    252
#define  Lnp1  153
extern const ExtReal EpsLn;
extern const ExtReal EpsLnRel;
extern const ExtReal EpsLnAbs;
extern const ExtReal MaxArgLnp1; /* Maximales Argument lnp1 DOMAIN Coproc */
#ifdef LINT_ARGS
extern int _s_ln  (const ExtReal *arg, ExtReal *res, int *rndtyp);
extern int _s_lnp1(const ExtReal *arg, ExtReal *res);
#else
extern int _s_ln(), _s_lnp1();
#endif /* LINT_ARGS */

/* ------------------------------------------------------------- */
/* --- loga, pow --- */
#define  Pow   160
#define  IPow  260
#ifdef LINT_ARGS
extern int powsub(const ExtReal *bas, const ExtReal *exp,
                  ExtReal *res, ExtReal *ri);
#else
extern int powsub();
#endif /* LINT_ARGS */
extern const ExtReal EpsPow;

/* ------------------------------------------------------------- */
/* sonstige Funktionen                                           */
/* ------------------------------------------------------------- */
#define EToI  180
#define EToL  181
#define LToE  182

#define Div   190
#define IDiv  290

#ifdef LINT_ARGS
extern int xtrexpe(const ExtReal *arg, int *e);

extern int _s_etoi(const ExtReal *arg, int *res);
extern int _s_etol(const ExtReal *arg, LongReal *res);
extern int _s_ltoe(const LongReal *arg, ExtReal *res);
extern int _s_etoie(const ExtReal *arg, ExtReal *res);
#else
extern int xtrexpe();
extern int _s_etoi(), _s_etol(), _s_ltoe(), _s_etoie();
#endif /* LINT_ARGS */

extern const ExtReal IntMax;      /* = INT_MAX = 32767 = (int)0x7fff   */
extern const ExtReal IntMin;      /* = INT_MIN = (int)0x8000           */
extern const ExtReal LongRealMax; /* = (double)DBL_MAX in float.h      */
extern const ExtReal LongRealMin; /* = (double)DBL_MIN in float.h      */
extern const ExtReal LongRealDenormMin; /* kleinste pos. Zahl          */

/* ------------------------------------------------------------------- */
/* sonstige Funktionen                                                 */
/* ------------------------------------------------------------------- */

#ifdef LINT_ARGS
/*--- modulo --- */
extern int mod2e(const ExtReal *arg);
extern int mod4e(const ExtReal *arg);

/* Abbruch */
extern void ieee_abortr1(int rc, void *arg);
extern void ieee_abortr2(int rc, void *arg1, void *arg2);
extern void ieee_aborti1(int rc, void *ai);
extern void ieee_aborti2(int rc, void *ai, void *bi);

extern int cond_to_exc(int condition);
#else
extern int mod2e(), mod4e();
extern void ieee_abortr1(), ieee_abortr2();
extern void ieee_aborti1(), ieee_aborti2();
extern int cond_to_exc();
#endif /* LINT_ARGS */

/* ---------------------------------------------------------------- */
/*                                                                  */
/* ---------------------------------------------------------------- */

/* --- macros ----------------------------------------------------- */
/* !!!!!!!!!!!!!!!!!!!!! dirty but efficient !!!!!!!!!!!!!!!!!!!!!! */

/*#define Bound(rnd, arg) (ExtReal*)arg+(((unsigned)rnd)>>15)*/

/* returns arg->l if rnd == DOWN, arg->u otherwise
   assumes u to be first element in struct IExtReal */
/* !!!!!!!!!!!!!!!!!!!!! dirty but efficient !!!!!!!!!!!!!!!!!!!!!!!! */

/* setzen deklaration "int ret;" voraus. "ret" ist Fehler-return-code */
/* ret == 0: fehlerfrei, andere noch nicht festgelegt */

#define SQRT(a,b) if(ret=isqrtee(a,b)) return ret;

#define SIN(a,b)  if(ret=isinee(a,b))  return ret;
#define COS(a,b)  if(ret=icosee(a,b))  return ret;
#define TAN(a,b)  if(ret=itanee(a,b))  return ret;
#define COT(a,b)  if(ret=icotee(a,b))  return ret;

#define ASIN(a,b) if(ret=iasinee(a,b)) return ret;
#define ACOS(a,b) if(ret=iacosee(a,b)) return ret;
#define ATAN(a,b) if(ret=iatanee(a,b)) return ret;
#define ACOT(a,b) if(ret=iacotee(a,b)) return ret;

#define LN(a,b)    if(ret=ilnee    (a,b)) return ret;
#define LNP1(a,b)  if(ret=ilnp1ee  (a,b)) return ret;
#define EXP(a,b)   if(ret=iexpee   (a,b)) return ret;
#define SINH(a,b)  if(ret=isinhee  (a,b)) return ret;
#define COSH(a,b)  if(ret=icoshee  (a,b)) return ret;
#define TANH(a,b)  if(ret=itanhee  (a,b)) return ret;
#define COTH(a,b)  if(ret=icothee  (a,b)) return ret;

#define ASINH(a,b) if(ret=iasinhee(a,b)) return ret;
#define ACOSH(a,b) if(ret=iacoshee(a,b)) return ret;
#define ATANH(a,b) if(ret=iatanhee(a,b)) return ret;
#define ACOTH(a,b) if(ret=iacothee(a,b)) return ret;

#define POW(a,b,c)  if(ret=ipowee(a,b,c))  return ret;
#define LOGA(a,b,c) if(ret=ilogaee(a,b,c)) return ret;

/* ------------------------------------------------------------- */

#ifdef ANSI_C
#else
#undef const
/*
#undef extern
*/
#endif

/* ---------------------------------------------------------------- */
/* Routines to force linkage                                        */
/* ---------------------------------------------------------------- */
#if VAX_VMS_C
#ifdef LINT_ARGS
extern void t_cnst(void);
extern void t_ecst(void);
extern void t_glbl(void);
#else
extern void t_glbl(), t_cnst(), t_ecst();
#endif
#endif

#ifdef LINT_ARGS
#ifdef ANSI_C
int round_ln(int rnd, int rr, const ExtReal *arg, ExtReal *res);
#else  /* NOT ANSI_C */
int round_ln(int rnd, int rr, ExtReal *arg, ExtReal *res);
#endif /* ANSI_C */
#else  /* NOT LINT_ARGS */
int round_ln();
#endif /* LINT_ARGS */

#ifdef ANSI_C
#ifdef LINT_ARGS
int t_10ex(const ExtReal *arg, ExtReal *res);
#else
int t_10ex();
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int t_10ex(ExtReal *arg, ExtReal *res);
#else
int t_10ex();
#endif /* LINT_ARGS */
#endif /* ANSI_C */
 
#ifdef ANSI_C
#ifdef LINT_ARGS
int t_2exp(const ExtReal *arg, ExtReal *res);
#else
int t_2exp();
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int t_2exp(ExtReal *arg, ExtReal *res);
#else
int t_2exp();
#endif /* LINT_ARGS */
#endif /* ANSI_C */
/* -------------------------------------------------------------- */





