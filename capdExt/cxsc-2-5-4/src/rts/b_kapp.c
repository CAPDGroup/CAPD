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

/* CVS $Id: b_kapp.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*      Filename        : b_kapp.c                              */
/*                                                              */
/*      Description     : Software approximations to standard   */
/*                        functions. These routines are used    */
/*                        as alternative for definition in      */
/*                        macros B_SQRT,B_LOG_,B_ASIN,B_ATAN.   */
/*                        Refer to b_lari.h.                    */
/*                                                              */
/****************************************************************/
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif

extern a_real *r_1o2_;
extern a_real *r_zero;
extern a_real *r_mone;
extern a_real *r_one_;
extern a_real *r_two_;
extern a_real *r_pio2;

#define a__log  b_klog
#define a__atan b_katn
#define a__sqrt b_ksqt
#define a__asin b_kasn

typedef union k_type { a_real r; a_btyp u[D_U_RATIO]; } k_type;

static int not_initialized = 1;
static k_type log1, log2, log3, log4, log5, log6, log7, log8, log9;
static k_type atn1, atn2, atn3, atn4, atn5;
static k_type sqt1, sqt2, sqt3, sqt4, sqt5, sqt6, sqt7;
static k_type asn1, asn2, asn3, asn4, asn5, asn6;

#ifdef LINT_ARGS
static void k_init(void)
#else
static void k_init()
#endif
{
   not_initialized = 0;

   /* --- decimal constants:                                --- */
   log1.u[B_HPART] = 0x3FEFAE14l;
   log1.u[B_LPART] = 0x7AE147AEl;   /* 9.900000000000000E-001  1 */
   log2.u[B_HPART] = 0x3FF0CCCCl;
   log2.u[B_LPART] = 0xCCCCCCCDl;   /* 1.050000000000000E+000  2 */
   log3.u[B_HPART] = 0x3FEFFF2El;
   log3.u[B_LPART] = 0x48E8A71El;   /* 9.999000000000000E-001  3 */
   log4.u[B_HPART] = 0x400558C0l;
   log4.u[B_LPART] = 0x3302D91Bl;   /* 2.668335341000000E+000  4 */
   log5.u[B_HPART] = 0x3FEF9B9El;
   log5.u[B_LPART] = 0x86B24407l;   /* 9.877464896000000E-001  5 */
   log6.u[B_HPART] = 0x3FFAE3B8l;
   log6.u[B_LPART] = 0xF44A36A7l;   /* 1.680596308000000E+000  6 */
   log7.u[B_HPART] = 0x4006CBD5l;
   log7.u[B_LPART] = 0xD3644933l;   /* 2.849528934000000E+000  7 */
   log8.u[B_HPART] = 0x3FDFFEFFl;
   log8.u[B_LPART] = 0x5BABA782l;   /* 4.999388118000000E-001  8 */
   log9.u[B_HPART] = 0x3FE62E42l;
   log9.u[B_LPART] = 0xFEFFBB3Cl;   /* 6.931471806000000E-001  9 */

   atn1.u[B_HPART] = 0x3FCED634l;
   atn1.u[B_LPART] = 0x113D8BC5l;   /* 2.409119686300000E-001  10 */
   atn2.u[B_HPART] = 0x400E47E8l;
   atn2.u[B_LPART] = 0xECBC4C0El;   /* 3.785112237450000E+000  11 */
   atn3.u[B_HPART] = 0x4016B552l;
   atn3.u[B_LPART] = 0x630D6339l;   /* 5.677072093670000E+000  12 */
   atn4.u[B_HPART] = 0x4016B58Al;
   atn4.u[B_LPART] = 0x4FF4BCDAl;   /* 5.677285432160000E+000  13 */
   atn5.u[B_HPART] = 0x4016B553l;
   atn5.u[B_LPART] = 0x144F6F60l;   /* 5.677074735020000E+000  14 */

   sqt1.u[B_HPART] = 0x4010F1BBl;
   sqt1.u[B_LPART] = 0x24631EA4l;   /* 4.236065453100000E+000  15 */
   sqt2.u[B_HPART] = 0x40180D84l;
   sqt2.u[B_LPART] = 0xB4AAF905l;   /* 6.013201544700000E+000  16 */
   sqt3.u[B_HPART] = 0x3FDAE89Fl;
   sqt3.u[B_LPART] = 0xA90BA82El;   /* 4.204482222400000E-001  17 */
   sqt4.u[B_HPART] = 0x401C9A90l;
   sqt4.u[B_LPART] = 0x8AE1A1C4l;   /* 7.150942010900000E+000  18 */
   sqt5.u[B_HPART] = 0x40042675l;
   sqt5.u[B_LPART] = 0xF0BF67ACl;   /* 2.518779641000000E+000  19 */
   sqt6.u[B_HPART] = 0x3FF6A09El;
   sqt6.u[B_LPART] = 0x6665983El;   /* 1.414213562000000E+000  20 */
   sqt7.u[B_HPART] = 0x3FE6A09El;
   sqt7.u[B_LPART] = 0x669C91FEl;   /* 7.071067814000001E-001  21 */

   asn1.u[B_HPART] = 0x3FE6A11Cl;
   asn1.u[B_LPART] = 0xB039EF0Fl;   /* 7.071670000000000E-001  22 */
   asn2.u[B_HPART] = 0x3FE0BA37l;
   asn2.u[B_LPART] = 0x7BB63C8Fl;   /* 5.227315346000000E-001  23 */
   asn3.u[B_HPART] = 0x4010DE56l;
   asn3.u[B_LPART] = 0x376F2D09l;   /* 4.217125765000000E+000  24 */
   asn4.u[B_HPART] = 0x40130944l;
   asn4.u[B_LPART] = 0x7FB50CAAl;   /* 4.759050364900000E+000  25 */
   asn5.u[B_HPART] = 0x40140A7Al;
   asn5.u[B_LPART] = 0x36B073D9l;   /* 5.010231833000000E+000  26 */
   asn6.u[B_HPART] = 0x40130944l;
   asn6.u[B_LPART] = 0x4292304Bl;   /* 4.759049453900000E+000  27 */
   return;
}

#ifdef LINT_ARGS
a_real a__log( a_real x )
#else
a_real a__log( x )
a_real x;
#endif
{
  a_real m, y;
  int ex;  /* Exponent zur Basis 2 */
  
  if (not_initialized) k_init();

  if ( r_lt(log1.r,x) && r_lt(x,log2.r) ) {
    R_ASSIGN(y,r_addd(x,*r_mone)) ;
    return ( r_divd(r_addd(y,y),r_addd(*r_two_,y)) );
  } 
  
  R_ASSIGN(m,x);
  ex= 0;
  if ( r_lt(m, *r_1o2_) ) { /* Exponentenberechnung sollte ueber    */
    for (ex= -1; ;ex--) {    /* Exponentenzugriff erfolgen!       */
      R_ASSIGN(m, r_addd(m,m));
      /* m *= 2.0; */
      if ( r_gt(m, *r_1o2_) ) break;
    } 
  }  
  else if ( r_gt(m, *r_one_) ) {
    for (ex=1; ;ex++) {
      R_ASSIGN(m, r_muld(m, *r_1o2_));
      /* m *= 0.5;  */
      if ( r_lt(m, *r_one_) ) break;
    }      
  }
  /* Die Mantisse 0.5 <= m < 1 ist bzgl. Basis 2 normalisiert */ 
  /* ex ist der zugehoerige Exponent                          */ 
  if ( r_gt(m, log3.r) ) {
    R_ASSIGN(y, r_subd(m, *r_one_));
    /* y= m - 1; */
  }  
  else {
    R_ASSIGN(y, r_divd(
       r_subd(r_muld(r_subd(r_muld(log4.r,m), log5.r), m), log6.r),
       r_addd(r_muld(r_addd(m, log7.r), m),  log8.r) ) );   /* 2622*/
  }
  return r_addd(y, r_muld(r_flot((a_intg)ex), log9.r)) ;
}

#ifdef LINT_ARGS
a_real a__atan( a_real x )
#else
a_real a__atan( x )
a_real x;
#endif
{
  a_real y, yy;
  
  if (not_initialized) k_init();

  R_ASSIGN(y, x);
  if ( r_lt(x, *r_zero) ) { R_ASSIGN(y, r_umin(y)); }
   
  if ( r_gt(y, *r_one_) ) {
    R_ASSIGN(y, r_divd(*r_one_, y));
    /* y= 1.0 / y; */
  }
  R_ASSIGN(yy, r_muld(y,y) );
  R_ASSIGN(y, r_divd(
   r_muld(y, r_addd(r_muld(r_addd(r_muld(atn1.r, yy), atn2.r), yy),atn3.r)),
   r_addd(r_muld(r_addd(yy,atn4.r), yy), atn5.r) ) );       /*5091*/
              
  if ( r_lt(x, *r_mone) || r_gt(x, *r_one_) ) {
    R_ASSIGN(y, r_subd(*r_pio2, y));
    /*y= 1.570796327 - y; */
  }
  if ( r_lt(x, *r_zero) ) { R_ASSIGN(y, r_umin(y)); }
  return y;
}  


#ifdef LINT_ARGS
a_real a__sqrt( a_real x )
#else
a_real a__sqrt( x )
a_real x;
#endif
{
  a_real m, y;
  int ex, ex2, i;
  
  if (not_initialized) k_init();

  R_ASSIGN(m, x);
  ex= 0;

  if ( r_eq(x, *r_zero) || r_eq(x, *r_one_) ) { return(x); } /* !!!
*/
  if ( r_lt(m,*r_1o2_) ) {    /* Exponentenberechnung sollte ueber */
    for (ex= -1; ;ex--) {    /* Exponentenzugriff erfolgen!       */
      R_ASSIGN(m, r_addd(m,m));
      /* m *= 2.0; */
      if ( r_gt(m, *r_1o2_) ) break;
    } 
  }  
  else if ( r_ge(m, *r_one_) ) {
    for (ex=1; ;ex++) {
      R_ASSIGN(m, r_muld(m, *r_1o2_));
      /* m *= 0.5; */
      if ( r_lt(m,*r_one_) ) break;
    }      
  }
  /* Die Mantisse 0.5 <= m < 1 ist bzgl. Basis 2 normalisiert */ 
  /* ex ist der zugehoerige gerade oder ungerade Exponent     */  
  R_ASSIGN(y, r_divd(
     r_addd(r_muld(r_addd(r_muld(sqt1.r, m), sqt2.r), m), sqt3.r),
     r_addd(r_muld(r_addd(m, sqt4.r), m), sqt5.r) ) );
     
  ex2 = ex/2;
  if ( ex > 0 ) {
    for (i=1; i<=ex2; i++) {
      R_ASSIGN(y, r_addd(y,y));
      /* y*= 2; */
    }
    if ( 2*ex2 != ex ) {
      R_ASSIGN(y, r_muld(y,sqt6.r));
      /* y*= 1.414213562;     y = y * sqrt(2)    */
    }    
  }  
  else if (ex < 0 ) {
    for (i= -1; i>=ex2; i--) {
      R_ASSIGN(y, r_muld(y, *r_1o2_));
      /* y*= 0.5; */
    }
    if ( 2*ex2 != ex ) {
      R_ASSIGN(y, r_muld(y, sqt7.r));
      /* y*= 0.7071067814;    y = y*( 1/sqrt(2) ) */
    }    
  }    
  
  return y;
}  

#ifdef LINT_ARGS
a_real a__asin( a_real x )
#else
a_real a__asin( x )
a_real x;
#endif
{
  a_real y, yy;
  int flag;

  if (not_initialized) k_init();

  flag= 0;  /* gibt an, ob abs(x) > 1/sqrt(2) = 0.707..  */
  R_ASSIGN(y, x);
  if ( r_lt(y, *r_zero) ) R_ASSIGN(y, r_umin(y));
  if ( r_gt(y, asn1.r) ) {
    flag = 1;
    R_ASSIGN(y, a__sqrt( r_muld(r_subd(*r_one_,y), r_addd(*r_one_,y)) ));
  }

  R_ASSIGN(yy, r_muld(y, y));
  R_ASSIGN( y, r_divd(
   r_muld(y,
     r_addd(r_muld(r_subd(r_muld(asn2.r, yy), asn3.r), yy), asn4.r)),
   r_addd(r_muld(r_subd(yy, asn5.r), yy), asn6.r) ) );            /*4723*/

  if ( flag ) {
    R_ASSIGN(y, r_subd(*r_pio2, y));
    /*y= 1.570796327 - y;   arcsin|x| = pi/2 - arcsin(sqrt(1-x*x)) */
  }
  if ( r_lt(x, *r_zero) ) R_ASSIGN(y, r_umin(y));
  return y;
}

/*
main( )
{ 
  a_real x;
  extern a_real log(),sqrt(), atan(), asin();
  for(;;)
  {
    printf("log(x) wird berechnet! x=?\n");
    scanf("%le", &x);
    printf("Eingelesen wurde %25.19le \n", x);
    printf("log(x)          %25.19le \n", a__log(x));
    printf("log(x) aus LIB  %25.19le \n\n", log(x));
    printf("sqrt(x)          %25.19le \n", a__sqrt(x));
    printf("sqrt(x) aus LIB  %25.19le \n\n", sqrt(x));
    printf("atan(x)          %25.19le \n", a__atan(x));
    printf("atan(x) aus LIB  %25.19le \n\n", atan(x));
    printf("asin(x)          %25.19le \n", a__asin(x));
    printf("asin(x) aus LIB  %25.19le \n\n", asin(x));
  }
  return 0; }  
*/





