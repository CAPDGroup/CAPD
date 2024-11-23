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

/* CVS $Id: p88rts.h,v 1.29 2014/01/30 17:24:11 cxsc Exp $ */

/***********************************************************/
/*                CONFIGURE P88RTS.H here !!!              */
/***********************************************************/
/* define   switches according to your system defaults     */
/* USING.. valid in this file                              */
/***********************************************************/

#include "cxscconf.h"
#ifdef SUN4_FORTE
#undef SUN4_GNU_C
#define SUN4_GNU_C 1  
#endif

/* new style prototyping compilers must have 1, old one 0 */
/* #define USING_LINT_ARGS 1 */
#if SUN4_CPP_C+SUN4_GNU_C+SUN4_GNU_C+SUN_OS5_GNU_C+IBM_LINUX_C+DEC_ALPHA_C+CXSC_PPC64
#define USING_LINT_ARGS 1
#else
#define USING_LINT_ARGS 0
#endif

/* use 1 for ATARI, VAX VMS , 0 for others */
#define USING_AREALSTRU 0

/* only 64 bit systems must have 1, others 0 */
#if DEC_ALPHA_C+GNU_X86_64+CXSC_PPC64
#define USING_64_BITSYS 1
#else
#define USING_64_BITSYS 0 
#endif

/***********************************************************/
/*              end of configuring section                 */
/***********************************************************/
#define t_ep10 t_ex10
                        /* provisorisch */
#if USING_LINT_ARGS
/* new-style prototyping */
#define LINT_ARGS 
#endif
#undef USING_LINT_ARGS

#ifndef LINT_ARGS
#define IBMva
                 /* old style variable argument lists */
#endif
                 /* must be defined for AIX PS2 */
/* #define IBMva */
                 /* IBM style variable argument lists */


/*--------------------------------------------------------------*/
/* PROTOTYPE manages the difficulties in prototyping of K & R   */
/* and new style ANSI compilers                                 */
/*--------------------------------------------------------------*/
#ifndef _PROTOTYPE
#ifdef LINT_ARGS
#define _PROTOTYPE(function,params) function params
#else 
#define _PROTOTYPE(function,params) function()
#endif
#endif

#define a_RETC (a_intg)0  /* value returned by main (default) */

#define a_local static
		  /* indicates denested routines */

/*---------------------------------------------------*/
/*  global pointer variable, used by var-traversing  */
/*---------------------------------------------------*/
extern char *tp_000;
extern char *tp_001;


/*----------------------------------------------*/
/* Number of characters forming a set value     */
/*----------------------------------------------*/
#define s_SIZE          32


/*--------------------------------------*/
/*   standard data types of Pascal-SC   */
/*--------------------------------------*/
typedef unsigned char a_byte ;

#ifdef LINT_ARGS
typedef void * a_VOID ;
#else
typedef char * a_VOID ;
#endif

#if USING_64_BITSYS
typedef int             a_intg ;
#if USING_AREALSTRU
typedef struct { unsigned int a[ 2 ]; } a_real;
#else
typedef double          a_real ;   
#endif
typedef unsigned char   a_char ;
typedef a_byte          a_bool ;
typedef unsigned int    a_btyp ;
typedef struct a_long { unsigned aa[ 3 ]; } a_long; 
/* end of using 64 bit */
#else
typedef long int        a_intg ;   /* adjust a_imax, if you modify a_intg */
#if USING_AREALSTRU
typedef struct a_real { a_byte a[ 8 ]; } a_real; 
#else
typedef double          a_real ;   
#endif
typedef unsigned char   a_char ;
typedef a_byte          a_bool ;
typedef unsigned long   a_btyp ;
typedef struct a_long { a_byte aa[ 12 ]; } a_long; 
#endif
#undef USING_64_BITSYS
#undef USING_AREALSTRU

typedef struct s_trng { char *ptr;
			size_t alen;
			size_t clen;
			unsigned int fix  : 1;
			unsigned int suba : 1;
			unsigned int tmp  : 1; } s_trng;
typedef a_char          s_etof [s_SIZE];
typedef a_btyp   *      d_otpr ;
typedef struct a_cmpx { a_real RE , IM; } a_cmpx ;
typedef struct a_intv { a_real INF, SUP;} a_intv ;
typedef struct a_cinv { a_intv RE , IM ;} a_cinv ;
typedef struct a_ilng { a_long INF, SUP;} a_ilng ;

typedef struct d_otpc { d_otpr RE , IM; } d_otpc ;
typedef struct d_otpi { d_otpr INF, SUP;} d_otpi ;
typedef struct d_otpz { d_otpi RE , IM; } d_otpz ;


/*----------------------------------------*/
/*   error handling, exception handling   */
/*----------------------------------------*/
extern int e_line ;
#ifdef LINT_ARGS
#ifdef IBMva
extern void e_trap () ;
#else
extern void e_trap (a_btyp code, int e_argc, ...) ;
#endif
#else
extern void e_trap () ;
#endif
_PROTOTYPE(extern void e_push, (char * functionname, char * filename)) ;
_PROTOTYPE(extern void e_popp, (void)) ;
_PROTOTYPE(extern void a_exit, (a_intg)) ;    /* closes files and exits */


/*-------------------------------------*/
/*   standard constants of Pascal-SC   */
/*-------------------------------------*/
#define a_true 1
#define a_flse 0
#define a_imax 0x7fffffff      /* = MAXINT */

/*--------------------------*/
/* Real constant conversion */
/*--------------------------*/

_PROTOTYPE(extern a_real r_cnsd, (char*)) ;
_PROTOTYPE(extern a_real r_cnst, (char*)) ;
_PROTOTYPE(extern a_real r_cnsu, (char*)) ;


/*-----------------------*/
/*   Memory management   */
/*-----------------------*/
#define a_asgn(s1,s2,c)    memcpy((s1),(s2),(c))
#define a_clrm(s,c)        memset((s),0,(c))

_PROTOTYPE(extern a_VOID  a_lloc, (size_t)) ;   /* malloc + errorhandling */
_PROTOTYPE(extern void    a_free, (char **)) ;  /* free + errorhandling   */
_PROTOTYPE(extern void    a_2psh, (a_VOID,a_VOID)) ;   /* push pointers */
_PROTOTYPE(extern void    a_cpsh, (a_VOID,a_VOID)) ;   /* push pointers */
_PROTOTYPE(extern void    a_2pop, (void)) ;         /* set/reset display  */
				       /* and pop stack  */

/*--------------------*/
/*   Set Operations   */
/*--------------------*/
_PROTOTYPE(extern a_bool  s_etin, (a_intg, s_etof)) ;     /* operator "in" */
/* empty set constructor [] */
_PROTOTYPE(extern a_VOID /* s_etof */  s_zero, (s_etof)) ;

#ifdef LINT_ARGS
#ifdef IBMva
extern a_VOID /* s_etof */  s_cons () ;         /* set constructor */
#else
extern a_VOID   s_cons (s_etof res, a_char * ctrl,...);  /* set constructor */
#endif
#else
extern a_VOID  s_cons () ;              /* set constructor */
#endif

/* set = set + set */
_PROTOTYPE(extern a_VOID /* s_etof */ s_add, (s_etof, s_etof, s_etof)); 
/* set = set - set */
_PROTOTYPE(extern a_VOID /* s_etof */ s_sub, (s_etof, s_etof, s_etof)); 
/* set = set * set */
_PROTOTYPE(extern a_VOID /* s_etof */ s_mul, (s_etof, s_etof, s_etof)); 

_PROTOTYPE(extern a_bool  s_eteq, (s_etof, s_etof));       /* set == set */
_PROTOTYPE(extern a_bool  s_etne, (s_etof, s_etof));       /* set != set */
_PROTOTYPE(extern a_bool  s_etge, (s_etof, s_etof));       /* set >= set */
_PROTOTYPE(extern a_bool  s_etgt, (s_etof, s_etof));       /* set >  set */
_PROTOTYPE(extern a_bool  s_etle, (s_etof, s_etof));       /* set <= set */
_PROTOTYPE(extern a_bool  s_etlt, (s_etof, s_etof));       /* set <  set */



/*-----------------------*/
/*   String Operations   */
/*-----------------------*/
#define s_inxn( s, i)  (*((s).ptr + (i) - 1))  /* index access without check */
#define s_asta( d, s, l)  s_asgn ((s_trng*)(d), s_stat ((s), (l)))
			      /* assigns a static array to a dynamic array */

    /* comparisons of mixed static & dynamic strings are impl. as macros */
#define s_saeq(x,y,n)  (s_sseq ((x), s_stat((y), (n))))
#define s_aseq(x,n,y)  (s_sseq (s_stat((x), (n)), (y)))
#define s_sane(x,y,n)  (s_ssne ((x), s_stat((y), (n))))
#define s_asne(x,n,y)  (s_ssne (s_stat((x), (n)), (y)))
#define s_sage(x,y,n)  (s_ssge ((x), s_stat((y), (n))))
#define s_asge(x,n,y)  (s_ssge (s_stat((x), (n)), (y)))
#define s_sagt(x,y,n)  (s_ssgt ((x), s_stat((y), (n))))
#define s_asgt(x,n,y)  (s_ssgt (s_stat((x), (n)), (y)))
#define s_sale(x,y,n)  (s_ssle ((x), s_stat((y), (n))))
#define s_asle(x,n,y)  (s_ssle (s_stat((x), (n)), (y)))
#define s_salt(x,y,n)  (s_sslt ((x), s_stat((y), (n))))
#define s_aslt(x,n,y)  (s_sslt (s_stat((x), (n)), (y)))

_PROTOTYPE(extern a_intg s_ixch, (a_intg index, size_t length)) ;
/* index access with check */
_PROTOTYPE(extern a_char * s_inxc, (s_trng, a_intg) );        
    /*   (*((s).ptr + s_ixch((i),(s).alen)))   */
    /* s_inxc must not be a macro, because the string variable s may */
    /*   contain indices with function calls.  */
    /* s_inxc must return a pointer, because it must be a left-value. */
/* index access with check */
_PROTOTYPE(extern a_char * s_ixcn, (s_trng*, a_intg)) ;        
    /* s_ixcn must have a var parameter, because "resize" might be necessary */

_PROTOTYPE(extern a_bool  s_aaeq, (a_char[], a_intg, a_char[], a_intg ));
_PROTOTYPE(extern a_bool  s_aceq, (a_char[], a_intg, a_char   ));
_PROTOTYPE(extern a_bool  s_caeq, (a_char          , a_char[], a_intg ));
_PROTOTYPE(extern a_bool  s_aane, (a_char[], a_intg, a_char[], a_intg ));
_PROTOTYPE(extern a_bool  s_acne, (a_char[], a_intg, a_char   ));
_PROTOTYPE(extern a_bool  s_cane, (a_char          , a_char[], a_intg ));
_PROTOTYPE(extern a_bool  s_aage, (a_char[], a_intg, a_char[], a_intg ));
_PROTOTYPE(extern a_bool  s_acge, (a_char[], a_intg, a_char   ));
_PROTOTYPE(extern a_bool  s_cage, (a_char          , a_char[], a_intg ));
_PROTOTYPE(extern a_bool  s_aagt, (a_char[], a_intg, a_char[], a_intg ));
_PROTOTYPE(extern a_bool  s_acgt, (a_char[], a_intg, a_char   ));
_PROTOTYPE(extern a_bool  s_cagt, (a_char          , a_char[], a_intg ));
_PROTOTYPE(extern a_bool  s_aale, (a_char[], a_intg, a_char[], a_intg ));
_PROTOTYPE(extern a_bool  s_acle, (a_char[], a_intg, a_char   ));
_PROTOTYPE(extern a_bool  s_cale, (a_char          , a_char[], a_intg ));
_PROTOTYPE(extern a_bool  s_aalt, (a_char[], a_intg, a_char[], a_intg ));
_PROTOTYPE(extern a_bool  s_aclt, (a_char[], a_intg, a_char   ));
_PROTOTYPE(extern a_bool  s_calt, (a_char          , a_char[], a_intg ));

_PROTOTYPE(extern a_bool  s_sseq, (s_trng, s_trng));
_PROTOTYPE(extern a_bool  s_sceq, (s_trng, a_char));
_PROTOTYPE(extern a_bool  s_cseq, (a_char, s_trng ));
_PROTOTYPE(extern a_bool  s_ssne, (s_trng, s_trng ));
_PROTOTYPE(extern a_bool  s_scne, (s_trng, a_char));
_PROTOTYPE(extern a_bool  s_csne, (a_char, s_trng ));
_PROTOTYPE(extern a_bool  s_ssge, (s_trng, s_trng ));
_PROTOTYPE(extern a_bool  s_scge, (s_trng, a_char));
_PROTOTYPE(extern a_bool  s_csge, (a_char, s_trng ));
_PROTOTYPE(extern a_bool  s_ssgt, (s_trng, s_trng ));
_PROTOTYPE(extern a_bool  s_scgt, (s_trng, a_char));
_PROTOTYPE(extern a_bool  s_csgt, (a_char, s_trng ));
_PROTOTYPE(extern a_bool  s_ssle, (s_trng, s_trng ));
_PROTOTYPE(extern a_bool  s_scle, (s_trng, a_char));
_PROTOTYPE(extern a_bool  s_csle, (a_char, s_trng ));
_PROTOTYPE(extern a_bool  s_sslt, (s_trng, s_trng ));
_PROTOTYPE(extern a_bool  s_sclt, (s_trng, a_char));
_PROTOTYPE(extern a_bool  s_cslt, (a_char, s_trng ));
_PROTOTYPE(extern void    s_vlcp, (s_trng *));
_PROTOTYPE(extern void    s_utmp, (s_trng *));
_PROTOTYPE(extern void    s_temp, (s_trng *));
_PROTOTYPE(extern s_trng  s_stat, (a_char [], a_intg));
_PROTOTYPE(extern void    s_init, (s_trng *,size_t));
_PROTOTYPE(extern void    s_free, (s_trng *));
_PROTOTYPE(extern s_trng  s_conc, (s_trng, s_trng));
_PROTOTYPE(extern s_trng  s_char, (a_char));
_PROTOTYPE(extern void    s_asgn, (s_trng *,s_trng));

_PROTOTYPE(extern s_trng  s_date, (s_trng));
/*-------------------------------------*/
/*   Scalar product for static arrays  */
/*-------------------------------------*/
_PROTOTYPE(extern a_real  r_scps, (a_real x[], a_real y[], a_intg size, 
							a_intg rdmode)) ;

/*------------------------------*/
/*  dotprecision routines       */
/*------------------------------*/
_PROTOTYPE(extern void d_init, (d_otpr * )) ;
_PROTOTYPE(extern void d_vlcp, (d_otpr * )) ;
_PROTOTYPE(extern void d_free, (d_otpr * )) ;
_PROTOTYPE(extern void d_temp, (d_otpr * )) ;
_PROTOTYPE(extern void d_utmp, (d_otpr * )) ;

/*------------------------------*/
/*  routines for #-expressions  */
/*------------------------------*/
_PROTOTYPE(void   c_cadd, (d_otpr *cr,d_otpr *ci,a_cmpx a));
_PROTOTYPE(void   c_csub, (d_otpr *cr,d_otpr *ci,a_cmpx a));
_PROTOTYPE(void   c_dadd, (d_otpr *cr,d_otpr *ci,d_otpc a));
_PROTOTYPE(d_otpc c_dsta, (d_otpr cr,d_otpr ci));
_PROTOTYPE(void   c_dsub, (d_otpr *cr,d_otpr *ci,d_otpc a));
_PROTOTYPE(void   c_padd, (d_otpr *cr,d_otpr *ci,a_cmpx a,a_cmpx b));
_PROTOTYPE(void   c_psub, (d_otpr *cr,d_otpr *ci,a_cmpx a,a_cmpx b));
_PROTOTYPE(void   c_rcad, (d_otpr *cr,d_otpr *ci,a_real a,a_cmpx b));
_PROTOTYPE(void   c_rcsb, (d_otpr *cr,d_otpr *ci,a_real a,a_cmpx b));
_PROTOTYPE(a_cmpx c_stad, (d_otpr cr,d_otpr ci));
_PROTOTYPE(a_cmpx c_stan, (d_otpr cr,d_otpr ci));
_PROTOTYPE(a_cmpx c_stau, (d_otpr cr,d_otpr ci));
_PROTOTYPE(void   d_ass, (d_otpr *a,d_otpr b));
_PROTOTYPE(void   d_assc, (d_otpc *a,d_otpc b));
_PROTOTYPE(void   d_assi, (d_otpi *a,d_otpi b));
_PROTOTYPE(void   d_assz, (d_otpz *a,d_otpz b));
_PROTOTYPE(void   d_clr, (d_otpr *a));
_PROTOTYPE(void   d_dadd, (d_otpr *a,d_otpr b));
_PROTOTYPE(void   d_dsub, (d_otpr *a,d_otpr b));
_PROTOTYPE(a_bool d_eq, (d_otpr a,d_otpr b));
_PROTOTYPE(void   d_free, (d_otpr *a));
_PROTOTYPE(a_bool d_ge, (d_otpr a,d_otpr b));
_PROTOTYPE(a_bool d_gt, (d_otpr a,d_otpr b));
_PROTOTYPE(void   d_init, (d_otpr *a));
_PROTOTYPE(a_bool d_le, (d_otpr a,d_otpr b));
_PROTOTYPE(a_bool d_lt, (d_otpr a,d_otpr b));
_PROTOTYPE(a_bool d_ne, (d_otpr a,d_otpr b));
_PROTOTYPE(void   d_padd, (d_otpr*c,a_real a,a_real b));
_PROTOTYPE(void   d_psub, (d_otpr*c,a_real a,a_real b));
_PROTOTYPE(void   d_radd, (d_otpr*c,a_real a));
_PROTOTYPE(void   d_rsub, (d_otpr*c,a_real a));
_PROTOTYPE(a_real d_stad, (d_otpr a));
_PROTOTYPE(a_real d_stan, (d_otpr a));
_PROTOTYPE(a_real d_stau, (d_otpr a));
_PROTOTYPE(void   i_dadd, (d_otpr*cl,d_otpr *cu,d_otpi a));
_PROTOTYPE(d_otpi i_dsta, (d_otpr cl,d_otpr  cu));
_PROTOTYPE(void   i_dsub, (d_otpr*cl,d_otpr *cu,d_otpi a));
_PROTOTYPE(void   i_iadd, (d_otpr*cl,d_otpr *cu,a_intv a));
_PROTOTYPE(a_intv i_ista, (d_otpr cl,d_otpr  cu));
_PROTOTYPE(void   i_isub, (d_otpr*cl,d_otpr *cu,a_intv a));
_PROTOTYPE(void   i_padd, (d_otpr*cl,d_otpr *cu,a_intv a,a_intv b));
_PROTOTYPE(void   i_psub, (d_otpr*cl,d_otpr *cu,a_intv a,a_intv b));
_PROTOTYPE(void   i_riad, (d_otpr*cl,d_otpr *cu,a_real a,a_intv b));
_PROTOTYPE(void   i_risb, (d_otpr*cl,d_otpr *cu,a_real a,a_intv b));
_PROTOTYPE(void   z_ciad, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  a_cmpx a,a_intv b));
_PROTOTYPE(void   z_cisb, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  a_cmpx a,a_intv b));
_PROTOTYPE(void   z_czad, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  a_cmpx a,a_cinv b));
_PROTOTYPE(void   z_czsb, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  a_cmpx a,a_cinv b));
_PROTOTYPE(void   z_dadd, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  d_otpz a));
_PROTOTYPE(d_otpz z_dsta, (d_otpr crl,d_otpr cil,d_otpr cru,d_otpr ciu));
_PROTOTYPE(void   z_dsub, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  d_otpz a));
_PROTOTYPE(void   z_izad, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  a_intv a,a_cinv b));
_PROTOTYPE(void   z_izsb, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  a_intv a,a_cinv b));
_PROTOTYPE(void   z_padd, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  a_cinv a,a_cinv b));
_PROTOTYPE(void   z_psub, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  a_cinv a,a_cinv b));
_PROTOTYPE(void   z_rzad, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  a_real a,a_cinv b));
_PROTOTYPE(void   z_rzsb, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  a_real a,a_cinv b));
_PROTOTYPE(void   z_zadd, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  a_cinv a));
_PROTOTYPE(a_cinv z_zsta, (d_otpr crl,d_otpr cil,d_otpr cru,d_otpr ciu));
_PROTOTYPE(void   z_zsub, (d_otpr*crl,d_otpr*cil,d_otpr*cru,d_otpr*ciu,
					  a_cinv a));


/*--------------------*/
/*   Runtime Checks   */
/*--------------------*/
/* Index Check */
_PROTOTYPE(extern a_btyp a_ixch, (a_intg index, a_intg lb, a_intg ub)) ;  
						       /* returns index - lb */
/* Pointer Check */
_PROTOTYPE(extern a_VOID a_nilc, (a_VOID pointer)) ;    
							  /* returns pointer */


/*--------------------*/
/*   Dynamic Arrays   */
/*--------------------*/
				       /* Array descriptor */
typedef struct {
	 a_intg lbound, ubound;
	 size_t stride ;
	}   y_bnds ;
#define y_arof( typ, dim) \
   struct                 \
    {  typ   * array ;    \
       a_byte  subarr ;   \
       a_byte  destroy ;  \
       a_byte  numdim ;   \
       size_t  elsize ;   \
       size_t  elnum ;    \
       y_bnds  fd [dim] ; \
    }

#ifdef LINT_ARGS

typedef  y_arof (void, 255)   y_desc ;
/* The typedef may cause an error message */
/* typedef  a_VOID y_dscp ; */
typedef  y_desc * y_dscp ;
     /*       some C compiler report: suspicious pointer conversion */

#else
typedef  y_arof(char, 255)   y_desc ;
/* typedef  a_VOID y_dscp ; */
typedef  y_desc * y_dscp ;
#endif

						/* Stride initialization */
_PROTOTYPE(extern void y_init, (y_dscp d, a_byte dim, size_t elsize)) ;
_PROTOTYPE(extern void y_inid, (y_dscp d, a_byte dim, size_t elsize)) ;
_PROTOTYPE(extern void y_free, (y_dscp d)) ;

#define y_lbnd(d,n)  (y_alck((y_dscp)(d))->fd[(n)-1].lbound)
#define y_ubnd(d,n)  (y_alck((y_dscp)(d))->fd[(n)-1].ubound)

				       /* Array access without indexcheck*/
#define y_inx1( d, i) \
 ((d).array[((i)-(d).fd[0].lbound)*(d).fd[0].stride])
							/* vector */

#define y_inx2( d, i, j)  \
 ((d).array[((i)-(d).fd[0].lbound)*(d).fd[0].stride \
	   +((j)-(d).fd[1].lbound)*(d).fd[1].stride])
							/* matrix */

#define y_inx3( d, i, j, k) \
 ((d).array[((i)-(d).fd[0].lbound)*(d).fd[0].stride \
	   +((j)-(d).fd[1].lbound)*(d).fd[1].stride \
	   +((k)-(d).fd[2].lbound)*(d).fd[2].stride])
							/* tensor */

  /* later, if a dynamic array may be component variable, then */
  /*  side effects in d become dangerous. d has to be a value parameter. */


				       /* Array access with indexcheck*/
#define y_ixc1( d, i) \
 ((d).array[y_ixch((i),y_alck((y_dscp)&(d))->fd[0])])
							/* vector */

#define y_ixc2( d, i, j) \
 ((d).array[y_ixch((i),y_alck((y_dscp)&(d))->fd[0])+y_ixch((j),(d).fd[1])])
							/* matrix */

#define y_ixc3( d, i, j, k) \
 ((d).array[y_ixch((i),y_alck((y_dscp)&(d))->fd[0])+y_ixch((j),(d).fd[1]) \
	   +y_ixch((k),(d).fd[2])])
							/* tensor */


/*   Array access for #-expressions    */

				       /* Array access without indexcheck*/
#define y_ynx1( d, i) \
 ((d).array[(i)*(d).fd[0].stride])
							/* vector */

#define y_ynx2( d, i, j) \
 ((d).array[(i)*(d).fd[0].stride \
	   +(j)*(d).fd[1].stride])
							/* matrix */

#define y_ynx3( d, i, j, k) \
 ((d).array[(i)*(d).fd[0].stride \
	   +(j)*(d).fd[1].stride \
	   +(k)*(d).fd[2].stride])
							/* tensor */

  /* later, if a dynamic array may be component variable, then */
  /*  side effects in d become dangerous. d has to be a value parameter. */


				       /* Array access with indexcheck*/
#define y_yxc1( d, i) \
 ((d).array[y_yxch ((i),y_alck((y_dscp)&(d))->fd[0])])
							/* vector */

#define y_yxc2( d, i, j) \
 ((d).array[y_yxch((i),y_alck((y_dscp)&(d))->fd[0])+y_yxch((j),(d).fd[1])])
							/* matrix */

#define y_yxc3( d, i, j, k) \
 ((d).array[y_yxch((i),y_alck((y_dscp)&(d))->fd[0])+y_yxch((j),(d).fd[1]) \
	   +y_yxch((k),(d).fd[2])])
							/* tensor */

#define y_aled(d) ((a_bool) (((y_dscp)(d))->array != NULL))
						    /* check if dyn. array is allocated */

				       /* Index check */
_PROTOTYPE(extern a_btyp y_ixch, (a_intg index, y_bnds bounds)) ;
_PROTOTYPE(extern a_btyp y_yxch, (a_intg index, y_bnds bounds)) ;

				       /* Dynamic array assignment */
_PROTOTYPE(extern void y_asgn, (y_dscp target, y_dscp source)) ;
_PROTOTYPE(extern void y_vlcp, (y_dscp  d)) ;
_PROTOTYPE(extern void y_temp, (y_dscp  d)) ;
_PROTOTYPE(extern void y_utmp, (y_dscp  d)) ;

_PROTOTYPE(extern y_dscp y_alck, (y_dscp d )) ;

#ifdef LINT_ARGS

#ifdef IBMva
extern a_VOID y_inxn( ) ;
extern a_VOID y_ixcn( ) ;
extern a_VOID y_suba( ) ;
extern a_VOID y_stat( ) ;
extern a_VOID y_ynxn( ) ;
extern a_VOID y_yxcn( ) ;
extern void   y_new ( ) ;
#else
				       /*--------------*/
				       /* Array access */
extern a_VOID y_inxn (y_dscp ,...) ; /*  without Index check */
extern a_VOID y_ixcn (y_dscp ,...) ; /*  with Index check */

				       /* Subarray access: */
				       /* Subarray descriptor generation */
extern a_VOID y_suba (a_VOID md, a_VOID sd,
		      a_char * mode ,... ) ;

				       /* convert static to dynamic array */
extern a_VOID y_stat (y_dscp d, a_VOID statarray, size_t elsize,
		      a_byte dim, ...) ;
extern a_VOID y_ynxn (y_dscp d,...) ; /*  without Index check */
extern a_VOID y_yxcn (y_dscp d,...) ; /*  with Index check */
extern void   y_new  (y_dscp d,...) ; /*  with Index check */
#endif

#else
extern a_VOID y_inxn () ;
extern a_VOID y_ixcn () ;
extern a_VOID y_stat () ;
extern a_VOID y_ynxn () ;
extern a_VOID y_yxcn () ;
extern a_VOID y_suba () ;
extern void   y_new  () ;
#endif


/*-----------------*/
/*  File Handling  */
/*-----------------*/

#define f_fnsz 64
						/* Filedescriptor */

#define f_ilof(typ)  struct {\
 FILE * fp ;\
 unsigned eof  : 1 ;\
 unsigned eoln : 1 ;\
 unsigned text : 1 ;\
 unsigned infl : 1 ;\
 unsigned outf : 1 ;\
 unsigned stdi : 1 ;\
 unsigned stdo : 1 ;\
 unsigned asgd : 1 ;\
 unsigned err  : 1 ;\
 unsigned temp : 1 ;\
 unsigned pp   : 1 ;\
 size_t   ellen ;\
 char   name[f_fnsz] ;\
 char *org;\
 a_VOID next;\
 union {\
   double f; /* forces alignm */\
   typ dow ;\
   char ch [sizeof (typ)] ;\
 } win ;\
}
typedef f_ilof (a_char) f_text ;

extern f_text f_inpu, f_outp ;

_PROTOTYPE(extern void p_init, (int argc, char ** argv) );
_PROTOTYPE(extern void f_assg, (f_text * filevar, char * filename, size_t ellen) );
_PROTOTYPE(extern void f_eofp, (void) );
_PROTOTYPE(extern void f_rset, (f_text*, a_char *, a_char *) );
_PROTOTYPE(extern void f_rwri, (f_text*, a_char *, a_char *) );
_PROTOTYPE(extern void f_get_, (f_text*) );
_PROTOTYPE(extern void f_put_, (f_text*) );
_PROTOTYPE(extern void f_read, (f_text*, a_VOID) );
_PROTOTYPE(extern void f_writ, (f_text*, a_VOID) );
_PROTOTYPE(extern void f_wrc1, (f_text*, a_char *, a_intg) );
_PROTOTYPE(extern void f_wrc2, (f_text*, a_char *, a_intg, a_intg) );
_PROTOTYPE(extern void f_free, (f_text*)) ;


/*-----------------------------------------------*/
/*   standard dynamic array types of Pascal-SC   */
/*-----------------------------------------------*/
typedef y_arof (a_real,1) a_rvty ;
typedef y_arof (a_cmpx,1) a_cvty ;
typedef y_arof (a_intv,1) a_ivty ;
typedef y_arof (a_cinv,1) a_civt ;
typedef y_arof (a_real,2) a_rmty ;
typedef y_arof (a_cmpx,2) a_cmty ;
typedef y_arof (a_intv,2) a_imty ;
typedef y_arof (a_cinv,2) a_cimt ;

/*-----------------------*/
/*  Standard functions:  */
/*-----------------------*/

#define a_chr_     (a_char)
#define a_chrx(x) ((a_char)a_ixch((x),(a_intg)0,(a_intg)255))
#define a_odd_(x)  (1&(x))
#define s_ytos(x)  ((x).ptr)





