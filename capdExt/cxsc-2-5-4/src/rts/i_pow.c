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

/* CVS $Id: i_pow.c,v 1.22 2014/01/30 17:24:09 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : i_pow.c                  b**e         */
/*                                                              */
/*      Entries         : a_intv i_pow(b, e)                    */
/*                        a_intv b, e;                          */
/*                                                              */
/*      Arguments       : b = interval argument base            */
/*                      : e = interval argument exponent        */
/*                                                              */
/*      Description     : Interval power function               */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_one_,*r_two_,*r_zero;
#endif

#ifdef LINT_ARGS
static a_btyp typ1(a_intv b,a_intv e,a_intv *r);
static a_btyp typ2(a_intv b,a_real ei,a_intv *r);
static a_btyp minmax(a_intv b,a_intv e,a_intv *r);
#else
static a_btyp typ1(), typ2(), minmax();
#endif

#ifdef LINT_ARGS
local a_intv i_pow( a_intv b, a_intv e )
#else
local a_intv i_pow( b, e )
a_intv b;
a_intv e;
#endif

{
        a_btyp        rc;
        a_real            ei;
        a_intv          res;

        E_SPUSH("i_pow")

        if (r_sign(e.INF)<0)
            R_ASSIGN(ei,r_umin(r_flor(r_abs(e.INF))));
        else
            R_ASSIGN(ei,r_flor(e.INF));

        if (0==i_iv(b) || 0==i_iv(e)) rc = 1; /* error no interval */

        /* --- 2 types --- */
        else
        if (r_sign(e.INF)>0 && i_point(e) && r_eq(ei,e.INF) &&
            r_sign(b.INF)<=0 && r_sign(b.SUP)>=0)
            rc = typ2(b, ei, &res);
        else if (!(r_sign(b.INF)<=0 && r_sign(b.SUP)>=0))
            rc = typ1(b, e, &res);
        else rc = 1;

        /* --- error --- */
        if (rc)
            e_trap(INV_ARG,8,E_TDBL+E_TEXT(5),&b.INF,
                             E_TDBL+E_TEXT(6),&b.SUP,
                             E_TDBL+E_TEXT(5),&e.INF,
                             E_TDBL+E_TEXT(6),&e.SUP);

        /* --- exit --- */
        E_SPOPP("i_pow")
        return(res);
}

/****************************************************************/

#ifdef LINT_ARGS
static a_btyp typ1(a_intv b,a_intv e,a_intv *r)
#else
static a_btyp typ1(b, e, r)
a_intv  b;  /* base */
a_intv  e;  /* exp  */
a_intv  *r; /* result */
#endif
{
        a_btyp          rc;
        a_real          bl=0, bu=0, el=0, eu=0;
	/* Initialisation to avoid warning 'may be used uninitialized in this function' */
        a_real          dummy;
        int             cmpbas, cmpexp;

        /* --- 1 --- */
        if (r_le(b.SUP,*r_one_)) {
            cmpbas = -1;
            R_ASSIGN(eu,e.INF);
            R_ASSIGN(el,e.SUP);
        }
        else if (r_ge(b.INF,*r_one_)) {
            cmpbas = 1;
            R_ASSIGN(eu,e.SUP);
            R_ASSIGN(el,e.INF);
        }
        else
            cmpbas = 0;

        /* --- 2 --- */
        if (r_sign(e.SUP)<=0) {
            cmpexp = -1;
            R_ASSIGN(bu,b.INF);
            R_ASSIGN(bl,b.SUP);
        }
        else if (r_sign(e.INF)>=0) {
            cmpexp = 1;
            R_ASSIGN(bu,b.SUP);
            R_ASSIGN(bl,b.INF);
        }
        else {
            cmpexp = 0;
            switch(cmpbas) {
            case -1:
               R_ASSIGN(bu,b.INF);
               R_ASSIGN(bl,b.INF);
               break;
            case  0:
               break;
            case  1:
               R_ASSIGN(bu,b.SUP);
               R_ASSIGN(bl,b.SUP);
               break;
            } /* switch */
        }

        /* --- 3 --- */
        if (cmpbas == 0)
            switch(cmpexp) {
            case -1:
               R_ASSIGN(eu,e.INF);
               R_ASSIGN(el,e.INF);
               break;
            case  0:
               break;
            case  1:
               R_ASSIGN(eu,e.SUP);
               R_ASSIGN(el,e.SUP);
               break;
            } /* switch */

        /* --- case double function call --- */
        if (cmpexp == 0 && cmpbas == 0) {
            rc =  minmax(b, e, r);
        }
        /* --- single function call --- */
        else {
            rc =  i_inv2(Lpower,bl,el,&(r->INF),&dummy);
            rc += i_inv2(Lpower,bu,eu,&dummy,&(r->SUP));
        }

        return rc;
} /* typ1() */

/****************************************************************/

#ifdef LINT_ARGS
static a_btyp typ2(a_intv b,a_real ei,a_intv *r)
#else
static a_btyp typ2(b, ei, r)
a_intv  b;
a_real    ei;
a_intv  *r;
#endif
{
        a_btyp          rc;
        a_real          bl, bu, eihalf, eih;
        a_real          dummy;

        R_ASSIGN(eihalf,r_divn(ei,*r_two_));
        if (r_sign(eihalf)<0)
            R_ASSIGN(eih,r_umin(r_flor(r_abs(eihalf))));
        else
            R_ASSIGN(eih,r_flor(eihalf));

        if (r_eq(eihalf,eih)) { /* even(ei) */
            R_ASSIGN(bl,*r_zero);
            R_ASSIGN(bu,(r_gt(r_abs(b.INF),r_abs(b.SUP)) ?
                        r_abs(b.INF) : r_abs(b.SUP)));

        /*  bu = max(r_abs(b.INF),r_abs(b.SUP)); */
        }
        else { /* odd(ei) */
            R_ASSIGN(bl,b.INF);
            R_ASSIGN(bu,b.SUP);
        }

        rc =  i_inv2(Lpower,bl,ei,&(r->INF),&dummy);
        rc += i_inv2(Lpower,bu,ei,&dummy,&(r->SUP));

        return rc;
} /* typ2() */

/****************************************************************/

#ifdef LINT_ARGS
static a_btyp minmax(a_intv b,a_intv e,a_intv *r)
#else
static a_btyp minmax(b, e, r)
a_intv  b;  /* base */
a_intv  e;  /* exp  */
a_intv  *r; /* result */
#endif
{
        a_btyp        rc;
        a_real            r1, r2;
        a_real            dummy;

        rc =  i_inv2(Lpower,b.SUP,e.INF,&r1,&dummy);
        rc += i_inv2(Lpower,b.INF,e.SUP,&r2,&dummy);
        if (r_lt(r1,r2))
           R_ASSIGN(r->INF,r1);
        else
           R_ASSIGN(r->INF,r2);
       /*r->INF = min(r1, r2);*/

        rc += i_inv2(Lpower,b.SUP,e.SUP,&dummy,&r1);
        rc += i_inv2(Lpower,b.INF,e.INF,&dummy,&r2);
        if (r_gt(r1,r2))
           R_ASSIGN(r->SUP,r1);
        else
           R_ASSIGN(r->SUP,r2);
       /*r->SUP = max(r1, r2);*/

        return rc;
} /* minmax() */





