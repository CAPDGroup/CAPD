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

/* CVS $Id: t_drea.h,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

#define DMantLen 24

typedef struct DMant{Digit digit[DMantLen+1];} DMant;
/* +1, da hoechstes digit immer 0 fuer carry */

typedef int DExp;
typedef int DSgn;
typedef struct DReal{
   DMant m;
   DExp  e;
   DSgn  s;
}  DReal;

/* --- Intervalle --- */
typedef struct IDReal{
   DReal u;
   DReal l;
} IDReal;


#ifdef LINT_ARGS

#ifdef ANSI_C

extern int initd(DReal *d);
extern int normd(const DReal *arg, DReal *res);

extern int cmpdd(DReal *arg, DReal *res);
extern int cmpabsdd(/*const*/ DReal *arg, /*const*/ DReal *res);
extern int copydd(const DReal *arg, DReal *res);

extern int adddd(const DReal *arg1, const DReal *arg2, DReal *res);
extern int subdd(const DReal *arg1, const DReal *arg2, DReal *res);
extern int muldd(const DReal *arg1, const DReal *arg2, DReal *res);
extern int muled(const ExtReal *arg1, const ExtReal *arg2, DReal *res);

extern int mulLdE(const ExtReal *arg1, DReal *res);
extern int mulendpi(const ExtReal *arg, int two_or_four, DReal *res);

extern int dreal_to_extreal(const DReal *arg, ExtReal *res);
extern int dreal_to_2extreal(const DReal *arg, ExtReal *resu, ExtReal *resl);
extern int extreal_to_dreal(const ExtReal *arg, DReal *res);

extern const DReal DOne;

#else  /* NOT ANSI_C */

extern int initd(DReal *d);
extern int normd(DReal *arg, DReal *res);

extern int cmpdd(DReal *arg, DReal *res);
extern int cmpabsdd(DReal *arg, DReal *res);
extern int copydd(DReal *arg, DReal *res);

extern int adddd(DReal *arg1, DReal *arg2, DReal *res);
extern int subdd(DReal *arg1, DReal *arg2, DReal *res);
extern int muldd(DReal *arg1, DReal *arg2, DReal *res);
extern int muled(ExtReal *arg1, ExtReal *arg2, DReal *res);

extern int mulLdE(ExtReal *arg1, DReal *res);
extern int mulendpi(ExtReal *arg, int two_or_four, DReal *res);

extern int dreal_to_extreal(DReal *arg, ExtReal *res);
extern int dreal_to_2extreal(DReal *arg, ExtReal *resu, ExtReal *resl);
extern int extreal_to_dreal(ExtReal *arg, DReal *res);

extern DReal DOne;

#endif /* ANSI_C */

#else  /* NOT LINT_ARGS */

extern int initd(), normd();
extern int cmpdd(), cmpabsdd(), copydd();
extern int adddd(), subdd(), muldd(), muled();
extern int mulLdE(), mulendpi();

extern int dreal_to_extreal(), dreal_to_2extreal(), extreal_to_dreal();

extern DReal DOne;

#endif /* LINT_ARGS */





