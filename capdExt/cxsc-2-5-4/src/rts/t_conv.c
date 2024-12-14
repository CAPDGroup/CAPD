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

/* CVS $Id: t_conv.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */


#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/*--------------------------------------------------------------*
 | conv                                                         |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int extreal_to_longreal (const ExtReal *arg, LongReal *res)
#else
int extreal_to_longreal (arg, res)
const ExtReal   *arg;
      LongReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int extreal_to_longreal (ExtReal *arg, LongReal *res)
#else
int extreal_to_longreal (arg, res)
ExtReal   *arg;
LongReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal  r;
   int      ret;

   /* --- ArgPruefung --- */
   if(arg_check==On) ret=chk_arg_1(EToL, arg, &r);
   else { copyee(arg, &r); ret = NoErr; }

   /* --- falls ArgPruefung positiv: etol --- */
   if(NoErr==ret || ExcPInf == ret || ExcMInf == ret)
      (void)_s_etol(arg, res);

   /* --- falls Fehler, versuch's nochmal --- */
   if(NoErr!=ret)
      if (ret!=ExcPInf && ret!=ExcMInf)
/*    if(NoErr==(ret = exc_handle_1(EToL, ret, arg, &r))) */
         (void)_s_etol(&r, res);

   return ret;

} /* extreal_to_longreal() */

/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int longreal_to_extreal (const LongReal *arg, ExtReal *res)
#else
int longreal_to_extreal (arg, res)
const LongReal *arg;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int longreal_to_extreal (LongReal *arg, ExtReal *res)
#else
int longreal_to_extreal (arg, res)
LongReal *arg;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int      ret;

   if (NoErr!=(ret = _s_ltoe(arg, res)))

      /* --- FehlerHandle --- */
      if(NoErr==(ret = exc_handle_1(LToE, ret, res, res)))
                ;       /* eigentlich kein Fehler moeglich */

   return ret;

} /* longreal_to_extreal() */

/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int extreal_to_int (const ExtReal *arg, int *res)
#else
int extreal_to_int (arg, res)
const ExtReal  *arg;
      int      *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int extreal_to_int (ExtReal *arg, int *res)
#else
int extreal_to_int (arg, res)
ExtReal  *arg;
int      *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal  r;
   int      ret;

   /* --- ArgPruefung --- */
   if(arg_check==On) ret=chk_arg_1(EToI, arg, &r);
   else ret = NoErr;

   /* --- falls ArgPruefung positiv: etoi --- */
   if(NoErr==ret) ret = _s_etoi(arg, res);

   /* --- falls Fehler, versuch's nochmal --- */
   if(NoErr!=ret)
      if(NoErr==(ret = exc_handle_1(EToI, ret, arg, &r)))
         ret=_s_etoi(&r, res);

   return ret;
} /* extreal_to_int() */

/* ----------- REST OF FILE IS UNUSED ----------- */

/*
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif
*/
/*--------------------------------------------------------------*
 | iconv                                                        |
 *--------------------------------------------------------------*/
/*
#ifdef ANSI_C
#ifdef LINT_ARGS
int ilongreal_to_extreal (const ILongReal *arg, IExtReal *res)
#else
int ilongreal_to_extreal (arg, res)
const ILongReal  *arg;
      IExtReal   *res;
#endif
#else
#ifdef LINT_ARGS
int ilongreal_to_extreal (ILongReal *arg, IExtReal *res)
#else
int ilongreal_to_extreal (arg, res)
ILongReal  *arg;
IExtReal   *res;
#endif
#endif
{
   int   ret;

   ret = longreal_to_extreal(&(arg->u), &(res->u));
   ret = longreal_to_extreal(&(arg->l), &(res->l));

   return ret;
}
*/
/* ilongreal_to_extreal() */

/*--------------------------------------------------------------*/
/*
#ifdef ANSI_C
#ifdef LINT_ARGS
int iextreal_to_longreal (const IExtReal *arg, ILongReal *res)
#else
int iextreal_to_longreal (arg, res)
const IExtReal   *arg;
      ILongReal  *res;
#endif
#else
#ifdef LINT_ARGS
int iextreal_to_longreal (IExtReal *arg, ILongReal *res)
#else
int iextreal_to_longreal (arg, res)
IExtReal   *arg;
ILongReal  *res;
#endif
#endif
{
   int   rnd;
   int   retu, retl;

   rnd = getrndmode();
   setrndmode(UP);
   retu = extreal_to_longreal(&(arg->u), &(res->u));
   setrndmode(DOWN);
   retl = extreal_to_longreal(&(arg->l), &(res->l));
   setrndmode(rnd);

   return max(retu, retl);
}
*/
/* iextreal_to_longreal() */






