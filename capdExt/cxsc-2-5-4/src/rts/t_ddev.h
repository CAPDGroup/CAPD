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

/* CVS $Id: t_ddev.h,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

#ifdef LINT_ARGS

#ifdef ANSI_C

extern int dadd(const DReal *arg1, const DReal *arg2, DReal *res);
extern int dsub(const DReal *arg1, const DReal *arg2, DReal *res);
extern int msplitee(const ExtReal *arg, ExtReal *ures, ExtReal *lres);
extern int dmshift(const DExp e, const DMant *arg, DMant *res);
extern int dmadjust(const DMant *m, const int len, DMant *mr, DExp *shiftr);

#else  /* NOT ANSI_C */

extern int dadd(DReal *arg1, DReal *arg2, DReal *res);
extern int dsub(DReal *arg1, DReal *arg2, DReal *res);
extern int msplitee(ExtReal *arg, ExtReal *ures, ExtReal *lres);
extern int dmshift(DExp e, DMant *arg, DMant *res);
extern int dmadjust(DMant *m, int len, DMant *mr, DExp *shiftr);

#endif /* ANSI_C */

#else  /* NOT LINT_ARGS */

extern int dadd();
extern int dsub();
extern int msplitee(), dmshift(), dmadjust();

#endif /* LINT_ARGS */





