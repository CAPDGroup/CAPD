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

/* CVS $Id: o_fcth.h,v 1.21 2014/01/30 17:24:10 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : o_fcth.h                              */
/*                                                              */
/*      Description     : function prototypes                   */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/int/a_fcth.h"
#include "/u/p88c/runtime/base/b_fcth.h"
#include "/u/p88c/runtime/dot/d_fcth.h"
#include "/u/p88c/runtime/dyn/y_fcth.h"
#include "/u/p88c/runtime/error/e_fcth.h"
#include "/u/p88c/runtime/file/f_fcth.h"
#include "/u/p88c/runtime/multi/l_fcth.h"
#include "/u/p88c/runtime/real/r_fcth.h"
#include "/u/p88c/runtime/set/s_fcth.h"
#ifndef DEC_ARITH
#include "/u/p88c/runtime/tbyte/t_fcth.h"
#endif
#else
#include "a_fcth.h"
#include "b_fcth.h"
#include "d_fcth.h"
#include "e_fcth.h"
#include "f_fcth.h"
#include "l_fcth.h"
#include "r_fcth.h"
#include "s_fcth.h"
#ifndef DEC_ARITH
#include "t_fcth.h"
#endif
#include "y_fcth.h"
#endif





