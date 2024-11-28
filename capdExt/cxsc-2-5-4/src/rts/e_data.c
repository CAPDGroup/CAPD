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

/* CVS $Id: e_data.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_data.c                              */
/*                                                              */
/*      Description     : Declaration of constant data for      */
/*                        error handling routines.              */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

local aentry e_aset =
    { INDEX_RANGE    , TRUE , TRUE , EXIT_CONT, TRUE , TRUE,
      (aentry *)NULL , e_xset };
local aentry e_aarg =
    { INV_ARG        , TRUE , TRUE , EXIT_CONT, TRUE , TRUE ,
      &e_aset , e_xarg };
local aentry e_aiob =
    { I_O_BUFFER     , TRUE , TRUE , EXIT_CONT, TRUE , TRUE ,
      &e_aarg , e_xiob };
local aentry e_aios =
    { I_O_ERROR      , TRUE , TRUE , EXIT_CONT, TRUE , TRUE ,
      &e_aiob , e_xios };
local aentry e_aall =
    { ALLOCATION     , TRUE , TRUE , EXIT_CONT, TRUE , TRUE ,
      &e_aios , e_xall };
local aentry e_aine =
    { INEXACT        , TRUE , TRUE , TRUE,      TRUE , TRUE ,
      &e_aall , e_xine };
local aentry e_aufl =
    { UNDERFLOW      , TRUE , TRUE , TRUE,      TRUE , TRUE ,
      &e_aine , e_xufl };
local aentry e_aofl =
    { OVERFLOW       , TRUE , TRUE , EXIT_CONT, TRUE , TRUE ,
      &e_aufl , e_xofl };
local aentry e_adbz =
    { DIV_BY_ZERO    , TRUE , TRUE , EXIT_CONT, TRUE , TRUE ,
      &e_aofl , e_xdbz };
local aentry e_aiop =
    { INV_OP         , TRUE , TRUE , EXIT_CONT, TRUE , TRUE ,
      &e_adbz , e_xiop };
local aentry e_anor =
    { NO_ERROR       , TRUE , TRUE , TRUE , TRUE , TRUE ,
      &e_aiop , e_xnor };

#if VAX_VMS_C
#ifdef LINT_ARGS
local void e_data(void)
#else
local void e_data()
#endif
        {
        }
#endif





