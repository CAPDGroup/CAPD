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

/* CVS $Id: t_exc_.h,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/* ------------------------------------------------------------- */
/* Exception Return Codes der IEEE Math Funktionen               */
/* ------------------------------------------------------------- */

#define NoErr           0       /* no error                          */

/* --- in math.h definierte:  */
/* math.h aus ieee.h entfernt, da Namenskonflikte mit PXSC-Runtime */
/* bei Bezeichnern OVERFLOW, UNDERFLOW und (hier nicht verwendet) MAXINT. */
#define DOMAIN          1       /* argument domain error             */
#define SING            2       /* argument singularity              */
#define OVER_FLOW       3       /* overflow range error              */
#define UNDER_FLOW      4       /* underflow range error             */
#define TLOSS           5       /* total loss of precision           */
#define PLOSS           6       /* partial loss of precision         */
/* --- ende math.h --- */

/* --- exception codes (abgeleitet aus Koprozessor Exceptions) --- */
#define ExcPUnorm       101     /* +Unormal                          */
#define ExcPNAN         102     /* +NAN                              */
#define ExcMUnorm       103     /* -Unormal                          */
#define ExcMNAN         104     /* -NAN                              */
#define ExcPNorm        105     /* +Normal                           */
#define ExcPInf         106     /* +Infinity                         */
#define ExcMNorm        107     /* -Normal                           */
#define ExcMInf         108     /* -Infinity                         */
#define ExcPZero        109     /* +0                                */
#define ExcMZero        111     /* -0                                */
#define ExcEmpty        112     /* empty                             */
#define ExcPDenorm      113     /* +Denormal                         */
#define ExcMDenorm      115     /* -Denormal                         */

/* --- weitere --- */
#define ExcInvalid      220     /* ungueltiges Argument              */
#define ExcNoI          250     /* Intervallgrenzen vertauscht       */
#define ExcISing        260     /* Intervallargument enthaelt Pol    */
#define ExcUnknown      999     /* falls unbestimmt                  */

/* 241090 Baumhof */
#define ExcDBZ          270     /* Division durch Null               */
#define ExcDIZ          280     /* Division durch Intervall mit 0    */





