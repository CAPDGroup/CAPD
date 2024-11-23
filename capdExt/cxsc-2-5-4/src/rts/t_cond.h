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

/* CVS $Id: t_cond.h,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

/* -------------------------------------------------------------*/
/* Condition Codes fuer Pruefroutinen                           */
/* -------------------------------------------------------------*/

/* --- in Anlehnung an Koprozessor Condition Codes --- */
           /*0000000000000001B; +Unnormal, not used on the 80387 */
#define PUnorm  0x0001
#define PNAN    0x0002 /*0000000000000010B; +NAN */
           /*0000000000000100B; -Unnormal, not used on the 80387 */
#define MUnorm  0x0004
#define MNAN    0x0008 /*0000000000001000B; -NAN */
#define PNorm   0x0010 /*0000000000010000B; +Normal */
#define PInf    0x0020 /*0000000000100000B; +Infinity */
#define MNorm   0x0040 /*0000000001000000B; -Normal */
#define MInf    0x0080 /*0000000010000000B; -Infinity */
#define PZero   0x0100 /*0000000100000000B; +0 */
#define Empty1  0x0200 /*0000001000000000B; empty */
#define MZero   0x0400 /*0000010000000000B; -0 */
#define Empty2  0x0800 /*0000100000000000B; empty */
#define PDenorm 0x1000 /*0001000000000000B; +Denormal */
/*                       0010000000000000B; empty, not used on the 80387 */
#define MDenorm 0x4000 /*0100000000000000B; -Denormal */
/*                       1000000000000000B; empty, not used on the 80387 */

/* --- DefinitionsBereiche (DB) fuer Funktionen --- */
#define DBPosExtReal   (PNorm)
#define DBNegExtReal   (MNorm)
#define DBPos0ExtReal  (DBPosExtReal|PZero|MZero)
#define DBNeg0ExtReal  (DBNegExtReal|MZero|PZero)

#define DBPosDExtReal  (PNorm|PDenorm)
#define DBNegDExtReal  (MNorm|MDenorm)
#define DBPosD0ExtReal (DBPosDExtReal|PZero)
#define DBNegD0ExtReal (DBNegDExtReal|MZero)

#define DBExtReal      (DBPos0ExtReal|DBNeg0ExtReal)
#define DBInfExtReal   (DBExtReal|PInf|MInf)





