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

/* CVS $Id: rtsrmath.h,v 1.22 2014/01/30 17:23:48 cxsc Exp $ */

#ifndef _CXSC_RTSRMATH_H_INCLUDED
#define _CXSC_RTSRMATH_H_INCLUDED


#ifdef __cplusplus
extern "C" {
#endif

// Ten-Byte Arithmetic
#define  ExtReal        TByte
#define  LongReal       double

#define  DOWN (-UP)
#define  NEAR 0
#define  UP   1
#define  CHOP 2

extern int t_ltoe(const LongReal *arg, ExtReal *res); // Convert LongReal to Extreal
extern int t_etol(const ExtReal *arg, LongReal *res); // Convert Extreal to Longreal

extern int t_grnd(void); // Get Round Mode
extern int t_srnd(int); // Set Round Mode

extern int t_sqte(const  ExtReal *arg,  ExtReal *res); // Sqrt
extern int t_sqme(const  ExtReal *arg,  ExtReal *res); // Sqrtm1

extern int t_sine(const  ExtReal *arg,  ExtReal *res); // Sin
extern int t_cose(const  ExtReal *arg,  ExtReal *res); // Cos
extern int t_tane(const  ExtReal *arg,  ExtReal *res); // Tan
extern int t_cote(const  ExtReal *arg,  ExtReal *res); // Cot

extern int t_asne(const  ExtReal *arg,  ExtReal *res); // ASin
extern int t_acse(const  ExtReal *arg,  ExtReal *res); // ACos
extern int t_atne(const  ExtReal *arg,  ExtReal *res); // ATan
extern int t_acte(const  ExtReal *arg,  ExtReal *res); // ACot

extern int t_exme(const  ExtReal *arg,  ExtReal *res); // Expm1 
extern int t_lnpe(const  ExtReal *arg,  ExtReal *res); // Lnp1

extern int t_expe(const  ExtReal *arg,  ExtReal *res); // Exp 
extern int t_lnee(const  ExtReal *arg,  ExtReal *res); // Ln

extern int t_snhe(const  ExtReal *arg,  ExtReal *res); // SinH
extern int t_cshe(const  ExtReal *arg,  ExtReal *res); // CosH
extern int t_tnhe(const  ExtReal *arg,  ExtReal *res); // TanH
extern int t_cthe(const  ExtReal *arg,  ExtReal *res); // CotH

extern int t_ashe(const  ExtReal *arg,  ExtReal *res); // ASinH
extern int t_ache(const  ExtReal *arg,  ExtReal *res); // ACosH
extern int t_anhe(const  ExtReal *arg,  ExtReal *res); // ATanH
extern int t_athe(const  ExtReal *arg,  ExtReal *res); // ACotH

extern int t_powe(const  ExtReal *bas,  const ExtReal *exp,  ExtReal *res); // Pow

#ifdef __cplusplus
}
#endif

#endif // _CXSC_RTSRMATH_H_INCLUDED
