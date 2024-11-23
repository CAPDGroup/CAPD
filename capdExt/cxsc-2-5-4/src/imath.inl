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

/* CVS $Id: imath.inl,v 1.8 2014/01/30 17:23:45 cxsc Exp $ */

#define CXSC_INCLUDE
#undef LINT_ARGS
#include <fi_lib.hpp>
#undef CXSC_INCLUDE

namespace cxsc{
using namespace fi_lib;
inline interval sqrt   (const interval &a)         { return j_sqrt(a); }

inline interval sin    (const interval &a) throw() { return j_sin(a); }
inline interval cos    (const interval &a) throw() { return j_cos(a); }
inline interval tan    (const interval &a) throw() { return j_tan(a); }
inline interval cot    (const interval &a) throw() { return j_cot(a); }

inline interval asin   (const interval &a)         { return j_asin(a); }
inline interval acos   (const interval &a)         { return j_acos(a); }
inline interval atan   (const interval &a)         { return j_atan(a); }
inline interval acot   (const interval &a)         { return j_acot(a); }

inline interval exp    (const interval &a) throw() { return j_exp(a); }
inline interval ln     (const interval &a)         { return j_log(a); }
inline interval log2   (const interval &a)         { return j_log2(a); }
inline interval log10  (const interval &a)         { return j_lg10(a); }

inline interval sinh   (const interval &a) throw() { return j_sinh(a); }
inline interval cosh   (const interval &a) throw() { return j_cosh(a); }
inline interval tanh   (const interval &a) throw() { return j_tanh(a); }
inline interval coth   (const interval &a) throw() { return j_coth(a); }

inline interval asinh  (const interval &a)         { return j_asnh(a); }
inline interval acosh  (const interval &a)         { return j_acsh(a); }
inline interval atanh  (const interval &a)         { return j_atnh(a); }
inline interval acoth  (const interval &a)         { return j_acth(a); }

} // namespace cxsc
