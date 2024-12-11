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

/* CVS $Id: xscclass.hpp,v 1.33 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_XSCCLASS_HPP_INCLUDED
#define _CXSC_XSCCLASS_HPP_INCLUDED

#include "compiler.h"

// #ifndef CXSC_INDEX_CHECK // 4.10.00
// #  define CXSC_INDEX_CHECK 1
// #endif
#ifndef CXSC_INDEX_CHECK
  #define CXSC_INDEX_CHECK 0
#endif

#ifdef OLD_CXSC

//#define index CXSC_index


namespace cxsc {

class index   // for backwards compatibility
{
private:
  int ind;
public:
  index() {}
  index(const int i) { ind=i; }
  int _int() const { return ind; }
};

} // namespace cxsc 

//typedef class _index index;

#endif

#define _CXSC_INLINE
#ifdef _CXSC_INLINE
#  define INLINE inline
#  define TINLINE inline   // templates in vector.inl, matrix.inl
#  define _CXSC_INCL_INL
#  define _CXSC_FRIEND_TPL // define friends in *vector.hpp, *matrix.hpp
#else
#  define INLINE
#  define TINLINE
#endif

#undef _CXSC_CPP
#ifdef _CXSC_CPP
#  undef INLINE
#  define INLINE
#  undef TINLINE
#  define TINLINE inline
#  define _CXSC_FRIEND_TPL
#  undef _CXSC_INCL_INL
#endif


// wegen gcc 2.95.1 <----------------------------------------================
// #undef TINLINE
// #define TINLINE

#include <string>

namespace cxsc {

using std::string; // this does the job for all files it is included

class real;             //     real.hpp
class interval;         // interval.hpp
class complex;          //  complex.hpp
class cinterval;        // cinterva.hpp   complex interval

class l_real;           //   l_real.hpp   staggered real
class l_interval;       // l_interv.hpp   staggered interval
class l_complex;
class l_cinterval;

class lx_interval;   // lx_interval.hpp   extended staggered interval
class lx_real;       // lx_real.hpp       extended staggered real
class lx_cinterval;  // lx_cinterval.hpp  extended staggered complex interval
class lx_complex;    // lx_complex.hpp    extended staggered complex

class dotprecision;     //      dot.hpp   
class idotprecision;    //     idot.hpp   interval dotprecision
class cdotprecision;    //     cdot.hpp   complex dotprecision
class cidotprecision;   //    cidot.hpp   complex interval dotprecision

class intvector;
class intvector_slice;
class rvector;
class rvector_slice;
class ivector;
class ivector_slice;
class cvector;
class cvector_slice;
class civector;
class civector_slice;
class l_rvector;
class l_rvector_slice;
class l_ivector;
class l_ivector_slice;

class intmatrix;
class intmatrix_slice;
class intmatrix_subv;
class rmatrix;
class rmatrix_slice;
class rmatrix_subv;
class imatrix;
class imatrix_slice;
class imatrix_subv;
class cmatrix;
class cmatrix_slice;
class cmatrix_subv;
class cimatrix;
class cimatrix_slice;
class cimatrix_subv;
class l_rmatrix;
class l_rmatrix_slice;
class l_rmatrix_subv;
class l_imatrix;
class l_imatrix_slice;
class l_imatrix_subv;

class srvector;
class sivector;
class scvector;
class scivector;
class srvector_slice;
class sivector_slice;
class scvector_slice;
class scivector_slice;

class srmatrix;
class simatrix;
class scmatrix;
class scimatrix;
class srmatrix_slice;
class simatrix_slice;
class scmatrix_slice;
class scimatrix_slice;
class srmatrix_subv;
class simatrix_subv;
class scmatrix_subv;
class scimatrix_subv;

inline string nameof(bool)    { return "bool"; }
inline string nameof(char)    { return "char"; }
inline string nameof(int)     { return "int"; }
inline string nameof(long)    { return "long"; }

inline string nameof(const float &)   { return "float"; }
inline string nameof(const double &)  { return "double"; }

inline string nameof(const real &)     { return "real"; }
inline string nameof(const interval &) { return "interval"; }
inline string nameof(const complex &)  { return "complex"; }
inline string nameof(const cinterval &){ return "cinterval"; }

inline string nameof(const l_real &)      { return "l_real"; }
inline string nameof(const l_interval &)  { return "l_interval"; }
inline string nameof(const l_complex &)   { return "l_complex"; }
inline string nameof(const l_cinterval &) { return "l_cinterval"; }

inline string nameof(const lx_real &)     { return "lx_real"; }
inline string nameof(const lx_interval &) { return "lx_interval"; }
inline string nameof(const lx_cinterval &){ return "lx_cinterval"; }

inline string nameof(const dotprecision &)   { return "dotprecision"; }
inline string nameof(const idotprecision &)  { return "idotprecision"; }
inline string nameof(const cdotprecision &)  { return "cdotprecision"; }
inline string nameof(const cidotprecision &) { return "cidotprecision"; }

inline string nameof(const intvector &) { return "intvector"; }
inline string nameof(const intvector_slice &) { return "intvector_slice"; }
inline string nameof(const rvector &) { return "rvector"; }
inline string nameof(const rvector_slice &) { return "rvector_slice"; }
inline string nameof(const ivector &) { return "ivector"; }
inline string nameof(const ivector_slice &) { return "ivector_slice"; }
inline string nameof(const cvector &) { return "cvector"; }
inline string nameof(const cvector_slice &) { return "cvector_slice"; }
inline string nameof(const civector &) { return "civector"; }
inline string nameof(const civector_slice &) { return "civector_slice"; }
inline string nameof(const l_rvector &) { return "l_rvector"; }
inline string nameof(const l_rvector_slice &) { return "l_rvector_slice"; }
inline string nameof(const l_ivector &) { return "l_ivector"; }
inline string nameof(const l_ivector_slice &) { return "l_ivector_slice"; }

inline string nameof(const intmatrix &) { return "intmatrix"; }
inline string nameof(const intmatrix_slice &) { return "intmatrix_slice"; }
inline string nameof(const rmatrix &) { return "rmatrix"; }
inline string nameof(const rmatrix_slice &) { return "rmatrix_slice"; }
inline string nameof(const rmatrix_subv &) { return "rmatrix_subv"; }
inline string nameof(const imatrix &) { return "imatrix"; }
inline string nameof(const imatrix_slice &) { return "imatrix_slice"; }
inline string nameof(const imatrix_subv &) { return "imatrix_subv"; }
inline string nameof(const cmatrix &) { return "cmatrix"; }
inline string nameof(const cmatrix_slice &) { return "cmatrix_slice"; }
inline string nameof(const cmatrix_subv &) { return "cmatrix_subv"; }
inline string nameof(const cimatrix &) { return "cimatrix"; }
inline string nameof(const cimatrix_slice &) { return "cimatrix_slice"; }
inline string nameof(const cimatrix_subv &) { return "cimatrix_subv"; }
inline string nameof(const l_rmatrix &) { return "l_rmatrix"; }
inline string nameof(const l_rmatrix_slice &) { return "l_rmatrix_slice"; }
inline string nameof(const l_rmatrix_subv &) { return "l_rmatrix_subv"; }
inline string nameof(const l_imatrix &) { return "l_imatrix"; }
inline string nameof(const l_imatrix_slice &) { return "l_imatrix_slice"; }
inline string nameof(const l_imatrix_subv &) { return "l_imatrix_subv"; }

inline string nameof(const srvector& ) { return "srvector"; }
inline string nameof(const sivector& ) { return "sivector"; }
inline string nameof(const scvector& ) { return "scvector"; }
inline string nameof(const scivector& ) { return "scivector"; }
inline string nameof(const srvector_slice& ) { return "srvector_slice"; }
inline string nameof(const sivector_slice& ) { return "sivector_slice"; }
inline string nameof(const scvector_slice& ) { return "scvector_slice"; }
inline string nameof(const scivector_slice& ) { return "scivector_slice"; }

inline string nameof(const srmatrix& ) { return "srmatrix"; }
inline string nameof(const simatrix& ) { return "simatrix"; }
inline string nameof(const scmatrix& ) { return "scmatrix"; }
inline string nameof(const scimatrix& ) { return "scimatrix"; }
inline string nameof(const srmatrix_slice& ) { return "srmatrix_slice"; }
inline string nameof(const simatrix_slice& ) { return "simatrix_slice"; }
inline string nameof(const scmatrix_slice& ) { return "scmatrix_slice"; }
inline string nameof(const scimatrix_slice& ) { return "scimatrix_slice"; }
inline string nameof(const srmatrix_subv& ) { return "srmatrix_subv"; }
inline string nameof(const simatrix_subv& ) { return "simatrix_subv"; }
inline string nameof(const scmatrix_subv& ) { return "scmatrix_subv"; }
inline string nameof(const scimatrix_subv& ) { return "scimatrix_subv"; }
} // namespace cxsc 

#endif // _CXSC_XSCCLASS_HPP_INCLUDED
