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

/* CVS $Id: sparsematrix.hpp,v 1.11 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_SPARSEMATRIX_HPP_INCLUDED
#define _CXSC_SPARSEMATRIX_HPP_INCLUDED

namespace cxsc {

//Helper class for construction of sparse matrices

template<class T>
class triplet_store {
  public:
    int col;
    int row;
    T val;
 
    triplet_store(int i, int j, T v) : col(j), row(i),  val(v) { }

    bool operator<(const triplet_store& t) const {
      if(col < t.col) {
        return true;
      } else if(col > t.col) {
        return false;
      } else {
        if(row < t.row)
          return true;
        else
          return false;
      }
    }
};

#if(CXSC_INDEX_CHECK)
template<class TA, class Tx, class Tres, class TDot, class TElement>
inline Tres spsl_mv_mult(const TA&, const Tx&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class Tx, class Tres, class TDot, class TElement>
inline Tres spsl_mv_mult(const TA&, const Tx&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class Tx, class Tres, class TDot, class TElement>
inline Tres spsp_mv_mult(const TA&, const Tx&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class Tx, class Tres, class TDot, class TElement>
inline Tres spsp_mv_mult(const TA&, const Tx&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class Tx, class Tres, class TDot>
inline Tres spf_mv_mult(const TA&, const Tx&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class Tx, class Tres, class TDot>
inline Tres spf_mv_mult(const TA&, const Tx&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class Tx, class Tres, class TDot>
inline Tres fsp_mv_mult(const TA&, const Tx&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class Tx, class Tres, class TDot>
inline Tres fsp_mv_mult(const TA&, const Tx&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class Tx, class Tres, class TDot>
inline Tres fsl_mv_mult(const TA&, const Tx&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class Tx, class Tres, class TDot>
inline Tres fsl_mv_mult(const TA&, const Tx&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres, class TDot, class TElement>
inline Tres spsp_mm_mult(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres, class TDot, class TElement>
inline Tres spsp_mm_mult(const TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres, class TDot>
inline Tres fsp_mm_mult(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres, class TDot>
inline Tres fsp_mm_mult(const TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres, class TDot>
inline Tres spf_mm_mult(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres, class TDot>
inline Tres spf_mm_mult(const TA&, const TB&) throw();
#endif

template<class TA, class Ts, class Tres>
inline Tres sp_ms_div(const TA&, const Ts&);

template<class TA, class Ts, class Tres>
inline Tres sp_ms_mult(const TA&, const Ts&);

template<class Ts, class TA, class Tres>
inline Tres sp_sm_mult(const Ts&, const TA&);

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres, class TElement>
inline Tres spsp_mm_add(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres, class TElement>
inline Tres spsp_mm_add(const TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres>
inline Tres spf_mm_add(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres>
inline Tres spf_mm_add(const TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres>
inline Tres fsp_mm_add(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres>
inline Tres fsp_mm_add(const TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres, class TElement>
inline Tres spsp_mm_sub(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres, class TElement>
inline Tres spsp_mm_sub(const TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres>
inline Tres spf_mm_sub(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres>
inline Tres spf_mm_sub(const TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres>
inline Tres fsp_mm_sub(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres>
inline Tres fsp_mm_sub(const TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres, class TElement>
inline Tres spsp_mm_hull(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres, class TElement>
inline Tres spsp_mm_hull(const TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres>
inline Tres spf_mm_hull(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres>
inline Tres spf_mm_hull(const TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres>
inline Tres fsp_mm_hull(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres>
inline Tres fsp_mm_hull(const TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres, class TElement>
inline Tres spsp_mm_intersect(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres, class TElement>
inline Tres spsp_mm_intersect(const TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres>
inline Tres spf_mm_intersect(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres>
inline Tres spf_mm_intersect(const TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class Tres>
inline Tres fsp_mm_intersect(const TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class Tres>
inline Tres fsp_mm_intersect(const TA&, const TB&) throw();
#endif

template<class TA, class TB>
inline bool spsp_mm_comp(const TA&, const TB&);

template<class TA, class TB>
inline bool spf_mm_comp(const TA&, const TB&);

template<class TA, class TB>
inline bool fsp_mm_comp(const TA&, const TB&);

template<class TA, class TB, class TType>
inline bool spsp_mm_less(const TA&, const TB&);

template<class TA, class TB, class TType>
inline bool spf_mm_less(const TA&, const TB&);

template<class TA, class TB, class TType>
inline bool fsp_mm_less(const TA&, const TB&);

template<class TA, class TB, class TType>
inline bool spsp_mm_leq(const TA&, const TB&);

template<class TA, class TB, class TType>
inline bool spf_mm_leq(const TA&, const TB&);

template<class TA, class TB, class TType>
inline bool fsp_mm_leq(const TA&, const TB&);

template<class TA, class TB, class TType>
inline bool spsp_mm_greater(const TA&, const TB&);

template<class TA, class TB, class TType>
inline bool spf_mm_greater(const TA&, const TB&);

template<class TA, class TB, class TType>
inline bool fsp_mm_greater(const TA&, const TB&);

template<class TA, class TB, class TType>
inline bool spsp_mm_geq(const TA&, const TB&);

template<class TA, class TB, class TType>
inline bool spf_mm_geq(const TA&, const TB&);

template<class TA, class TB, class TType>
inline bool fsp_mm_geq(const TA&, const TB&);

template<class TA, class Tres>
inline Tres sp_m_negative(const TA&);

template<class TA, class TType>
inline std::ostream& sp_m_output(std::ostream&, const TA&);

template<class TA, class TType>
inline std::istream& sp_m_input(std::istream&, TA&);

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class TElement>
inline TA& slsp_mm_assign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class TElement>
inline TA& slsp_mm_assign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class TElement, class TType>
inline TA& slf_mm_assign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class TElement, class TType>
inline TA& slf_mm_assign(TA&, const TB&) throw();
#endif

template<class TA, class TB, class TType>
inline TA& spf_mm_assign(TA&, const TB&);

template<class TA, class Ts>
inline TA& sp_ms_divassign(TA&, const Ts&);

template<class TA, class Ts>
inline TA& sp_ms_multassign(TA&, const Ts&);

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class TDot, class TElement>
inline TA& spsp_mm_multassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class TDot, class TElement>
inline TA& spsp_mm_multassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class TDot, class TFull>
inline TA& spf_mm_multassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class TDot, class TFull>
inline TA& spf_mm_multassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class TDot, class TFull>
inline TA& fsp_mm_multassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class TDot, class TFull>
inline TA& fsp_mm_multassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB>
inline TA& fsp_mm_addassign(TA& A, const TB& B) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB>
inline TA& fsp_mm_addassign(TA& A, const TB& B) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class TFull>
inline TA& spf_mm_addassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class TFull>
inline TA& spf_mm_addassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class TElement>
inline TA& spsp_mm_addassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class TElement>
inline TA& spsp_mm_addassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB>
inline TA& spsp_mm_addassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB>
inline TA& spsp_mm_addassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB>
inline TA& fsp_mm_subassign(TA& A, const TB& B) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB>
inline TA& fsp_mm_subassign(TA& A, const TB& B) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class TFull>
inline TA& spf_mm_subassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class TFull>
inline TA& spf_mm_subassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class TElement>
inline TA& spsp_mm_subassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class TElement>
inline TA& spsp_mm_subassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB>
inline TA& spsp_mm_subassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB>
inline TA& spsp_mm_subassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class TFull>
inline TA& spf_mm_hullassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class TFull>
inline TA& spf_mm_hullassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB>
inline TA& fsp_mm_hullassign(TA& A, const TB& B)  throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB>
inline TA& fsp_mm_hullassign(TA& A, const TB& B)  throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class TElement>
inline TA& spsp_mm_hullassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class TElement>
inline TA& spsp_mm_hullassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB>
inline TA& spsp_mm_hullassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB>
inline TA& spsp_mm_hullassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB>
inline TA& fsp_mm_intersectassign(TA& A, const TB& B)  throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB>
inline TA& fsp_mm_intersectassign(TA& A, const TB& B)  throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class TFull>
inline TA& spf_mm_intersectassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class TFull>
inline TA& spf_mm_intersectassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB, class TElement>
inline TA& spsp_mm_intersectassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB, class TElement>
inline TA& spsp_mm_intersectassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class TB>
inline TA& spsp_mm_intersectassign(TA&, const TB&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class TB>
inline TA& spsp_mm_intersectassign(TA&, const TB&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
inline Tx& svsp_vv_assign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
inline Tx& svsp_vv_assign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
inline Tx& svsl_vv_assign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
inline Tx& svsl_vv_assign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
inline Tx& svf_vv_assign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
inline Tx& svf_vv_assign(Tx&, const Ty&) throw();
#endif

template<class TA, class Ts, class TType>
inline TA& sp_ms_assign(TA&, const Ts&);

template<class TA, class Ts, class TElement, class TType>
inline TA& sl_ms_assign(TA&, const Ts&);

template<class Tx, class Ts>
inline Tx& sv_vs_assign(Tx&, const Ts&);

template<class TA>
inline bool sp_m_not(const TA&);

template<class Tx>
inline bool sv_v_not(const Tx&);

template <class TA>
inline void sp_m_resize(TA& A) throw();

#if(CXSC_INDEX_CHECK)
template <class TA>
inline void sp_m_resize(TA &A,const int &m, const int &n) throw(WRONG_BOUNDARIES);
#else
template <class TA>
inline void sp_m_resize(TA &A,const int &m, const int &n) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA>
inline void sp_m_resize(TA &A,const int &m1, const int &m2,const int &n1,const int &n2) throw(WRONG_BOUNDARIES);
#else
template<class TA>
inline void sp_m_resize(TA &A,const int &m1, const int &m2,const int &n1,const int &n2) throw();
#endif


}

#endif
