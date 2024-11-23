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

/* CVS $Id: vector_friend_declarations.inl,v 1.10 2014/01/30 17:23:49 cxsc Exp $ */


#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres spsp_vv_mult(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres spsp_vv_mult(const Tx&, const Ty&) throw();
#endif


#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres slsp_vv_mult(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres slsp_vv_mult(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres spsl_vv_mult(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres spsl_vv_mult(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres spf_vv_mult(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres spf_vv_mult(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres fsp_vv_mult(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres fsp_vv_mult(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres slf_vv_mult(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres slf_vv_mult(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres fsl_vv_mult(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres fsl_vv_mult(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres slsl_vv_mult(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TDot>
friend inline Tres slsl_vv_mult(const Tx&, const Ty&) throw();
#endif

template<class Tv, class Ts, class Tres>
friend inline Tres sp_vs_div(const Tv&, const Ts&);

template<class Tv, class Ts, class Tres>
friend inline Tres sp_vs_mult(const Tv&, const Ts&);

template<class Ts, class Tv, class Tres>
friend inline Tres sp_sv_mult(const Ts&, const Tv&);

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsp_vv_add(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsp_vv_add(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsp_vv_add(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsp_vv_add(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsl_vv_add(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsl_vv_add(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres spf_vv_add(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres spf_vv_add(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres fsp_vv_add(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres fsp_vv_add(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres slf_vv_add(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres slf_vv_add(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres fsl_vv_add(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres fsl_vv_add(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsl_vv_add(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsl_vv_add(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsp_vv_sub(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsp_vv_sub(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsp_vv_sub(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsp_vv_sub(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsl_vv_sub(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsl_vv_sub(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres spf_vv_sub(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres spf_vv_sub(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres fsp_vv_sub(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres fsp_vv_sub(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres slf_vv_sub(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres slf_vv_sub(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres fsl_vv_sub(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres fsl_vv_sub(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsl_vv_sub(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsl_vv_sub(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsp_vv_hull(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsp_vv_hull(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsp_vv_hull(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsp_vv_hull(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsl_vv_hull(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsl_vv_hull(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres spf_vv_hull(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres spf_vv_hull(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres fsp_vv_hull(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres fsp_vv_hull(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres slf_vv_hull(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres slf_vv_hull(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres fsl_vv_hull(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres fsl_vv_hull(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsl_vv_hull(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsl_vv_hull(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsp_vv_intersect(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsp_vv_intersect(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsp_vv_intersect(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsp_vv_intersect(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsl_vv_intersect(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres spsl_vv_intersect(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres spf_vv_intersect(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres spf_vv_intersect(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres fsp_vv_intersect(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres fsp_vv_intersect(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres slf_vv_intersect(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres slf_vv_intersect(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres>
friend inline Tres fsl_vv_intersect(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres>
friend inline Tres fsl_vv_intersect(const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsl_vv_intersect(const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class Tres, class TType>
friend inline Tres slsl_vv_intersect(const Tx&, const Ty&) throw();
#endif

template<class Tx, class Ty>
friend inline bool spsp_vv_comp(const Tx&, const Ty&);

template<class Tx, class Ty>
friend inline bool slsp_vv_comp(const Tx&, const Ty&);

template<class Tx, class Ty>
friend inline bool spsl_vv_comp(const Tx&, const Ty&);

template<class Tx, class Ty>
friend inline bool spf_vv_comp(const Tx&, const Ty&);

template<class Tx, class Ty>
friend inline bool fsp_vv_comp(const Tx&, const Ty&);

template<class Tx, class Ty>
friend inline bool slf_vv_comp(const Tx&, const Ty&);

template<class Tx, class Ty>
friend inline bool fsl_vv_comp(const Tx&, const Ty&);

template<class Tx, class Ty>
friend inline bool slsl_vv_comp(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool spsp_vv_less(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool slsp_vv_less(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool spsl_vv_less(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool spf_vv_less(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool fsp_vv_less(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool slf_vv_less(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool fsl_vv_less(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool slsl_vv_less(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool spsp_vv_leq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool slsp_vv_leq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool spsl_vv_leq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool spf_vv_leq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool fsp_vv_leq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool slf_vv_leq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool fsl_vv_leq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool slsl_vv_leq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool spsp_vv_greater(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool slsp_vv_greater(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool spsl_vv_greater(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool spf_vv_greater(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool fsp_vv_greater(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool slf_vv_greater(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool fsl_vv_greater(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool slsl_vv_greater(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool spsp_vv_geq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool slsp_vv_geq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool spsl_vv_geq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool spf_vv_geq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool fsp_vv_geq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool slf_vv_geq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool fsl_vv_geq(const Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline bool slsl_vv_geq(const Tx&, const Ty&);

template<class Tx, class Ts, class TType>
friend inline Tx& sp_vs_assign(Tx&, const Ts&);

template<class Tx, class Ts, class TType, class TIt>
friend inline Tx& sl_vs_assign(Tx&, const Ts&);

template<class Tx, class Ty, class TType>
friend inline Tx& spf_vv_assign(Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline Tx& spsl_vv_assign(Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline Tx& fsp_vv_assign(Tx&, const Ty&);

template<class Tx, class Ty, class TType>
friend inline Tx& fsl_vv_assign(Tx&, const Ty&);

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class TType>
friend inline Tx& fssp_vv_assign(Tx& v1, const Ty& v2) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class TType>
friend inline Tx& fssp_vv_assign(Tx& v1, const Ty& v2) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class TType>
friend inline Tx& fssl_vv_assign(Tx& v1, const Ty& v2) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class TType>
friend inline Tx& fssl_vv_assign(Tx& v1, const Ty& v2) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class TType, class TIt>
friend inline Tx& slsl_vv_assign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class TType, class TIt>
friend inline Tx& slsl_vv_assign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class TType, class TIt>
friend inline Tx& slsp_vv_assign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class TType, class TIt>
friend inline Tx& slsp_vv_assign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class TType, class TIt>
friend inline Tx& slf_vv_assign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class TType, class TIt>
friend inline Tx& slf_vv_assign(Tx&, const Ty&) throw();
#endif

template<class Tx, class TType>
friend inline std::ostream& sp_v_output(std::ostream&, const Tx&);

template<class Tx, class TType>
friend inline std::ostream& sl_v_output(std::ostream&, const Tx&);

template<class Tx, class TType>
friend inline std::istream& sp_v_input(std::istream&, Tx&);

template<class Tx, class TType>
friend inline std::istream& sl_v_input(std::istream&, Tx&);

template<class Tx>
friend inline Tx sp_v_negative(const Tx&);

template<class Tx, class Tres>
friend inline Tres sl_v_negative(const Tx&);

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& spf_vv_addassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& spf_vv_addassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& spsp_vv_addassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& spsp_vv_addassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& spsl_vv_addassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& spsl_vv_addassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& fsp_vv_addassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& fsp_vv_addassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& fsl_vv_addassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& fsl_vv_addassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& spf_vv_subassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& spf_vv_subassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& spsp_vv_subassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& spsp_vv_subassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& spsl_vv_subassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& spsl_vv_subassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& fsp_vv_subassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& fsp_vv_subassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& fsl_vv_subassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& fsl_vv_subassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& spf_vv_hullassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& spf_vv_hullassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& spsp_vv_hullassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& spsp_vv_hullassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& spsl_vv_hullassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& spsl_vv_hullassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& fsp_vv_hullassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& fsp_vv_hullassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& fsl_vv_hullassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& fsl_vv_hullassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& spf_vv_intersectassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& spf_vv_intersectassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& spsp_vv_intersectassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& spsp_vv_intersectassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& spsl_vv_intersectassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& spsl_vv_intersectassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& fsp_vv_intersectassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& fsp_vv_intersectassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& fsl_vv_intersectassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& fsl_vv_intersectassign(Tx&, const Ty&) throw();
#endif

template<class Tx, class Ty>
friend inline Tx& sp_vs_multassign(Tx& v, const Ty& s);

template<class Tx, class Ty>
friend inline Tx& sp_vs_divassign(Tx& v, const Ty& s);

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class TType>
friend inline Tx& slf_vv_addassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class TType>
friend inline Tx& slf_vv_addassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& slsp_vv_addassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& slsp_vv_addassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& slsl_vv_addassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& slsl_vv_addassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class TType>
friend inline Tx& slf_vv_subassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class TType>
friend inline Tx& slf_vv_subassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& slsp_vv_subassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& slsp_vv_subassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& slsl_vv_subassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& slsl_vv_subassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class TType>
friend inline Tx& slf_vv_hullassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class TType>
friend inline Tx& slf_vv_hullassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& slsp_vv_hullassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& slsp_vv_hullassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& slsl_vv_hullassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& slsl_vv_hullassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty, class TType>
friend inline Tx& slf_vv_intersectassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty, class TType>
friend inline Tx& slf_vv_intersectassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& slsp_vv_intersectassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& slsp_vv_intersectassign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& slsl_vv_intersectassign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& slsl_vv_intersectassign(Tx&, const Ty&) throw();
#endif

template<class Tx, class Ty>
friend inline Tx& sl_vs_multassign(Tx&, const Ty&);

template<class Tx, class Ty>
friend inline Tx& sl_vs_divassign(Tx&, const Ty&);

#if(CXSC_INDEX_CHECK)
template<class TA, class Tx, class Tres, class TDot, class TElement>
friend inline Tres spsl_mv_mult(const TA&, const Tx&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class Tx, class Tres, class TDot, class TElement>
friend inline Tres spsl_mv_mult(const TA&, const Tx&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class Tx, class Tres, class TDot, class TElement>
friend inline Tres spsp_mv_mult(const TA&, const Tx&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class Tx, class Tres, class TDot, class TElement>
friend inline Tres spsp_mv_mult(const TA&, const Tx&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class Tx, class Tres, class TDot>
friend inline Tres spf_mv_mult(const TA&, const Tx&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class Tx, class Tres, class TDot>
friend inline Tres spf_mv_mult(const TA&, const Tx&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class Tx, class Tres, class TDot>
friend inline Tres fsp_mv_mult(const TA&, const Tx&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class Tx, class Tres, class TDot>
friend inline Tres fsp_mv_mult(const TA&, const Tx&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TA, class Tx, class Tres, class TDot>
friend inline Tres fsl_mv_mult(const TA&, const Tx&) throw(OP_WITH_WRONG_DIM);
#else
template<class TA, class Tx, class Tres, class TDot>
friend inline Tres fsl_mv_mult(const TA&, const Tx&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& svsp_vv_assign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& svsp_vv_assign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& svsl_vv_assign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& svsl_vv_assign(Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx, class Ty>
friend inline Tx& svf_vv_assign(Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class Tx, class Ty>
friend inline Tx& svf_vv_assign(Tx&, const Ty&) throw();
#endif

template<class Tx, class Ts>
friend inline Tx& sv_vs_assign(Tx&, const Ts&);

template<class Tx>
friend inline bool sv_v_not(const Tx&);

template<class Tx>
friend inline bool sp_v_not(const Tx&);

template<class Tx>
friend inline bool sl_v_not(const Tx&);

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void spsp_vv_accu(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void spsp_vv_accu(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void spf_vv_accu(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void spf_vv_accu(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void fsp_vv_accu(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void fsp_vv_accu(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void slsl_vv_accu(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void slsl_vv_accu(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void spsl_vv_accu(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void spsl_vv_accu(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void slsp_vv_accu(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void slsp_vv_accu(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void slf_vv_accu(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void slf_vv_accu(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void fsl_vv_accu(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void fsl_vv_accu(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void spsp_vv_accuapprox(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void spsp_vv_accuapprox(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void spf_vv_accuapprox(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void spf_vv_accuapprox(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void fsp_vv_accuapprox(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void fsp_vv_accuapprox(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void slsl_vv_accuapprox(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void slsl_vv_accuapprox(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void spsl_vv_accuapprox(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void spsl_vv_accuapprox(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void slsp_vv_accuapprox(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void slsp_vv_accuapprox(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void slf_vv_accuapprox(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void slf_vv_accuapprox(TDot&, const Tx&, const Ty&) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void fsl_vv_accuapprox(TDot&, const Tx&, const Ty&) throw(OP_WITH_WRONG_DIM);
#else
template<class TDot, class Tx, class Ty, class TSparseDot>
friend inline void fsl_vv_accuapprox(TDot&, const Tx&, const Ty&) throw();
#endif

template<class Tx>
friend inline void sp_v_resize(Tx &v) throw();

#if(CXSC_INDEX_CHECK)
template <class Tx>
friend inline void sp_v_resize(Tx &v, const int &len) throw(WRONG_BOUNDARIES);
#else
template <class Tx>
friend inline void sp_v_resize(Tx &v, const int &len) throw();
#endif

#if(CXSC_INDEX_CHECK)
template<class Tx>
friend inline void sp_v_resize(Tx &v, const int &lb, const int &ub) throw(WRONG_BOUNDARIES);
#else
template<class Tx>
friend inline void sp_v_resize(Tx &v, const int &lb, const int &ub) throw();
#endif
