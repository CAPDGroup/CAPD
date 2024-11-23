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

/* CVS $Id: liveclrmat.hpp,v 1.25 2014/01/30 17:23:47 cxsc Exp $ */

// Here are definitions for l_ivector x l_rmatrix-Functions
#ifndef _CXSC_LIVECLRMAT_HPP_INCLUDED
#define _CXSC_LIVECLRMAT_HPP_INCLUDED

namespace cxsc {

	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_ivector _l_ivector(const l_rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_ivector _l_ivector(const l_rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const l_ivector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Sets the infimum of the vector
	INLINE void SetInf(l_ivector &iv,const l_rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
        //! Sets the supremum of the vector
	INLINE void SetSup(l_ivector &iv,const l_rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Sets the infimum of the vector
	INLINE void SetInf(l_ivector_slice &iv,const l_rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
        //! Sets the supremum of the vector
	INLINE void SetSup(l_ivector_slice &iv,const l_rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

        //! Sets the unchecked infimum of the vector
	INLINE void UncheckedSetInf(l_ivector &iv,const l_rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
        //! Sets the unchecked supremum of the vector
	INLINE void UncheckedSetSup(l_ivector &iv,const l_rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
        //! Sets the unchecked infimum of the vector
	INLINE void UncheckedSetInf(l_ivector_slice &iv,const l_rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
        //! Sets the unchecked supremum of the vector
	INLINE void UncheckedSetSup(l_ivector_slice &iv,const l_rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif


	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_rmatrix &m,const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_rmatrix_slice &ms,const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector &v,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector &operator *=(l_ivector &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector &operator *=(l_ivector &v,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector_slice &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const ivector &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const ivector &v,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const ivector_slice &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_rmatrix &m,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_rmatrix_slice &ms,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

} // namespace cxsc 

#endif

