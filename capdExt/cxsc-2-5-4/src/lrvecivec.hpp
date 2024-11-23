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

/* CVS $Id: lrvecivec.hpp,v 1.25 2014/01/30 17:23:47 cxsc Exp $ */

// Here are definitions for imatrix x cmatrix-Functions
#ifndef _CXSC_LRVECIVEC_HPP_INCLUDED
#define _CXSC_LRVECIVEC_HPP_INCLUDED

#include "l_interval.hpp"

namespace cxsc {

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const ivector & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp,const ivector_slice &sl,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp,const ivector &rv,const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const ivector_slice & sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_rvector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_rvector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const ivector & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const ivector_slice & sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

} // namespace cxsc 

#endif

