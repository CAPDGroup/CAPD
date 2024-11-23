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

/* CVS $Id: civecrmat.hpp,v 1.24 2014/01/30 17:23:44 cxsc Exp $ */

// Here are definitions for civector x rmatrix-Functions
#ifndef _CXSC_CIVECRMAT_HPP_INCLUDED
#define _CXSC_CIVECRMAT_HPP_INCLUDED

namespace cxsc {

	INLINE civector _civector(const rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	INLINE civector _civector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif

	void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	void accumulate(cidotprecision &dp, const civector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const civector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	void accumulate(cidotprecision &dp, const civector_slice & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif


	INLINE civector operator *(const rmatrix &m,const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	INLINE civector operator *(const rmatrix_slice &ms,const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	INLINE civector operator *(const civector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	INLINE civector operator *(const civector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	INLINE civector &operator *=(civector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	INLINE civector &operator *=(civector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	INLINE civector operator *(const civector_slice &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

} // namespace cxsc 

#endif

