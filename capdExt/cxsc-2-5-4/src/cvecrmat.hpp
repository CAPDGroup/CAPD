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

/* CVS $Id: cvecrmat.hpp,v 1.24 2014/01/30 17:23:44 cxsc Exp $ */

// Here are definitions for cvector x rmatrix-Functions
#ifndef _CXSC_CVECRMAT_HPP_INCLUDED
#define _CXSC_CVECRMAT_HPP_INCLUDED

namespace cxsc {

	INLINE cvector _cvector(const rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	INLINE cvector _cvector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif

	void accumulate(cdotprecision &dp, const rmatrix_subv & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	void accumulate(cdotprecision &dp, const cvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	void accumulate(cidotprecision &dp, const cvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	void accumulate(cdotprecision &dp, const rmatrix_subv & rv1, const cvector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	void accumulate(cdotprecision &dp, const cvector_slice & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const cvector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	void accumulate(cidotprecision &dp, const cvector_slice & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	INLINE void SetIm(cvector &iv,const rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	INLINE void SetRe(cvector &iv,const rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	INLINE void SetIm(cvector_slice &iv,const rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	INLINE void SetRe(cvector_slice &iv,const rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif


	INLINE cvector operator *(const rmatrix &m,const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	INLINE cvector operator *(const rmatrix_slice &ms,const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	INLINE cvector operator *(const cvector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	INLINE cvector operator *(const cvector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	INLINE cvector &operator *=(cvector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	INLINE cvector &operator *=(cvector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	INLINE cvector operator *(const cvector_slice &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

} // namespace cxsc 

#endif

