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

/* CVS $Id: lrvecimat.inl,v 1.23 2014/01/30 17:23:47 cxsc Exp $ */

#ifndef _CXSC_LRVECIMAT_INL_INCLUDED
#define _CXSC_LRVECIMAT_INL_INCLUDED

namespace cxsc {

	INLINE void accumulate(idotprecision &dp, const imatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<idotprecision,l_rvector,imatrix_subv>(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<idotprecision,l_rvector,imatrix_subv>(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const imatrix_subv & rv1, const l_rvector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<idotprecision,l_rvector,imatrix_subv>(dp,l_rvector(rv2),rv1); }
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<idotprecision,l_rvector,imatrix_subv>(dp,l_rvector(rv1),rv2); }

} // namespace cxsc

#endif

