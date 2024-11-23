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

/* CVS $Id: iveccmat.hpp,v 1.25 2014/01/30 17:23:45 cxsc Exp $ */

// Here are definitions for ivector x cmatrix-Functions
#ifndef _CXSC_IVECCMAT_HPP_INCLUDED
#define _CXSC_IVECCMAT_HPP_INCLUDED

namespace cxsc {

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const ivector & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const ivector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const ivector_slice & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

} // namespace cxsc 
	
#endif

