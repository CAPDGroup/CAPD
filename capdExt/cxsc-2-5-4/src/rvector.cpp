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

/* CVS $Id: rvector.cpp,v 1.26 2014/01/30 17:23:48 cxsc Exp $ */

#define _CXSC_CPP

#include "rvector.hpp"
#include "vector.inl"
#include "rvector.inl"

#include "dotk.inl"


namespace cxsc {


	void accumulate(dotprecision &dp, const rvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(dotprecision&, const rvector &, const rvector &)"));
#endif
		addDot(dp,rv1,rv2); 
	}

	void accumulate_approx(dotprecision &dp, const rvector & rv1, const rvector &rv2) {
		addDot_op(dp,rv1,rv2);
	}


//	INLINE void accumulate(dotprecision &dp, const rvector & rv1, const rmatrix_subv &rv2) throw(OP_WITH_WRONG_DIM);
//	INLINE void accumulate(dotprecision &dp, const rmatrix_subv & rv1, const rvector &rv2) throw(OP_WITH_WRONG_DIM);
	void accumulate(dotprecision &dp,const rvector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(dotprecision&, const rvector_slice &, const rvector &)"));
#endif
		addDot(dp,sl,rv); 
	}

	void accumulate_approx(dotprecision &dp,const rvector_slice &sl,const rvector &rv) {
		addDot_op(dp,sl,rv);
	}

	void accumulate(dotprecision &dp,const rvector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(dotprecision&, const rvector &, const rvector_slice &)"));
#endif
		addDot(dp,sl,rv); 
	}

	void accumulate_approx(dotprecision &dp,const rvector &rv,const rvector_slice &sl) {
		addDot_op(dp,rv,sl);
	}

	void accumulate(dotprecision &dp, const rvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(dotprecision&, const rvector_slice &, const rvector_slice &)"));
#endif
		addDot(dp,sl1,sl2); 
	}

	void accumulate_approx(dotprecision &dp, const rvector_slice & sl1, const rvector_slice &sl2) {
		addDot_op(dp,sl1,sl2);
	}


	void accumulate(idotprecision &dp, const rvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision&, const rvector &, const rvector &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}

//	INLINE void accumulate(idotprecision &dp, const rvector & rv1, const rmatrix_subv &rv2) throw(OP_WITH_WRONG_DIM);
//	INLINE void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const rvector &rv2) throw(OP_WITH_WRONG_DIM);


	void accumulate(idotprecision &dp,const rvector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision&, const rvector_slice &, const rvector &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,sl,rv);
		dp += tmp;
	}


	void accumulate(idotprecision &dp,const rvector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision&, const rvector &, const rvector_slice &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv,sl);
		dp += tmp;
	}

	void accumulate(idotprecision &dp, const rvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision&, const rvector_slice &, const rvector_slice &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,sl1,sl2);
		dp += tmp;
	}

	void accumulate(cdotprecision &dp, const rvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const rvector &, const rvector &)"));
#endif
		addDot(Re(dp),rv1,rv2); 
	}

	void accumulate_approx(cdotprecision &dp, const rvector & rv1, const rvector &rv2)
	{ 
		addDot_op(Re(dp),rv1,rv2); 
	}

//	INLINE void accumulate(cdotprecision &dp, const rvector & rv1, const rmatrix_subv &rv2) throw(OP_WITH_WRONG_DIM);
//	INLINE void accumulate(cdotprecision &dp, const rmatrix_subv & rv1, const rvector &rv2) throw(OP_WITH_WRONG_DIM);


	void accumulate(cdotprecision &dp,const rvector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const rvector_slice &, const rvector &)"));
#endif
		addDot(Re(dp),sl,rv);
	}

	void accumulate_approx(cdotprecision &dp,const rvector_slice &sl,const rvector &rv)
	{ 
		addDot_op(Re(dp),sl,rv);
	}


	void accumulate(cdotprecision &dp,const rvector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const rvector &, const rvector_slice &)"));
#endif
		addDot(Re(dp),rv,sl);
	}

	void accumulate_approx(cdotprecision &dp,const rvector &rv,const rvector_slice &sl)
	{ 
		addDot_op(Re(dp),rv,sl);
	}


	void accumulate(cdotprecision &dp, const rvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const rvector_slice &, const rvector_slice &)"));
#endif
		addDot(Re(dp),sl1,sl2);
	}

	void accumulate_approx(cdotprecision &dp, const rvector_slice & sl1, const rvector_slice &sl2)
	{ 
		addDot_op(Re(dp),sl1,sl2);
	}


	void accumulate(cidotprecision &dp, const rvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector &, const rvector &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}

//	INLINE void accumulate(cidotprecision &dp, const rvector & rv1, const rmatrix_subv &rv2) throw(OP_WITH_WRONG_DIM);
//	INLINE void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const rvector &rv2) throw(OP_WITH_WRONG_DIM);


	void accumulate(cidotprecision &dp,const rvector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector_slice &, const rvector &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,sl,rv);
		dp += tmp;
	}


	void accumulate(cidotprecision &dp,const rvector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector &, const rvector_slice &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv,sl);
		dp += tmp;
	}


	void accumulate(cidotprecision &dp, const rvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector_slice &, const rvector_slice &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,sl1,sl2);
		dp += tmp;
	}


	//Summation accumulates
	void accumulate(dotprecision &dp, const rvector& v) {
                addSum(dp,v);
        }

	void accumulate(idotprecision &dp, const rvector& v) {
                dotprecision tmp(0.0);
                tmp.set_k(dp.get_k());   
                addSum(tmp,v);
                Inf(dp) += tmp;
                Sup(dp) += tmp;
        }

	void accumulate(cdotprecision &dp, const rvector& v) {
                addSum(Re(dp),v);
        }

	void accumulate(cidotprecision &dp, const rvector& v) {
                dotprecision tmp(0.0);
                tmp.set_k(dp.get_k());   
                addSum(tmp,v);
                InfRe(dp) += tmp;
                SupRe(dp) += tmp;
        }
}
