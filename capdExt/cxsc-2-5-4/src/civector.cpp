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

/* CVS $Id: civector.cpp,v 1.26 2014/01/30 17:23:44 cxsc Exp $ */

#define _CXSC_CPP

#include "civector.hpp"
#include "vector.inl"
#include "civector.inl"
#include "rmatrix.hpp"

#include "iveccvec.inl"

#include "cidotk.inl"


namespace cxsc {


	 void accumulate(cidotprecision &dp, const civector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector &, const civector &)"));
#endif
		addDot(dp,rv1,rv2);
	}

	void accumulate(cidotprecision &dp, const civector_slice & sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice &, const civector &)"));
#endif
		addDot(dp,sl,rv);
	}

	void accumulate(cidotprecision &dp, const civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector &, const civector_slice &)"));
#endif
		addDot(dp,rv,sl);
	}

	void accumulate(cidotprecision &dp, const civector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice &, const civector_slice &)"));
#endif
		addDot(dp,sl1,sl2);
	}

	void accumulate(cidotprecision &dp, const rvector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector &, const civector &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,rv1,Re(rv2));
		addDot(tmp_im,rv1,Im(rv2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp, const civector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector &, const rvector &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv1),rv2);
		addDot(tmp_im,Im(rv1),rv2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const rvector_slice & sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector_slice &, const civector &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,sl,Re(rv));
		addDot(tmp_im,sl,Im(rv));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp,const civector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice&, const rvector &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(sl),rv);
		addDot(tmp_im,Im(sl),rv);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const rvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector &, const civector_slice &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,rv,Re(sl));
		addDot(tmp_im,rv,Im(sl));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp,const civector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector &, const rvector_slice &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv),sl);
		addDot(tmp_im,Im(rv),sl);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const civector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice &, const rvector_slice &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(sl1),sl2);
		addDot(tmp_im,Im(sl1),sl2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const rvector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector_slice &, const civector_slice &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,sl1,Re(sl2));
		addDot(tmp_im,sl1,Im(sl2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const civector_slice & sl1, const rmatrix_subv &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice &, const rmatrix_subv &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(sl1),sl2);
		addDot(tmp_im,Im(sl1),sl2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const rmatrix_subv & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rmatrix_subv &, const civector_slice &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,sl1,Re(sl2));
		addDot(tmp_im,sl1,Im(sl2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp, const cvector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector &, const civector &)"));
#endif
		addDot(dp,rv1,rv2);
	}

	 void accumulate(cidotprecision &dp, const civector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector &, const cvector &)"));
#endif
		addDot(dp,rv1,rv2);
	}

	void accumulate(cidotprecision &dp, const cvector_slice & sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector_slice &, const civector &)"));
#endif
		addDot(dp,sl,rv);
	}

	void accumulate(cidotprecision &dp,const civector_slice &sl,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice &, const cvector &)"));
#endif
		addDot(dp,sl,rv);
	}

	void accumulate(cidotprecision &dp, const cvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector &, const civector_slice &)"));
#endif
		addDot(dp,rv,sl);
	}

	void accumulate(cidotprecision &dp,const civector &rv,const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector &, const cvector_slice &)"));
#endif
		addDot(dp,rv,sl);
	}

	void accumulate(cidotprecision &dp, const civector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice &, const cvector_slice &)"));
#endif
		addDot(dp,sl1,sl2);
	}

	 void accumulate(cidotprecision &dp, const cvector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector_slice &, const civector_slice &)"));
#endif
		addDot(dp,sl1,sl2);
	}

	void accumulate(cidotprecision &dp, const ivector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector &, const civector &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,rv1,Re(rv2));
		addDot(tmp_im,rv1,Im(rv2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const civector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector &, const ivector &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv1),rv2);
		addDot(tmp_im,Im(rv1),rv2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp, const ivector_slice & sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector_slice &, const civector &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,sl,Re(rv));
		addDot(tmp_im,sl,Im(rv));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp,const civector_slice &sl,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice &, const ivector &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(sl),rv);
		addDot(tmp_im,Im(sl),rv);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const ivector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector &, const civector_slice &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,rv,Re(sl));
		addDot(tmp_im,rv,Im(sl));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp,const civector &rv,const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector &, const ivector_slice &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv),sl);
		addDot(tmp_im,Im(rv),sl);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const civector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice &, const ivector_slice &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(sl1),sl2);
		addDot(tmp_im,Im(sl1),sl2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp, const ivector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector_slice &, const civector_slice &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,sl1,Re(sl2));
		addDot(tmp_im,sl1,Im(sl2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

        //Summation
	void accumulate(cidotprecision &dp, const civector& v) {
                addSum(InfRe(dp),InfRe(v));
                addSum(SupRe(dp),SupRe(v));
                addSum(InfIm(dp),InfIm(v));
                addSum(SupIm(dp),SupIm(v));
        }
} // namespace cxsc

