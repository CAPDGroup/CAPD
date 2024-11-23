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

/* CVS $Id: cvector.cpp,v 1.25 2014/01/30 17:23:44 cxsc Exp $ */

#define _CXSC_CPP

#include "cvector.hpp"
#include "vector.inl"
#include "cvector.inl"

#include "ivector.hpp"

#include "idotk.inl"
#include "cdotk.inl"

namespace cxsc {

	void accumulate(cdotprecision &dp, const cvector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const cvector &, const cvector &)"));
#endif
		addDot(dp,rv1,rv2);
	}

	void accumulate_approx(cdotprecision &dp, const cvector & rv1, const cvector &rv2)
	{ 
		addDot_op(dp,rv1,rv2);
	}


	void accumulate(cdotprecision &dp, const cvector_slice & sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const cvector_slice &, const cvector &)"));
#endif
		addDot(dp,sl,rv);
	}

	void accumulate_approx(cdotprecision &dp, const cvector_slice & sl, const cvector &rv)
	{ 
		addDot_op(dp,sl,rv);
	}


	void accumulate(cdotprecision &dp, const cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const cvector &, const cvector_slice &)"));
#endif
		addDot(dp,rv,sl);
	}

	void accumulate_approx(cdotprecision &dp, const cvector &rv, const cvector_slice &sl)
	{ 
		addDot_op(dp,rv,sl);
	}


	void accumulate(cdotprecision &dp, const cvector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const cvector_slice &, const cvector_slice &)"));
#endif
		addDot(dp,sl1,sl2);
	}

	void accumulate_approx(cdotprecision &dp, const cvector_slice & sl1, const cvector_slice &sl2)
	{ 
		addDot_op(dp,sl1,sl2);
	}

	void accumulate(cidotprecision &dp, const cvector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector &, const cvector &)"));
#endif
		cdotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const cvector_slice & sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector_slice &, const cvector &)"));
#endif
		cdotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,sl,rv);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector &, const cvector_slice &)"));
#endif
		cdotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv,sl);
		dp += tmp;
	}

	void accumulate(cidotprecision &dp, const cvector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector_slice &, const cvector_slice &)"));
#endif
		cdotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,sl1,sl2);
		dp += tmp;
	}

	void accumulate(cdotprecision &dp, const cvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const cvector &, const rvector &)"));
#endif	
		addDot(Re(dp), Re(rv1), rv2);
		addDot(Im(dp), Im(rv1), rv2);
	}

	void accumulate_approx(cdotprecision &dp, const cvector & rv1, const rvector &rv2)
	{ 
		addDot_op(Re(dp), Re(rv1), rv2);
		addDot_op(Im(dp), Im(rv1), rv2);
	}

	 void accumulate(cdotprecision &dp, const rvector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const rvector &, const cvector &)"));
#endif	
		addDot(Re(dp), rv1, Re(rv2));
		addDot(Im(dp), rv1, Im(rv2));
	}

	 void accumulate_approx(cdotprecision &dp, const rvector & rv1, const cvector &rv2)
	{ 
		addDot_op(Re(dp), rv1, Re(rv2));
		addDot_op(Im(dp), rv1, Im(rv2));
	}

	 void accumulate(cdotprecision &dp, const rvector_slice & sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const rvector_slice &, const cvector &)"));
#endif	
		addDot(Re(dp), sl, Re(rv));
		addDot(Im(dp), sl, Im(rv));
	}

	void accumulate_approx(cdotprecision &dp, const rvector_slice & sl, const cvector &rv)
	{ 
		addDot_op(Re(dp), sl, Re(rv));
		addDot_op(Im(dp), sl, Im(rv));
	}

	 void accumulate(cdotprecision &dp,const cvector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const cvector_slice &, const rvector &)"));
#endif	
		addDot(Re(dp), Re(sl), rv);
		addDot(Im(dp), Im(sl), rv);
	}

	void accumulate_approx(cdotprecision &dp,const cvector_slice &sl,const rvector &rv)
	{ 
		addDot_op(Re(dp), Re(sl), rv);
		addDot_op(Im(dp), Im(sl), rv);
	}

	void accumulate(cdotprecision &dp, const rvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const rvector &, const cvector_slice &)"));
#endif	
		addDot(Re(dp), rv, Re(sl));
		addDot(Im(dp), rv, Im(sl));
	}

	void accumulate_approx(cdotprecision &dp, const rvector &rv, const cvector_slice &sl)
	{ 
		addDot_op(Re(dp), rv, Re(sl));
		addDot_op(Im(dp), rv, Im(sl));
	}

	void accumulate(cdotprecision &dp,const cvector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const cvector &, const rvector_slice &)"));
#endif	
		addDot(Re(dp), Re(rv), sl);
		addDot(Im(dp), Im(rv), sl);
	}

	void accumulate_approx(cdotprecision &dp,const cvector &rv,const rvector_slice &sl)
	{ 
		addDot_op(Re(dp), Re(rv), sl);
		addDot_op(Im(dp), Im(rv), sl);
	}

	 void accumulate(cdotprecision &dp, const cvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const cvector_slice &, const rvector_slice &)"));
#endif	
		addDot(Re(dp), Re(sl1), sl2);
		addDot(Im(dp), Im(sl1), sl2);
	}

	void accumulate_approx(cdotprecision &dp, const cvector_slice & sl1, const rvector_slice &sl2)
	{ 
		addDot_op(Re(dp), Re(sl1), sl2);
		addDot_op(Im(dp), Im(sl1), sl2);
	}

	void accumulate(cdotprecision &dp, const rvector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const rvector_slice &, const cvector_slice &)"));
#endif	
		addDot(Re(dp), sl1, Re(sl2));
		addDot(Im(dp), sl1, Im(sl2));
	}

	void accumulate_approx(cdotprecision &dp, const rvector_slice & sl1, const cvector_slice &sl2)
	{ 
		addDot_op(Re(dp), sl1, Re(sl2));
		addDot_op(Im(dp), sl1, Im(sl2));
	}

	 void accumulate(cidotprecision &dp, const cvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector &, const rvector &)"));
#endif
		cdotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(Re(tmp),Re(rv1),rv2);
		addDot(Im(tmp),Im(rv1),rv2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const rvector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector &, const cvector &)"));
#endif
		cdotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(Re(tmp),rv1,Re(rv2));
		addDot(Im(tmp),rv1,Im(rv2));
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const rvector_slice & sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector_slice &, const cvector &)"));
#endif
		cdotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(Re(tmp),sl,Re(rv));
		addDot(Im(tmp),sl,Im(rv));
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp,const cvector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector_slice &, const rvector &)"));
#endif
		cdotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(Re(tmp),Re(sl),rv);
		addDot(Im(tmp),Im(sl),rv);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const rvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector &, const cvector_slice &)"));
#endif
		cdotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(Re(tmp),rv,Re(sl));
		addDot(Im(tmp),rv,Im(sl));
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp,const cvector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv)!=VecLen(sl)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector &, const rvector_slice &)"));
#endif
		cdotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(Re(tmp),Re(rv),sl);
		addDot(Im(tmp),Im(rv),sl);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const cvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector_slice &, const rvector_slice &)"));
#endif
		cdotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(Re(tmp),Re(sl1),sl2);
		addDot(Im(tmp),Im(sl1),sl2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const rvector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector_slice &, const cvector_slice &)"));
#endif
		cdotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(Re(tmp),sl1,Re(sl2));
		addDot(Im(tmp),sl1,Im(sl2));
		dp += tmp;
	}

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector_slice &, const ivector &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(sl),rv);
		addDot(tmp_im,Im(sl),rv);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector & rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector &, const cvector_slice &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,rv,Re(sl));
		addDot(tmp_im,rv,Im(sl));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector &, const ivector &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv1),rv2);
		addDot(tmp_im,Im(rv1),rv2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector &, const cvector &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,rv1,Re(rv2));
		addDot(tmp_im,rv1,Im(rv2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector_slice &, const cvector &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,sl,Re(rv));
		addDot(tmp_im,sl,Im(rv));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(rv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector &, const ivector_slice &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv),sl);
		addDot(tmp_im,Im(rv),sl);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector_slice &, const ivector_slice &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(sl1),sl2);
		addDot(tmp_im,Im(sl1),sl2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector_slice &, const cvector_slice &)"));
#endif
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,sl1,Re(sl2));
		addDot(tmp_im,sl1,Im(sl2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	//Summation accumulates
	void accumulate(cdotprecision &dp, const cvector& v) {
                addSum(Re(dp),Re(v));
                addSum(Im(dp),Im(v));
        }

	void accumulate(cidotprecision &dp, const cvector& v) {
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
                accumulate(tmp_re,Re(v));
                accumulate(tmp_im,Im(v));
                dp += cidotprecision(tmp_re,tmp_im);
        }


} // namespace cxsc

