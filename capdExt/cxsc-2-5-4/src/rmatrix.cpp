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

/* CVS $Id: rmatrix.cpp,v 1.28 2014/01/30 17:23:48 cxsc Exp $ */

#define _CXSC_CPP

#include "rmatrix.hpp"
#include "matrix.inl"
#include "rmatrix.inl"

#include "dotk.inl"

namespace cxsc {
  
    //Ostrowskis comparison matrix
    rmatrix CompMat ( const rmatrix& A) throw() {
      rmatrix M(Lb(A,1), Ub(A,1), Lb(A,2), Ub(A,2));

      for(int i=Lb(A,1) ; i<=Ub(A,1) ; i++) {
	for(int j=Lb(A,2) ; j<=Ub(A,2) ; j++) {
	  if(i-Lb(A,1) == j-Lb(A,2))
	    M[i][j] = abs(A[i][j]);
	  else
	    M[i][j] = -abs(A[i][j]);
	}
      }

  return M;
}

      rmatrix Id ( const rmatrix& A )                        // Real identity matrix
      {                                                      //---------------------
        int i,j;
        int lbi = Lb(A,1), ubi = Ub(A,1);
        int lbj = Lb(A,2), ubj = Ub(A,2);
        rmatrix B(lbi,ubi,lbj,ubj);
      
        for (i = lbi; i <= ubi; i++)
          for (j = lbj; j <= ubj; j++)
            B[i][j] = (i==j) ? 1.0 : 0.0;
        return B;
      }
      
      rmatrix transp ( const rmatrix& A )                       // Transposed matrix
      {                                                         //------------------
        int      n;
        rmatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));
      
        for (n = Lb(A,1); n <= Ub(A,1); n++) Col(res,n) = Row(A,n);
        return res;
      }

      void DoubleSize ( rmatrix& A )
      {
        int n = Lb(A,1);
        Resize(A,n,2*Ub(A,1)-n+1,Lb(A,2),Ub(A,2));
      }


	void accumulate(dotprecision &dp, const rmatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(dotprecision&, const rmatrix_subv &, const rmatrix_subv &)"));
#endif
		addDot(dp,rv1,rv2);
	}

	void accumulate_approx(dotprecision &dp, const rmatrix_subv & rv1, const rmatrix_subv &rv2) {
		addDot_op(dp,rv1,rv2);
	}


	 void accumulate(dotprecision &dp, const rvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(dotprecision&, const rvector &, const rmatrix_subv &)"));
#endif
		addDot(dp,rv1,rv2);
	}

	void accumulate_approx(dotprecision &dp, const rvector & rv1, const rmatrix_subv &rv2) {
		addDot_op(dp,rv1,rv2);
	}


	void accumulate(dotprecision &dp, const rmatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(dotprecision&, const rmatrix_subv &, const rvector &)"));
#endif
		addDot(dp,rv1,rv2);
	}

	void accumulate_approx(dotprecision &dp, const rmatrix_subv & rv1, const rvector &rv2) {
		addDot_op(dp,rv1,rv2);
	}


	 void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision&, const rmatrix_subv &, const rmatrix_subv &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}



	 void accumulate(idotprecision &dp, const rvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision&, const rvector&, const rmatrix_subv &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}



	 void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision&, const rmatrix_subv&, const rvector &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}



	 void accumulate(cdotprecision &dp, const rmatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const rmatrix_subv&, const rmatrix_subv &)"));
#endif
		addDot(Re(dp),rv1,rv2);
	}

	 void accumulate_approx(cdotprecision &dp, const rmatrix_subv & rv1, const rmatrix_subv &rv2)
	{ 
		addDot_op(Re(dp),rv1,rv2);
	}


	 void accumulate(cdotprecision &dp, const rvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const rvector&, const rmatrix_subv &)"));
#endif
		addDot(Re(dp),rv1,rv2);
	}

	 void accumulate_approx(cdotprecision &dp, const rvector & rv1, const rmatrix_subv &rv2)
	{ 
		addDot_op(Re(dp),rv1,rv2);
	}

	 void accumulate(cdotprecision &dp, const rmatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const rmatrix_subv&, const rvector &)"));
#endif
		addDot(Re(dp),rv1,rv2);
	}

	 void accumulate_approx(cdotprecision &dp, const rmatrix_subv & rv1, const rvector &rv2)
	{ 
		addDot_op(Re(dp),rv1,rv2);
	}

	 void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rmatrix_subv &, const rmatrix_subv &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const rvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector &, const rmatrix_subv &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rmatrix_subv &, const rvector &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}


	void accumulate(dotprecision &dp,const rvector_slice &sl,const rmatrix_subv &sv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl)!=VecLen(sv)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(dotprecision&, const rvector_slice &, const rmatrix_subv &)"));
#endif
		addDot(dp,sl,sv);
	}

	void accumulate_approx(dotprecision &dp,const rvector_slice &sl,const rmatrix_subv &sv) {
		addDot_op(dp,sl,sv);
	}


	 void accumulate(cdotprecision &dp, const rvector_slice & sl1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const rvector_slice&, const rmatrix_subv &)"));
#endif	
		addDot(Re(dp),sl1,rv2);
	}

	 void accumulate_approx(cdotprecision &dp, const rvector_slice & sl1, const rmatrix_subv &rv2)
	{ 
		addDot_op(Re(dp),sl1,rv2);
	}


	 void accumulate(idotprecision &dp, const rvector_slice & sl1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision&, const rvector_slice&, const rmatrix_subv &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,sl1,rv2);
		dp += tmp;
	}


	 void accumulate(cidotprecision &dp, const rvector_slice & sl1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector_slice &, const rmatrix_subv &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,sl1,rv2);
		dp += tmp;
	}


	void accumulate(dotprecision &dp,const rmatrix_subv &mv,const rvector_slice &vs)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(mv)!=VecLen(vs)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(dotprecision&, const rmatrix_subv &, const rvector_slice &)"));
#endif
		addDot(dp,vs,mv); 
	}

	void accumulate_approx(dotprecision &dp,const rmatrix_subv &mv,const rvector_slice &vs) {
		addDot_op(dp,vs,mv);
	}



	 void accumulate(cdotprecision &dp, const rmatrix_subv & rv1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cdotprecision&, const rmatrix_subv&, const rvector_slice &)"));
#endif	
		addDot(Re(dp),rv1,sl2);
	}

	 void accumulate_approx(cdotprecision &dp, const rmatrix_subv & rv1, const rvector_slice &sl2)
	{ 
		addDot_op(Re(dp),rv1,sl2);
	}


	 void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision&, const rmatrix_subv&, const rvector_slice &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,sl2);
		dp += tmp;
	}


	 void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rmatrix_subv &, const rvector_slice &)"));
#endif
		dotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,sl2);
		dp += tmp;
	}


} // namespace cxsc
