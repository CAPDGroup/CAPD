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

/* CVS $Id: imatrix.cpp,v 1.29 2014/01/30 17:23:45 cxsc Exp $ */

#define _CXSC_CPP

#include "imatrix.hpp"
#include "matrix.inl"
#include "imatrix.inl"
#include "ivecrmat.inl"

#include "idotk.inl"

namespace cxsc {

//Ostrowskis comparison matrix
rmatrix CompMat ( const imatrix& A) {
  rmatrix M(Lb(A,1), Ub(A,1), Lb(A,2), Ub(A,2));

  for(int i=Lb(A,1) ; i<=Ub(A,1) ; i++) {
    for(int j=Lb(A,2) ; j<=Ub(A,2) ; j++) {
      if(i-Lb(A,1) == j-Lb(A,2))
        M[i][j] = AbsMin(A[i][j]);
      else
        M[i][j] = -AbsMax(A[i][j]);
    }
  }

  return M;
}


imatrix Id ( const imatrix& A )                    // Interval identity matrix
{                                                  //-------------------------
  int i,j;
  int lbi = Lb(A,1), ubi = Ub(A,1);
  int lbj = Lb(A,2), ubj = Ub(A,2);
  imatrix B(lbi,ubi,lbj,ubj);

  for (i = lbi; i <= ubi; i++)
    for (j = lbj; j <= ubj; j++)
      B[i][j] = interval( (i==j) ? 1.0 : 0.0 );
  return B;
}

imatrix transp ( const imatrix& A )                       // Transposed matrix
{                                                         //------------------
  int      n;
  imatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));
      
  for (n = Lb(A,1); n <= Ub(A,1); n++) Col(res,n) = Row(A,n);
  return res;
}

real MaxRelDiam ( const imatrix_subv& v )                    // Maximum relative diameter
{                                                 //--------------------------
  real r;
  int  i, l=Lb(v), u=Ub(v);

  r = 0.0;
  for (i = l; i <= u; i++)
    if (RelDiam(v[i]) > r) r = RelDiam(v[i]);
  return r;
}

// The 'DoubleSize' functions double the number of rows of a matrix
// or double the length of a vector preserving existing components.
//------------------------------------------------------------------

void DoubleSize ( imatrix& A )
{
  int n = Lb(A,1);
  Resize(A,n,2*Ub(A,1)-n+1,Lb(A,2),Ub(A,2));
}


	 void accumulate(idotprecision &dp, const imatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const imatrix_subv &, imatrix_subv &)"));
#endif
		addDot(dp,rv1,rv2);	
	}

	 void accumulate(idotprecision &dp, const ivector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const ivector &, imatrix_subv &)"));
#endif
		addDot(dp,rv1,rv2);	
	}

	 void accumulate(idotprecision &dp, const imatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const imatrix_subv &, ivector &)"));
#endif
		addDot(dp,rv1,rv2);	
	}


	 void accumulate(idotprecision &dp, const ivector_slice & sl1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const ivector_slice &, imatrix_subv &)"));
#endif
		addDot(dp,sl1,rv2);	
	}

	 void accumulate(idotprecision &dp, const imatrix_subv & rv1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const imatrix_subv &, ivector_slice &)"));
#endif
		addDot(dp,rv1,sl2);	
	}

	 void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const imatrix_subv &, const imatrix_subv &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const ivector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector &, const imatrix_subv &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const imatrix_subv &, const ivector &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const ivector_slice & sl1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector_slice &, const imatrix_subv &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,sl1,rv2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const imatrix_subv &, const ivector_slice &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,sl2);
		dp += tmp;
	}

	 void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const rmatrix_subv &, imatrix_subv &)"));
#endif
		addDot(dp,rv1,rv2);	
	}

	 void accumulate(idotprecision &dp, const rvector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const rvector &, imatrix_subv &)"));
#endif
		addDot(dp,rv1,rv2);	
	}

	 void accumulate(idotprecision &dp, const rvector_slice & sl1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const rvector_slice &, imatrix_subv &)"));
#endif
		addDot(dp,sl1,rv2);	
	}

	 void accumulate(idotprecision &dp, const imatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const imatrix_subv &, rmatrix_subv &)"));
#endif
		addDot(dp,rv1,rv2);	
	}

	 void accumulate(idotprecision &dp, const imatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const imatrix_subv &, rvector &)"));
#endif
		addDot(dp,rv1,rv2);	
	}

	 void accumulate(idotprecision &dp, const imatrix_subv & rv1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const imatrix_subv &, rvector_slice &)"));
#endif
		addDot(dp,rv1,sl2);	
	}

	 void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rmatrix_subv &, const imatrix_subv &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const rvector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector &, const imatrix_subv &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const rvector_slice & sl1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector_slice&, const imatrix_subv &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,sl1,rv2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const imatrix_subv &, const rmatrix_subv &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const imatrix_subv &, const rvector &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;
	}

	 void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const imatrix_subv &, const rvector_slice &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,sl2);
		dp += tmp;
	}

	void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const rmatrix_subv &, ivector &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

	void accumulate(idotprecision &dp, const ivector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const ivector &, rmatrix_subv &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

	void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rmatrix_subv &, const ivector &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;	
	}

	void accumulate(cidotprecision &dp, const ivector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector &, const rmatrix_subv &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;	
	}

	void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const ivector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const rmatrix_subv &, ivector_slice &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

	void accumulate(idotprecision &dp, const ivector_slice & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(idotprecision &, const ivector_slice &, rmatrix_subv &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

	void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const ivector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rmatrix_subv &, const ivector_slice &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;	
	}

	void accumulate(cidotprecision &dp, const ivector_slice & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector_slice &, const rmatrix_subv &)"));
#endif
		idotprecision tmp(0.0);
		tmp.set_k(dp.get_k());
		addDot(tmp,rv1,rv2);
		dp += tmp;	
	}


} // namespace cxsc

