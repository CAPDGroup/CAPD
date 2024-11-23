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

/* CVS $Id: cimatrix.cpp,v 1.30 2014/01/30 17:23:43 cxsc Exp $ */

#define _CXSC_CPP

#include "cimatrix.hpp"
#include "matrix.inl"
#include "cimatrix.inl"
#include "civecrmat.inl"
#include "civeccmat.inl"
#include "civecimat.inl"

// they should also be included here
#include "cvecimat.inl"
#include "iveccmat.inl"
#include "cmatimat.inl"

#include "cidotk.inl"

namespace cxsc {

//Ostrowskis comparison matrix
rmatrix CompMat ( const cimatrix& A) {
  rmatrix M(Lb(A,1), Ub(A,1), Lb(A,2), Ub(A,2));

  for(int i=Lb(A,1) ; i<=Ub(A,1) ; i++) {
    for(int j=Lb(A,2) ; j<=Ub(A,2) ; j++) {
      if(i-Lb(A,1) == j-Lb(A,2))
        M[i][j] = Inf(abs(A[i][j]));
      else
        M[i][j] = -Sup(abs(A[i][j]));
    }
  }

  return M;
}


cimatrix Id ( const cimatrix& A )                  // Interval identity matrix
{                                                  //-------------------------
  int i,j;
  int lbi = Lb(A,1), ubi = Ub(A,1);
  int lbj = Lb(A,2), ubj = Ub(A,2);
  cimatrix B(lbi,ubi,lbj,ubj);

  for (i = lbi; i <= ubi; i++)
    for (j = lbj; j <= ubj; j++)
      B[i][j] = cinterval( (i==j) ? 1.0 : 0.0 );
  return B;
}

cimatrix transp ( const cimatrix& A )                     // Transposed matrix
{                                                         //------------------
  int      n;
  cimatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));
      
  for (n = Lb(A,1); n <= Ub(A,1); n++) Col(res,n) = Row(A,n);
  return res;
}

void DoubleSize ( cimatrix& A )
{
  int n = Lb(A,1);
  Resize(A,n,2*Ub(A,1)-n+1,Lb(A,2),Ub(A,2));
}


	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cimatrix_subv &, const cimatrix_subv &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

	 void accumulate(cidotprecision &dp, const civector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector &, const cimatrix_subv &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cimatrix_subv &, const civector &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

	 void accumulate(cidotprecision &dp, const civector_slice & sl1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice &, const cimatrix_subv &)"));
#endif	
		addDot(dp,sl1,rv2);
	}

	void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cimatrix_subv &, const civector_slice &)"));
#endif	
		addDot(dp,rv1,sl2);
	}
	
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cimatrix_subv &, const rmatrix_subv &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv1),rv2);
		addDot(tmp_im,Im(rv1),rv2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cimatrix_subv &, const rvector_slice &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv1),sl2);
		addDot(tmp_im,Im(rv1),sl2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cimatrix_subv &, const rvector &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv1),rv2);
		addDot(tmp_im,Im(rv1),rv2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rmatrix_subv &, const cimatrix_subv &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,rv1,Re(rv2));
		addDot(tmp_im,rv1,Im(rv2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const rvector_slice & sl1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector_slice &, const cimatrix_subv &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,sl1,Re(rv2));
		addDot(tmp_im,sl1,Im(rv2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp, const rvector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rvector &, const cimatrix_subv &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,rv1,Re(rv2));
		addDot(tmp_im,rv1,Im(rv2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

// complex
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cimatrix_subv &, const cmatrix_subv &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cimatrix_subv &, const cvector_slice &)"));
#endif	
		addDot(dp,rv1,sl2);
	}

	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cimatrix_subv &, const cvector &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cmatrix_subv &, const cimatrix_subv &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

	void accumulate(cidotprecision &dp, const cvector_slice & sl1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector_slice &, const cimatrix_subv &)"));
#endif	
		addDot(dp,sl1,rv2);
	}

	 void accumulate(cidotprecision &dp, const cvector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cvector &, const cimatrix_subv &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

// interval
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cimatrix_subv &, const imatrix_subv &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv1),rv2);
		addDot(tmp_im,Im(rv1),rv2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(sl2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cimatrix_subv &, const ivector_slice &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv1),sl2);
		addDot(tmp_im,Im(rv1),sl2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cimatrix_subv &, const ivector &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv1),rv2);
		addDot(tmp_im,Im(rv1),rv2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const imatrix_subv &, const cimatrix_subv &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,rv1,Re(rv2));
		addDot(tmp_im,rv1,Im(rv2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp, const ivector_slice & sl1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(sl1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector_slice &, const cimatrix_subv &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,sl1,Re(rv2));
		addDot(tmp_im,sl1,Im(rv2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	 void accumulate(cidotprecision &dp, const ivector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const ivector &, const cimatrix_subv &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,rv1,Re(rv2));
		addDot(tmp_im,rv1,Im(rv2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rmatrix_subv &, const civector &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,rv1,Re(rv2));
		addDot(tmp_im,rv1,Im(rv2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const civector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector &, const rmatrix_subv &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv1),rv2);
		addDot(tmp_im,Im(rv1),rv2);
		dp += cidotprecision(tmp_re,tmp_im);
	}
	
// 	void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const civector_slice &rv2)
// #if(CXSC_INDEX_CHECK)
// 	throw(OP_WITH_WRONG_DIM)
// #else
// 	throw()
// #endif
// 	{ 
// #if(CXSC_INDEX_CHECK)
// 		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const rmatrix_subv &, const civector_slice &)"));
// #endif	
// 		idotprecision tmp_re(0.0);
// 		idotprecision tmp_im(0.0);
// 		tmp_re.set_k(dp.get_k());
// 		tmp_im.set_k(dp.get_k());
// 		addDot(tmp_re,rv1,Re(rv2));
// 		addDot(tmp_im,rv1,Im(rv2));
// 		dp += cidotprecision(tmp_re,tmp_im);
// 	}
// 
// 	void accumulate(cidotprecision &dp, const civector_slice & rv1, const rmatrix_subv &rv2)
// #if(CXSC_INDEX_CHECK)
// 	throw(OP_WITH_WRONG_DIM)
// #else
// 	throw()
// #endif
// 	{ 
// #if(CXSC_INDEX_CHECK)
// 		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice &, const rmatrix_subv &)"));
// #endif	
// 		idotprecision tmp_re(0.0);
// 		idotprecision tmp_im(0.0);
// 		tmp_re.set_k(dp.get_k());
// 		tmp_im.set_k(dp.get_k());
// 		addDot(tmp_re,Re(rv1),rv2);
// 		addDot(tmp_im,Im(rv1),rv2);
// 		dp += cidotprecision(tmp_re,tmp_im);
// 	}

	void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const imatrix_subv &, const civector &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,rv1,Re(rv2));
		addDot(tmp_im,rv1,Im(rv2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const civector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector &, const imatrix_subv &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv1),rv2);
		addDot(tmp_im,Im(rv1),rv2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	
	void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const civector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice &, const rmatrix_subv &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,rv1,Re(rv2));
		addDot(tmp_im,rv1,Im(rv2));
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const civector_slice & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice &, const imatrix_subv &)"));
#endif	
		idotprecision tmp_re(0.0);
		idotprecision tmp_im(0.0);
		tmp_re.set_k(dp.get_k());
		tmp_im.set_k(dp.get_k());
		addDot(tmp_re,Re(rv1),rv2);
		addDot(tmp_im,Im(rv1),rv2);
		dp += cidotprecision(tmp_re,tmp_im);
	}

	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cmatrix_subv &, const civector &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

	void accumulate(cidotprecision &dp, const civector & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector &, const cmatrix_subv &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const civector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cmatrix_subv &, const civector_slice &)"));
#endif	
		addDot(dp,rv1,rv2);
	}

	void accumulate(cidotprecision &dp, const civector_slice & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ 
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const civector_slice &, const cmatrix_subv &)"));
#endif	
		addDot(dp,rv1,rv2);
	}
	
	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const cmatrix_subv &, const imatrix_subv &)"));
#endif
		idotprecision idot(0.0);
		idot.set_k(dp.get_k());

		addDot(idot,Re(rv1),rv2);
		SetRe(dp,Re(dp)+idot);

		idot = 0.0;
		addDot(idot,Im(rv1),rv2);
		SetIm(dp,Im(dp)+idot);
	}

	void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(VecLen(rv1)!=VecLen(rv2)) cxscthrow(OP_WITH_WRONG_DIM("void accumulate(cidotprecision&, const imatrix_subv &, const cmatrix_subv &)"));
#endif
		idotprecision idot(0.0);
		idot.set_k(dp.get_k());

		addDot(idot,rv1,Re(rv2));
		SetRe(dp,Re(dp)+idot);

		idot = 0.0;
		addDot(idot,rv1,Im(rv2));
		SetIm(dp,Im(dp)+idot);
	}

} // namespace cxsc

