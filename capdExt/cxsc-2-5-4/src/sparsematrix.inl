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

/* CVS $Id: sparsematrix.inl,v 1.25 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef CXSC_SPARSEMATRIX_INLINE
#define CXSC_SPARSEMATRIX_INLINE

#ifdef CXSC_SPARSE_ACCU_MULT_SAVE_MEMORY
#define SPMULT_SAVE_MEMORY true
#else
#define SPMULT_SAVE_MEMORY false
#endif

#include <sstream>
#include <cctype>

namespace cxsc {

class sparse_dot;
class sparse_cdot;

template<class T>
inline bool isReal() {
  return false;
}

template<class T>
inline bool isComplex() {
  return false;
}


template<> inline bool isReal<real>() {
  return true;
}

template<> inline bool isReal<sparse_dot>() {
  return true;
}

template<> inline bool isComplex<complex>() {
  return true;
}

template<> inline bool isComplex<sparse_cdot>() {
  return true;
}


template<class TA, class Tx, class Tres, class TDot, class TElement>
inline Tres spsl_mv_mult(const TA& A, const Tx& v) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(A.n!=v.n) cxscthrow(OP_WITH_WRONG_DIM("operator*(const " + nameof(A) + " &, const " + nameof(v) + " &)"));
#endif
  Tres x(A.m,A.get_nnz()+v.get_nnz());

  if(opdotprec == 1 && isReal<TElement>()) {
    std::vector<TElement> dot(A.m, TElement(0.0));

    for(int i=0 ; i<v.get_nnz() ; i++) {
      for(int k=A.p[v.p[i]-v.offset] ; k<A.p[v.p[i]-v.offset+1] ; k++) {
        dot[A.ind[k]] += A.x[k] * v.x[i];
      }
    }

    for(int i=0 ; i<A.m ; i++) {
      if(dot[i] != 0.0) x[i+1] = dot[i];
    }

  } else if(opdotprec == 0 && SPMULT_SAVE_MEMORY) {

    TA At(A.n,A.m);
    At = transp(A);
    TDot dot(opdotprec);

    for(int i=0 ; i<A.m ; i++) {
      dot.reset();

      for(int k=At.p[i] ; k<At.p[i+1] ; k++) {
        for(int j=v.start ; j<=v.end && v.p[j]-v.offset<=At.ind[k] ; j++) {
          if(v.p[j]-v.offset == At.ind[k])
            dot.add_dot(At.x[k], v.x[j]);
        }
      }

      x[i+1] = dot.result();
    }

  } else {

    std::vector<TDot> dot(A.m, TDot(opdotprec));

    for(int i=0 ; i<v.get_nnz() ; i++) {
      for(int k=A.p[v.p[i]-v.offset] ; k<A.p[v.p[i]-v.offset+1] ; k++) {
        dot[A.ind[k]].add_dot(A.x[k], v.x[i]);
      }
    }

    TElement tmp;
    for(int i=0 ; i<A.m ; i++) {
      tmp = dot[i].result();
      if(tmp != 0.0) x[i+1] = tmp;
    }

  }

  return x;
}

template<class TA, class Tx, class Tres, class TDot, class TElement>
inline Tres spsp_mv_mult(const TA& A, const Tx& v) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(A.n!=v.n) cxscthrow(OP_WITH_WRONG_DIM("operator*(const " + nameof(A) + " &, const " + nameof(v) + " &)"));
#endif
  Tres x(A.m,A.get_nnz()+v.get_nnz());

  if(opdotprec == 1 && isReal<TElement>()) {
    std::vector<TElement> dot(A.m, TElement(0.0));

    for(int i=0 ; i<v.get_nnz() ; i++) {
      for(int k=A.p[v.p[i]] ; k<A.p[v.p[i]+1] ; k++) {
        dot[A.ind[k]] += A.x[k] * v.x[i];
      }
    }

    for(int i=0 ; i<A.m ; i++) {
      if(dot[i] != 0.0) x[i+1] = dot[i];
    }

  } else if(opdotprec == 0 && SPMULT_SAVE_MEMORY) {

    TA At(A.n,A.m);
    At = transp(A);
    TDot dot(opdotprec);

    for(int i=0 ; i<A.m ; i++) {
      dot.reset();

      for(int k=At.p[i] ; k<At.p[i+1] ; k++) {
        for(int j=0 ; j<v.get_nnz() && v.p[j]<=At.ind[k] ; j++) {
          if(v.p[j] == At.ind[k])
            dot.add_dot(At.x[k], v.x[j]);
        }
      }

      x[i+1] = dot.result();
    }

  } else {
    std::vector<TDot> dot(A.m, TDot(opdotprec));

    for(int i=0 ; i<v.get_nnz() ; i++) {
      for(int k=A.p[v.p[i]] ; k<A.p[v.p[i]+1] ; k++) {
        dot[A.ind[k]].add_dot(A.x[k], v.x[i]);
      }
    }

    TElement tmp;
    for(int i=0 ; i<A.m ; i++) {
      tmp = dot[i].result();
      if(tmp != 0.0) x[i+1] = tmp;
    }
  }

  return x;
}

template<class TA, class Tx, class Tres, class TDot>
inline Tres spf_mv_mult(const TA& A, const Tx& r) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(A.n!=VecLen(r)) cxscthrow(OP_WITH_WRONG_DIM("operator*(const " + nameof(A) + " &, const " + nameof(r) + " &)"));
#endif

  Tres x(A.m);

  if(opdotprec == 1 && isReal<TDot>()) {
    x = 0.0;

    for(int j=0 ; j<A.n ; j++) {
      for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
        x[A.ind[k]+1] += A.x[k] * r[j+Lb(r)];
      }
    }

  } if(opdotprec == 0 && SPMULT_SAVE_MEMORY) {

    TA At(A.n,A.m);
    At = transp(A);
    TDot dot(opdotprec);

    for(int i=0 ; i<A.m ; i++) {
      dot.reset();

      for(int k=At.p[i] ; k<At.p[i+1] ; k++) {
        dot.add_dot(At.x[k], r[At.ind[k]+Lb(r)]);
      }

      x[i+1] = dot.result();
    }

  } else {
    std::vector<TDot> dot(A.m, TDot(opdotprec));

    for(int j=0 ; j<A.n ; j++) {
      for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
        dot[A.ind[k]].add_dot(A.x[k], r[j+Lb(r)]);
      }
    }

    for(int i=0 ; i<A.m ; i++)
      x[i+1] = dot[i].result();

  }

  return x;
}

template<class TA, class Tx, class Tres, class TDot>
inline Tres fsp_mv_mult(const TA& A, const Tx& v) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=v.n) cxscthrow(OP_WITH_WRONG_DIM("operator*(const " + nameof(A) + " &, const " + nameof(v) + " &)"));
#endif
  Tres x(v.n);
  int lb1 = Lb(A,1);
  int lb2 = Lb(A,2);

  if(opdotprec == 1 && isReal<TDot>()) {
    x = 0.0;

    for(int i=0 ; i<v.n ; i++) {
      for(unsigned int j=0 ; j<v.p.size() ; j++) {
        x[i+1] += A[i+lb1][v.p[j]+lb2] * v.x[j];
      }
    }

  } else {

    TDot dot(opdotprec);

    for(int i=0 ; i<v.n ; i++) {
      dot.reset();
      for(unsigned int j=0 ; j<v.p.size() ; j++) {
        dot.add_dot(A[i+lb1][v.p[j]+lb2], v.x[j]);
      }
      x[i+1] = dot.result();
    }
  }
  return x;
}

template<class TA, class Tx, class Tres, class TDot>
inline Tres fsl_mv_mult(const TA& A, const Tx& v) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=v.n) cxscthrow(OP_WITH_WRONG_DIM("operator*(const " + nameof(A) + " &, const " + nameof(v) + " &)"));
#endif
  Tres x(v.n);
  int lb1 = Lb(A,1);
  int lb2 = Lb(A,2);

  if(opdotprec == 1 && isReal<TDot>()) {
    x = 0.0;
    for(int i=0 ; i<v.n ; i++) {
      for(int j=v.start ; j<=v.end ; j++) {
        x[i+1] += A[i+lb1][v.p[j]-v.offset+lb2] * v.x[j];
      }
    }

  } else {

    TDot dot(opdotprec);
    for(int i=0 ; i<v.n ; i++) {
      dot.reset();
      for(int j=v.start ; j<=v.end ; j++) {
        dot.add_dot(A[i+lb1][v.p[j]-v.offset+lb2], v.x[j]);
      }
      x[i+1] = dot.result();
    }
  }

  return x;
}


//----------------------------------------------------------------------------

template<class Tx, class Ty, class Tz>
inline void fp_mult(const Tx& x, const Ty& y, Tz& z) {
  z = x * y;
}

template<>
inline void fp_mult<complex,complex,complex>(const complex& x, const complex& y, complex& z) {
  Re(z) = Re(x)*Re(y)-Im(x)*Im(y);
  Im(z) = Re(x)*Im(y)+Im(x)*Re(y);
}

template<class Tx, class Ty, class Tz>
inline void fp_multadd(const Tx& x, const Ty& y, Tz& z) {
  z += x * y;
}

template<>
inline void fp_multadd<complex,complex,complex>(const complex& x, const complex& y, complex& z) {
  Re(z) += Re(x)*Re(y)-Im(x)*Im(y);
  Im(z) += Re(x)*Im(y)+Im(x)*Re(y);
}

template<class TA, class TB, class Tres, class TDot, class TElement>
inline Tres spsp_mm_mult(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{  
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator*(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int m = ColLen(A);
  //int n = RowLen(A);
  int o = RowLen(B);

  Tres C(m,o,3*(A.get_nnz()+B.get_nnz()));

  if(opdotprec == 1 && isReal<TElement>()) {
    std::vector<TElement> dot(m,TElement(0.0));
    std::vector<int> w(m,-1);
    int nnz = 0;

    for(int i=0 ; i<o ; i++) {

      for(int k=B.p[i] ; k<B.p[i+1] ; k++) {
        for(int l=A.p[B.ind[k]] ; l<A.p[B.ind[k]+1] ; l++) {
          if(w[A.ind[l]] < i) {
            w[A.ind[l]] = i;
            C.ind.push_back(A.ind[l]);
            dot[A.ind[l]] = A.x[l] * B.x[k];
            nnz++;
          } else {
            dot[A.ind[l]] += A.x[l] * B.x[k];
          }
        }
      }

      sort(C.ind.begin()+C.p[i], C.ind.begin()+nnz);

      for(int j=C.p[i] ; j<nnz ; j++) {
         C.x.push_back(dot[C.ind[j]]);
      }

      C.p[i+1] = nnz;
    }

  } else if(opdotprec == 1 && isComplex<TElement>()) {

    std::vector<TElement> dot(m,TElement(0.0));
    std::vector<int> w(m,-1);
    int nnz = 0;

    for(int i=0 ; i<o ; i++) {

      for(int k=B.p[i] ; k<B.p[i+1] ; k++) {
        for(int l=A.p[B.ind[k]] ; l<A.p[B.ind[k]+1] ; l++) {
          if(w[A.ind[l]] < i) {
            w[A.ind[l]] = i;
            C.ind.push_back(A.ind[l]);
            //dot[A.ind[l]] = A.x[l] * B.x[k];
            fp_mult(A.x[l],B.x[k],dot[A.ind[l]]);
            nnz++;
          } else {
            //dot[A.ind[l]] += A.x[l] * B.x[k];
            fp_multadd(A.x[l],B.x[k],dot[A.ind[l]]);
          }
        }
      }

      sort(C.ind.begin()+C.p[i], C.ind.begin()+nnz);

      for(int j=C.p[i] ; j<nnz ; j++) {
         C.x.push_back(dot[C.ind[j]]);
      }

      C.p[i+1] = nnz;
    }

  } else if(opdotprec != 0 || !SPMULT_SAVE_MEMORY) {

    std::vector<TDot> dot(m,TDot(opdotprec));
    std::vector<int> w(m,-1);
    int nnz = 0;

    for(int i=0 ; i<o ; i++) {

      for(int k=B.p[i] ; k<B.p[i+1] ; k++) {
        for(int l=A.p[B.ind[k]] ; l<A.p[B.ind[k]+1] ; l++) {
          if(w[A.ind[l]] < i) {
            w[A.ind[l]] = i;
            C.ind.push_back(A.ind[l]);
            dot[A.ind[l]].reset();
            dot[A.ind[l]].add_dot(A.x[l],B.x[k]);
            nnz++;
          } else {
            dot[A.ind[l]].add_dot(A.x[l],B.x[k]);
          }
        }
      }

      sort(C.ind.begin()+C.p[i], C.ind.begin()+nnz);

      for(int j=C.p[i] ; j<nnz ; j++) {
         C.x.push_back(dot[C.ind[j]].result());
      }

      C.p[i+1] = nnz;
    }

  } else {
  
    TA At(A.n,A.m);
    At = transp(A);

    TDot dot(opdotprec);

    int nnz = 0;

    for(int i=0 ; i<o ; i++) {
      for(int j=0 ; j<m ; j++) {
        dot.reset();

        for(int k=At.p[j] ; k<At.p[j+1] ; k++) {
          for(int l=B.p[i] ; l<B.p[i+1] && At.ind[k] >= B.ind[l] ; l++) {
            if(At.ind[k] == B.ind[l]) {
              dot.add_dot(At.x[k], B.x[l]);
            }
          }
        }

        TElement result = dot.result();
        if(result != 0.0) {
          C.x.push_back(result);
          C.ind.push_back(j);
          nnz++;
        }
      }

      C.p[i+1] = nnz;
    }


  }

  return C;
}

template<class TA, class TB, class Tres, class TDot>
inline Tres fsp_mm_mult(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator*(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int lb1 = Lb(A,1);
  int lb2 = Lb(A,2);

  int m = ColLen(A);
  //int n = RowLen(A);
  int o = RowLen(B);

  Tres C(m,o);

  C = 0.0;

  if(opdotprec == 1 && isReal<TDot>()) {

    for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
        for(int k=B.p[j] ; k<B.p[j+1] ; k++) {
          C[i+1][j+1] += A[i+lb1][B.ind[k]+lb2] * B.x[k];
        }
      }
    }

  } else {

    TDot dot(opdotprec);

    for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
        dot.reset();
        for(int k=B.p[j] ; k<B.p[j+1] ; k++) {
          dot.add_dot(A[i+lb1][B.ind[k]+lb2], B.x[k]);
        }
        C[i+1][j+1] = dot.result();
      }
    }
  }

  return C;
}

template<class TA, class TB, class Tres, class TDot>
inline Tres spf_mm_mult(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator*(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int lb1 = Lb(B,1);
  int lb2 = Lb(B,2);

  int m = ColLen(A);
  //int n = RowLen(A);
  int o = RowLen(B);

  Tres C(m,o);

  C = 0.0;

  TA At(A.n,A.m);
  At = transp(A);

  if(opdotprec == 1 && isReal<TDot>()) {

    for(int i=0 ; i<o ; i++) {
      for(int j=0 ; j<m ; j++) {
        for(int k=At.p[j] ; k<At.p[j+1] ; k++) {
           C[j+1][i+1] += At.x[k] * B[At.ind[k]+lb1][i+lb2];
        }
      }
    }

  } else {

    TDot dot(opdotprec);

    for(int i=0 ; i<o ; i++) {
      for(int j=0 ; j<m ; j++) {
        dot.reset();
        for(int k=At.p[j] ; k<At.p[j+1] ; k++) {
           dot.add_dot(At.x[k], B[At.ind[k]+lb1][i+lb2]);
        }
        C[j+1][i+1] = dot.result();
      }
    }

  }

  return C;
}

//--------------------------------------------------------------------------

template<class TA, class Ts, class Tres>
inline Tres sp_ms_div(const TA& A, const Ts& s) 
{
  if(s == 0.0) {
    Tres C(A.m,A.n);
    return C;
  } else {
    Tres C(A);
    for(unsigned int i=0 ; i<A.x.size() ; i++) 
      C.x[i] /= s;
    return C;
  }
}

template<class TA, class Ts, class Tres>
inline Tres sp_ms_mult(const TA& A, const Ts& s) {
  if(s == 0.0) {
    Tres C(A.m,A.n);
    return C;
  } else {
    Tres C(A);
    for(unsigned int i=0 ; i<A.x.size() ; i++) 
      C.x[i] *= s;
    return C;
  }
}

template<class Ts, class TA, class Tres>
inline Tres sp_sm_mult(const Ts& s, const TA& A) {
  if(s == 0.0) {
    Tres C(A.m,A.n);
    return C;
  } else {
    Tres C(A);
    for(unsigned int i=0 ; i<A.x.size() ; i++) 
      C.x[i] *= s;
    return C;
  }
}


//--------------------------------------------------------------------------

template<class TA, class TB, class Tres, class TElement>
inline Tres spsp_mm_add(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator+(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int m = ColLen(A);
  int n = RowLen(A);

  Tres C(m, n, A.get_nnz()+B.get_nnz());

  int nnz = 0;

  for(int j=0 ; j<n ; j++) {

    int k = A.p[j];
    int l = B.p[j]; 

    while(k<A.p[j+1] && l<B.p[j+1]) {
      if(A.ind[k] == B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(A.x[k] + B.x[l]);
        k++; l++;
      } else if(A.ind[k] < B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(TElement(A.x[k]));
        k++;
      } else {
        C.ind.push_back(B.ind[l]);
        C.x.push_back(TElement(B.x[l]));
        l++;
      }
      nnz++;
    }

    for( ; k<A.p[j+1] ; k++) {
      C.ind.push_back(A.ind[k]);
      C.x.push_back(TElement(A.x[k]));
      nnz++;
    }

    for( ; l<B.p[j+1] ; l++) {
      C.ind.push_back(B.ind[l]);
      C.x.push_back(TElement(B.x[l]));
      nnz++;
    }

    C.p[j+1] = nnz;

  }

  return C;
}

template<class TA, class TB, class Tres>
inline Tres spf_mm_add(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator+(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  int n = RowLen(A);
  
  Tres C(B);
  SetLb(C,ROW,1);
  SetLb(C,COL,1);

  for(int j=0 ; j<n ; j++) {
    for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
      C[A.ind[k]+1][j+1] += A.x[k];
    }
  }

  return C;
}

template<class TA, class TB, class Tres>
inline Tres fsp_mm_add(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator+(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  int n = RowLen(A);
  
  Tres C(A);
  SetLb(C,ROW,1);
  SetLb(C,COL,1);

  for(int j=0 ; j<n ; j++) {
    for(int k=B.p[j] ; k<B.p[j+1] ; k++) {
      C[B.ind[k]+1][j+1] += B.x[k];
    }
  }

  return C;
}


//--------------------------------------------------------------------------

template<class TA, class TB, class Tres, class TElement>
inline Tres spsp_mm_sub(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator-(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int m = ColLen(A);
  int n = RowLen(A);

  Tres C(m, n, A.get_nnz()+B.get_nnz());

  int nnz = 0;

  for(int j=0 ; j<n ; j++) {

    int k = A.p[j];
    int l = B.p[j]; 

    while(k<A.p[j+1] && l<B.p[j+1]) {
      if(A.ind[k] == B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(A.x[k] - B.x[l]);
        k++; l++;
      } else if(A.ind[k] < B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(TElement(A.x[k]));
        k++;
      } else {
        C.ind.push_back(B.ind[l]);
        C.x.push_back(TElement(-B.x[l]));
        l++;
      }
      nnz++;
    }

    for( ; k<A.p[j+1] ; k++) {
      C.ind.push_back(A.ind[k]);
      C.x.push_back(TElement(A.x[k]));
      nnz++;
    }

    for( ; l<B.p[j+1] ; l++) {
      C.ind.push_back(B.ind[l]);
      C.x.push_back(TElement(-B.x[l]));
      nnz++;
    }

    C.p[j+1] = nnz;

  }

  return C;
}

template<class TA, class TB, class Tres>
inline Tres spf_mm_sub(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator-(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  int n = RowLen(A);
  
  Tres C(-B);
  SetLb(C,ROW,1);
  SetLb(C,COL,1);

  for(int j=0 ; j<n ; j++) {
    for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
      C[A.ind[k]+1][j+1] += A.x[k];
    }
  }

  return C;
}

template<class TA, class TB, class Tres>
inline Tres fsp_mm_sub(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator-(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  int n = RowLen(A);
  
  Tres C(A);
  SetLb(C,ROW,1);
  SetLb(C,COL,1);

  for(int j=0 ; j<n ; j++) {
    for(int k=B.p[j] ; k<B.p[j+1] ; k++) {
      C[B.ind[k]+1][j+1] -= B.x[k];
    }
  }

  return C;
}

//--------------------------------------------------------------------------

template<class TA, class TB, class Tres, class TElement>
inline Tres spsp_mm_hull(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator|(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int m = ColLen(A);
  int n = RowLen(A);

  Tres C(m, n, A.get_nnz()+B.get_nnz());

  int nnz = 0;

  for(int j=0 ; j<n ; j++) {

    int k = A.p[j];
    int l = B.p[j]; 

    while(k<A.p[j+1] && l<B.p[j+1]) {
      if(A.ind[k] == B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(A.x[k] | B.x[l]);
        k++; l++;
      } else if(A.ind[k] < B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(TElement(A.x[k]));
        k++;
      } else {
        C.ind.push_back(B.ind[l]);
        C.x.push_back(TElement(B.x[l]));
        l++;
      }
      nnz++;
    }

    for( ; k<A.p[j+1] ; k++) {
      C.ind.push_back(A.ind[k]);
      C.x.push_back(TElement(A.x[k]));
      nnz++;
    }

    for( ; l<B.p[j+1] ; l++) {
      C.ind.push_back(B.ind[l]);
      C.x.push_back(TElement(B.x[l]));
      nnz++;
    }

    C.p[j+1] = nnz;

  }

  return C;
}

template<class TA, class TB, class Tres>
inline Tres spf_mm_hull(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator|(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  int n = RowLen(A);
  
  Tres C(B);
  SetLb(C,ROW,1);
  SetLb(C,COL,1);

  for(int j=0 ; j<n ; j++) {
    for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
      C[A.ind[k]+1][j+1] |= A.x[k];
    }
  }

  return C;
}

template<class TA, class TB, class Tres>
inline Tres fsp_mm_hull(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator|(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  int n = RowLen(A);
  
  Tres C(A);
  SetLb(C,ROW,1);
  SetLb(C,COL,1);

  for(int j=0 ; j<n ; j++) {
    for(int k=B.p[j] ; k<B.p[j+1] ; k++) {
      C[B.ind[k]+1][j+1] |= B.x[k];
    }
  }

  return C;
}

//--------------------------------------------------------------------------

template<class TA, class TB, class Tres, class TElement>
inline Tres spsp_mm_intersect(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator&(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int m = ColLen(A);
  int n = RowLen(A);

  Tres C(m, n, A.get_nnz()+B.get_nnz());

  int nnz = 0;

  for(int j=0 ; j<n ; j++) {

    int k = A.p[j];
    int l = B.p[j]; 

    while(k<A.p[j+1] && l<B.p[j+1]) {
      if(A.ind[k] == B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(A.x[k] & B.x[l]);
        k++; l++;
      } else if(A.ind[k] < B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(TElement(A.x[k]));
        k++;
      } else {
        C.ind.push_back(B.ind[l]);
        C.x.push_back(TElement(B.x[l]));
        l++;
      }
      nnz++;
    }

    for( ; k<A.p[j+1] ; k++) {
      C.ind.push_back(A.ind[k]);
      C.x.push_back(TElement(A.x[k]));
      nnz++;
    }

    for( ; l<B.p[j+1] ; l++) {
      C.ind.push_back(B.ind[l]);
      C.x.push_back(TElement(B.x[l]));
      nnz++;
    }

    C.p[j+1] = nnz;

  }


  return C;
}

template<class TA, class TB, class Tres>
inline Tres spf_mm_intersect(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator&(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  int n = RowLen(A);
  
  Tres C(B);
  SetLb(C,ROW,1);
  SetLb(C,COL,1);

  for(int j=0 ; j<n ; j++) {
    for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
      C[A.ind[k]+1][j+1] &= A.x[k];
    }
  }

  return C;
}

template<class TA, class TB, class Tres>
inline Tres fsp_mm_intersect(const TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator&(const " + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  int n = RowLen(A);
  
  Tres C(A);
  SetLb(C,ROW,1);
  SetLb(C,COL,1);

  for(int j=0 ; j<n ; j++) {
    for(int k=B.p[j] ; k<B.p[j+1] ; k++) {
      C[B.ind[k]+1][j+1] &= B.x[k];
    }
  }

  return C;
}

//--------------------------------------------------------------------------

template<class TA, class TB>
inline bool spsp_mm_comp(const TA& A, const TB& B) {
  if(A.m != B.m  ||  A.n != B.n) return false;

  for(int j=0 ; j<A.n ; j++) {
    int k=A.p[j], l=B.p[j];
    while(k < A.p[j+1]  &&  l < B.p[j+1]) {
      if(A.ind[k] < B.ind[l]) {
        if(A.x[k] != 0.0) return false;
        k++;
      } else if(A.ind[k] == B.ind[l]) {
        if(A.x[k] != B.x[l]) return false;
        k++; l++;
      } else if(A.ind[k] > B.ind[l]) {
        if(B.x[l] != 0.0) return false;
        l++;
      }
    }

    for( ; k<A.p[j+1] ; k++)
      if(A.x[k] != 0.0) return false;

    for( ; l<B.p[j+1] ; l++)
      if(B.x[l] != 0.0) return false;

  }

  return true;
}

template<class TA, class TB>
inline bool spf_mm_comp(const TA& A, const TB& B) {
  if(A.m != ColLen(B)  ||  A.n != RowLen(B)) return false;
 
  int lb1 = Lb(B,1);
  int lb2 = Lb(B,2);

  for(int j=0 ; j<A.n ; j++) {
    int k=A.p[j];
    for(int i=0 ; i<A.m ; i++) {
      if(A.ind[k] == i) {
        if(A.x[k] != B[i+lb1][j+lb2]) return false;
        if(k<A.p[j+1]-1) k++;
      } else {
        if(B[i+lb1][j+lb2] != 0.0) return false;
      }
    }
  }

  return true;
}

template<class TA, class TB>
inline bool fsp_mm_comp(const TA& A, const TB& B) {
  if(ColLen(A) != B.m  ||  RowLen(A) != B.n) return false;
 
  int lb1 = Lb(A,1);
  int lb2 = Lb(A,2);

  for(int j=0 ; j<B.n ; j++) {
    int k=B.p[j];
    for(int i=0 ; i<B.m ; i++) {
      if(B.ind[k] == i) {
        if(B.x[k] != A[i+lb1][j+lb2]) return false;
        if(k<B.p[j+1]-1) k++;
      } else {
        if(A[i+lb1][j+lb2] != 0.0) return false;
      }
    }
  }

  return true;
}

//--------------------------------------------------------------------------

template<class TA, class TB, class TType>
inline bool spsp_mm_less(const TA& A, const TB& B) {
  if(A.m != B.m  ||  A.n != B.n) return false;

  for(int j=0 ; j<A.n ; j++) {
    int k=A.p[j], l=B.p[j];
    while(k < A.p[j+1]  &&  l < B.p[j+1]) {
      if(A.ind[k] < B.ind[l]) {
        if(!(A.x[k] < TType(0.0))) return false;
        k++;
      } else if(A.ind[k] == B.ind[l]) {
        if(A.x[k] < B.x[l]) return false;
        k++; l++;
      } else if(A.ind[k] > B.ind[l]) {
        if(!(TType(0.0) < B.x[l])) return false;
        l++;
      }
    }

    for( ; k<A.p[j+1] ; k++)
      if(!(A.x[k] < TType(0.0))) return false;

    for( ; l<B.p[j+1] ; l++)
      if(!(TType(0.0) < B.x[l])) return false;

  }

  return true;
}

template<class TA, class TB, class TType>
inline bool spf_mm_less(const TA& A, const TB& B) {
  if(A.m != ColLen(B)  ||  A.n != RowLen(B)) return false;
 
  int lb1 = Lb(B,1);
  int lb2 = Lb(B,2);

  for(int j=0 ; j<A.n ; j++) {
    int k=A.p[j];
    for(int i=0 ; i<A.m ; i++) {
      if(A.ind[k] == i) {
        if(!(A.x[k] < B[i+lb1][j+lb2])) return false;
        if(k<A.p[j+1]-1) k++;
      } else {
        if(!(TType(0.0) < B[i+lb1][j+lb2])) return false;
      }
    }
  }

  return true;
}

template<class TA, class TB, class TType>
inline bool fsp_mm_less(const TA& A, const TB& B) {
  if(ColLen(A) != B.m  ||  RowLen(A) != B.n) return false;
 
  int lb1 = Lb(A,1);
  int lb2 = Lb(A,2);

  for(int j=0 ; j<B.n ; j++) {
    int k=B.p[j];
    for(int i=0 ; i<B.m ; i++) {
      if(B.ind[k] == i) {
        if(!(A[i+lb1][j+lb2] < B.x[k])) return false;
        if(k<B.p[j+1]-1) k++;
      } else {
        if(A[i+lb1][j+lb2] < TType(0.0)) return false;
      }
    }
  }

  return true;
}

//--------------------------------------------------------------------------

template<class TA, class TB, class TType>
inline bool spsp_mm_leq(const TA& A, const TB& B) {
  if(A.m != B.m  ||  A.n != B.n) return false;

  for(int j=0 ; j<A.n ; j++) {
    int k=A.p[j], l=B.p[j];
    while(k < A.p[j+1]  &&  l < B.p[j+1]) {
      if(A.ind[k] < B.ind[l]) {
        if(!(A.x[k] <= TType(0.0))) return false;
        k++;
      } else if(A.ind[k] == B.ind[l]) {
        if(A.x[k] <= B.x[l]) return false;
        k++; l++;
      } else if(A.ind[k] > B.ind[l]) {
        if(!(TType(0.0) <= B.x[l])) return false;
        l++;
      }
    }

    for( ; k<A.p[j+1] ; k++)
      if(!(A.x[k] <= TType(0.0))) return false;

    for( ; l<B.p[j+1] ; l++)
      if(!(TType(0.0) <= B.x[l])) return false;

  }

  return true;
}

template<class TA, class TB, class TType>
inline bool spf_mm_leq(const TA& A, const TB& B) {
  if(A.m != ColLen(B)  ||  A.n != RowLen(B)) return false;
 
  int lb1 = Lb(B,1);
  int lb2 = Lb(B,2);

  for(int j=0 ; j<A.n ; j++) {
    int k=A.p[j];
    for(int i=0 ; i<A.m ; i++) {
      if(A.ind[k] == i) {
        if(!(A.x[k] <= B[i+lb1][j+lb2])) return false;
        if(k<A.p[j+1]-1) k++;
      } else {
        if(!(TType(0.0) <= B[i+lb1][j+lb2])) return false;
      }
    }
  }

  return true;
}

template<class TA, class TB, class TType>
inline bool fsp_mm_leq(const TA& A, const TB& B) {
  if(ColLen(A) != B.m  ||  RowLen(A) != B.n) return false;
 
  int lb1 = Lb(A,1);
  int lb2 = Lb(A,2);

  for(int j=0 ; j<B.n ; j++) {
    int k=B.p[j];
    for(int i=0 ; i<B.m ; i++) {
      if(B.ind[k] == i) {
        if(!(A[i+lb1][j+lb2] <= B.x[k])) return false;
        if(k<B.p[j+1]-1) k++;
      } else {
        if(A[i+lb1][j+lb2] <= TType(0.0)) return false;
      }
    }
  }

  return true;
}


//--------------------------------------------------------------------------

template<class TA, class TB, class TType>
inline bool spsp_mm_greater(const TA& A, const TB& B) {
  if(A.m != B.m  ||  A.n != B.n) return false;

  for(int j=0 ; j<A.n ; j++) {
    int k=A.p[j], l=B.p[j];
    while(k < A.p[j+1]  &&  l < B.p[j+1]) {
      if(A.ind[k] < B.ind[l]) {
        if(!(A.x[k] > TType(0.0))) return false;
        k++;
      } else if(A.ind[k] == B.ind[l]) {
        if(A.x[k] > B.x[l]) return false;
        k++; l++;
      } else if(A.ind[k] > B.ind[l]) {
        if(!(TType(0.0) > B.x[l])) return false;
        l++;
      }
    }

    for( ; k<A.p[j+1] ; k++)
      if(!(A.x[k] > TType(0.0))) return false;

    for( ; l<B.p[j+1] ; l++)
      if(!(TType(0.0) > B.x[l])) return false;

  }

  return true;
}

template<class TA, class TB, class TType>
inline bool spf_mm_greater(const TA& A, const TB& B) {
  if(A.m != ColLen(B)  ||  A.n != RowLen(B)) return false;
 
  int lb1 = Lb(B,1);
  int lb2 = Lb(B,2);

  for(int j=0 ; j<A.n ; j++) {
    int k=A.p[j];
    for(int i=0 ; i<A.m ; i++) {
      if(A.ind[k] == i) {
        if(!(A.x[k] > B[i+lb1][j+lb2])) return false;
        if(k<A.p[j+1]-1) k++;
      } else {
        if(!(TType(0.0) > B[i+lb1][j+lb2])) return false;
      }
    }
  }

  return true;
}

template<class TA, class TB, class TType>
inline bool fsp_mm_greater(const TA& A, const TB& B) {
  if(ColLen(A) != B.m  ||  RowLen(A) != B.n) return false;
 
  int lb1 = Lb(A,1);
  int lb2 = Lb(A,2);

  for(int j=0 ; j<B.n ; j++) {
    int k=B.p[j];
    for(int i=0 ; i<B.m ; i++) {
      if(B.ind[k] == i) {
        if(!(A[i+lb1][j+lb2] > B.x[k])) return false;
        if(k<B.p[j+1]-1) k++;
      } else {
        if(A[i+lb1][j+lb2] > TType(0.0)) return false;
      }
    }
  }

  return true;
}

//--------------------------------------------------------------------------

template<class TA, class TB, class TType>
inline bool spsp_mm_geq(const TA& A, const TB& B) {
  if(A.m != B.m  ||  A.n != B.n) return false;

  for(int j=0 ; j<A.n ; j++) {
    int k=A.p[j], l=B.p[j];
    while(k < A.p[j+1]  &&  l < B.p[j+1]) {
      if(A.ind[k] < B.ind[l]) {
        if(!(A.x[k] >= TType(0.0))) return false;
        k++;
      } else if(A.ind[k] == B.ind[l]) {
        if(A.x[k] >= B.x[l]) return false;
        k++; l++;
      } else if(A.ind[k] > B.ind[l]) {
        if(!(TType(0.0) >= B.x[l])) return false;
        l++;
      }
    }

    for( ; k<A.p[j+1] ; k++)
      if(!(A.x[k] >= TType(0.0))) return false;

    for( ; l<B.p[j+1] ; l++)
      if(!(TType(0.0) >= B.x[l])) return false;

  }

  return true;
}

template<class TA, class TB, class TType>
inline bool spf_mm_geq(const TA& A, const TB& B) {
  if(A.m != ColLen(B)  ||  A.n != RowLen(B)) return false;
 
  int lb1 = Lb(B,1);
  int lb2 = Lb(B,2);

  for(int j=0 ; j<A.n ; j++) {
    int k=A.p[j];
    for(int i=0 ; i<A.m ; i++) {
      if(A.ind[k] == i) {
        if(!(A.x[k] >= B[i+lb1][j+lb2])) return false;
        if(k<A.p[j+1]-1) k++;
      } else {
        if(!(TType(0.0) >= B[i+lb1][j+lb2])) return false;
      }
    }
  }

  return true;
}

template<class TA, class TB, class TType>
inline bool fsp_mm_geq(const TA& A, const TB& B) {
  if(ColLen(A) != B.m  ||  RowLen(A) != B.n) return false;
 
  int lb1 = Lb(A,1);
  int lb2 = Lb(A,2);

  for(int j=0 ; j<B.n ; j++) {
    int k=B.p[j];
    for(int i=0 ; i<B.m ; i++) {
      if(B.ind[k] == i) {
        if(!(A[i+lb1][j+lb2] >= B.x[k])) return false;
        if(k<B.p[j+1]-1) k++;
      } else {
        if(A[i+lb1][j+lb2] >= TType(0.0)) return false;
      }
    }
  }

  return true;
}

//--------------------------------------------------------------------------

template<class TA, class Tres>
inline Tres sp_m_negative(const TA& M) {
  Tres A(M);
  for(unsigned int i=0 ; i<A.x.size() ; i++)
    A.x[i] = -A.x[i];
  return A;
}

//--------------------------------------------------------------------------



template<class TA, class TType>
inline std::ostream& sp_m_output(std::ostream& os, const TA& A) {
  if(ioflags.isset(IOFlags::fullinout)) {

    for(int k=0 ; k<A.m ; k++) {
 
      for(unsigned int i=0 ; i<A.p.size()-1 ; i++) {
 
        bool found = false;
 
        for(int j=A.p[i] ; j<A.p[i+1] && !found ; j++) {
          if(A.ind[j] == k) {
            os << A.x[j] << " ";
            found = true;
          }
        }

        if(!found) os << (TType)0.0 << " ";

      } 

      os << std::endl;
    }

  } else if(ioflags.isset(IOFlags::sparseinout)) {

    os << A.m << " " << A.n << std::endl;
    os << A.get_nnz() << std::endl;

    for(unsigned int i=0 ; i<A.p.size() ; i++) {
      os << A.p[i] << " ";
    }
    os << std::endl;

    for(unsigned int i=0 ; i<A.ind.size() ; i++) {
      os << A.ind[i] << " ";
    }
    os << std::endl;

    for(unsigned int i=0 ; i<A.x.size() ; i++) {
      os << A.x[i] << " ";
    }
    os << std::endl;

  } else {

    //MatrixMarket  
    os << "%%MatrixMarket matrix coordinate " << ElementName(A) <<" general" << std::endl;
    os << "%Generated by C-XSC" << std::endl;
    os << A.m << " " << A.n << " " << A.get_nnz() << std::endl;
    for(unsigned int i=0 ; i<A.p.size()-1 ; i++) {
      for(int k=A.p[i] ; k<A.p[i+1] ; k++) {
        os << A.ind[k]+1 << " " << i+1 << " " 
           << MatrixMarketElement(A.x[k]) << std::endl;
      }
    }

  }

  return os;
}


inline void mm_add_element(std::vector<triplet_store<real> >& v, string& tmp, MM_FIELD& field, MM_SYMMETRY& symmetry) {
  int i,j;
  std::stringstream ss(tmp);
  real val;

  ss >> i >> j;
  i--; j--;

  if(symmetry == mm_general) {

      if(field != mm_pattern) {
        ss >> val;
        v.push_back(triplet_store<real>(i,j,val));
      } else {
        v.push_back(triplet_store<real>(i,j,1.0));
      }

  } else if(symmetry == mm_symmetric || symmetry == mm_hermitian) {

      if(field != mm_pattern) {
        ss >> val;
        v.push_back(triplet_store<real>(i,j,val));
        if(i!=j) v.push_back(triplet_store<real>(j,i,val));
      } else {
        v.push_back(triplet_store<real>(i,j,1.0));
        if(i!=j) v.push_back(triplet_store<real>(j,i,1.0));
      }

  } else if(symmetry == mm_skew_symmetric) {

      if(field != mm_pattern) {
        ss >> val;
        v.push_back(triplet_store<real>(i,j,val));
        if(i!=j) v.push_back(triplet_store<real>(j,i,-val));
      } else {
        v.push_back(triplet_store<real>(i,j,1.0));
        if(i!=j) v.push_back(triplet_store<real>(j,i,-1.0));
      }

  }

}

inline void mm_add_element(std::vector<triplet_store<complex> >& v, string& tmp, MM_FIELD& field, MM_SYMMETRY& symmetry) {
  int i,j;
  std::stringstream ss(tmp);
  real val_re, val_im;

  ss >> i >> j;
  i--; j--;

  if(symmetry == mm_general) {

      if(field == mm_real || field == mm_integer || field == mm_interval) {
        ss >> val_re;
        v.push_back(triplet_store<complex>(i,j,complex(val_re)));
      } else if(field == mm_complex || field == mm_cinterval) {
        ss >> val_re >> val_im;
        v.push_back(triplet_store<complex>(i,j,complex(val_re,val_im)));
      } else if(field == mm_pattern) {
        v.push_back(triplet_store<complex>(i,j,complex(1.0)));
      }

  } else if(symmetry == mm_symmetric) {

      if(field == mm_real || field == mm_integer || field == mm_interval) {
        ss >> val_re;
        v.push_back(triplet_store<complex>(i,j,complex(val_re)));
        if(i!=j) v.push_back(triplet_store<complex>(j,i,complex(val_re)));
      } else if(field == mm_complex || field == mm_cinterval) {
        ss >> val_re >> val_im;
        v.push_back(triplet_store<complex>(i,j,complex(val_re,val_im)));
        if(i!=j) v.push_back(triplet_store<complex>(j,i,complex(val_re,val_im)));
      } else if(field == mm_pattern) {
        v.push_back(triplet_store<complex>(i,j,complex(1.0)));
        if(i!=j) v.push_back(triplet_store<complex>(j,i,complex(1.0)));
      }

  } else if(symmetry == mm_skew_symmetric) {

      if(field == mm_real || field == mm_integer || field == mm_interval) {
        ss >> val_re;
        v.push_back(triplet_store<complex>(i,j,complex(val_re)));
        if(i!=j) v.push_back(triplet_store<complex>(j,i,complex(-val_re)));
      } else if(field == mm_complex || field == mm_cinterval) {
        ss >> val_re >> val_im;
        v.push_back(triplet_store<complex>(i,j,complex(val_re,val_im)));
        if(i!=j) v.push_back(triplet_store<complex>(j,i,-complex(val_re,val_im)));
      } else if(field == mm_pattern) {
        v.push_back(triplet_store<complex>(i,j,complex(1.0)));
        if(i!=j) v.push_back(triplet_store<complex>(j,i,complex(-1.0)));
      }

  } else if(symmetry == mm_hermitian) {

      if(field == mm_real || field == mm_integer || field == mm_interval) {
        ss >> val_re;
        v.push_back(triplet_store<complex>(i,j,complex(val_re)));
        if(i!=j) v.push_back(triplet_store<complex>(j,i,complex(val_re)));
      } else if(field == mm_complex || field == mm_cinterval) {
        ss >> val_re >> val_im;
        v.push_back(triplet_store<complex>(i,j,complex(val_re,val_im)));
        if(i!=j) v.push_back(triplet_store<complex>(j,i,complex(val_re,-val_im)));
      } else if(field == mm_pattern) {
        v.push_back(triplet_store<complex>(i,j,complex(1.0)));
        if(i!=j) v.push_back(triplet_store<complex>(j,i,complex(1.0)));
      }

  }

}


inline void mm_add_element(std::vector<triplet_store<interval> >& v, string& tmp, MM_FIELD& field, MM_SYMMETRY& symmetry) {
  int i,j;
  std::stringstream ss(tmp);
  real val_inf, val_sup;

  ss >> i >> j;
  i--; j--;

  if(symmetry == mm_general) {

      if(field == mm_interval || field == mm_cinterval) {
        ss >> val_inf >> val_sup;
        v.push_back(triplet_store<interval>(i,j,interval(val_inf,val_sup))); 
      } else if(field != mm_pattern) {
        ss >> val_inf;
        v.push_back(triplet_store<interval>(i,j,interval(val_inf)));
      } else {
        v.push_back(triplet_store<interval>(i,j,interval(1.0)));
      }

  } else if(symmetry == mm_symmetric || symmetry == mm_hermitian) {

      if(field == mm_interval || field == mm_cinterval) {
        ss >> val_inf >> val_sup;
        v.push_back(triplet_store<interval>(i,j,interval(val_inf,val_sup))); 
        if(i!=j) v.push_back(triplet_store<interval>(j,i,interval(val_inf,val_sup))); 
      } else if(field != mm_pattern) {
        ss >> val_inf;
        v.push_back(triplet_store<interval>(i,j,interval(val_inf)));
        if(i!=j) v.push_back(triplet_store<interval>(j,i,interval(val_inf)));
      } else {
        v.push_back(triplet_store<interval>(i,j,interval(1.0)));
        if(i!=j) v.push_back(triplet_store<interval>(i,j,interval(1.0)));
      }

  } else if(symmetry == mm_skew_symmetric) {

      if(field == mm_interval || field == mm_cinterval) {
        ss >> val_inf >> val_sup;
        v.push_back(triplet_store<interval>(i,j,interval(val_inf,val_sup))); 
        if(i!=j) v.push_back(triplet_store<interval>(j,i,-interval(val_inf,val_sup))); 
      } else if(field != mm_pattern) {
        ss >> val_inf;
        v.push_back(triplet_store<interval>(i,j,interval(val_inf)));
        if(i!=j) v.push_back(triplet_store<interval>(j,i,interval(-val_inf)));
      } else {
        v.push_back(triplet_store<interval>(i,j,interval(1.0)));
        if(i!=j) v.push_back(triplet_store<interval>(i,j,interval(-1.0)));
      }

  }
}


inline void mm_add_element(std::vector<triplet_store<cinterval> >& v, string& tmp, MM_FIELD& field, MM_SYMMETRY& symmetry) {
  int i,j;
  std::stringstream ss(tmp);
  real val_inf_re, val_sup_re, val_inf_im, val_sup_im;

  ss >> i >> j;
  i--; j--;

  if(symmetry == mm_general) {

      if(field == mm_real || field == mm_integer) {
        ss >> val_inf_re;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(val_inf_re)));
      } else if(field == mm_complex) {
        ss >> val_inf_re >> val_inf_im;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(complex(val_inf_re,val_inf_im))));
      } else if(field == mm_interval) {
        ss >> val_inf_re >> val_sup_re;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(interval(val_inf_re,val_sup_re))));
      } else if(field == mm_cinterval) {
        ss >> val_inf_re >> val_sup_re >> val_inf_im >> val_sup_im;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(interval(val_inf_re,val_sup_re),interval(val_inf_im,val_sup_im))));
      } else if(field == mm_pattern) {
        v.push_back(triplet_store<cinterval>(i,j,cinterval(1.0)));
      }

  } else if(symmetry == mm_symmetric) {

      if(field == mm_real || field == mm_integer) {
        ss >> val_inf_re;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(val_inf_re)));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(val_inf_re)));
      } else if(field == mm_complex) {
        ss >> val_inf_re >> val_inf_im;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(complex(val_inf_re,val_inf_im))));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(complex(val_inf_re,val_inf_im))));
      } else if(field == mm_interval) {
        ss >> val_inf_re >> val_sup_re;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(interval(val_inf_re,val_sup_re))));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(interval(val_inf_re,val_sup_re))));
      } else if(field == mm_cinterval) {
        ss >> val_inf_re >> val_sup_re >> val_inf_im >> val_sup_im;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(interval(val_inf_re,val_sup_re),interval(val_inf_im,val_sup_im))));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(interval(val_inf_re,val_sup_re),interval(val_inf_im,val_sup_im))));
      } else if(field == mm_pattern) {
        v.push_back(triplet_store<cinterval>(i,j,cinterval(1.0)));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(1.0)));
      }

  } else if(symmetry == mm_skew_symmetric) {

      if(field == mm_real || field == mm_integer) {
        ss >> val_inf_re;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(val_inf_re)));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(-val_inf_re)));
      } else if(field == mm_complex) {
        ss >> val_inf_re >> val_inf_im;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(complex(val_inf_re,val_inf_im))));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(-complex(val_inf_re,val_inf_im))));
      } else if(field == mm_interval) {
        ss >> val_inf_re >> val_sup_re;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(interval(val_inf_re,val_sup_re))));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(-interval(val_inf_re,val_sup_re))));
      } else if(field == mm_cinterval) {
        ss >> val_inf_re >> val_sup_re >> val_inf_im >> val_sup_im;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(interval(val_inf_re,val_sup_re),interval(val_inf_im,val_sup_im))));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,-cinterval(interval(val_inf_re,val_sup_re),interval(val_inf_im,val_sup_im))));
      } else if(field == mm_pattern) {
        v.push_back(triplet_store<cinterval>(i,j,cinterval(1.0)));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(-1.0)));
      }

  } else if(symmetry == mm_hermitian) {

      if(field == mm_real || field == mm_integer) {
        ss >> val_inf_re;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(val_inf_re)));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(val_inf_re)));
      } else if(field == mm_complex) {
        ss >> val_inf_re >> val_inf_im;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(complex(val_inf_re,val_inf_im))));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(complex(val_inf_re,-val_inf_im))));
      } else if(field == mm_interval) {
        ss >> val_inf_re >> val_sup_re;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(interval(val_inf_re,val_sup_re))));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(interval(val_inf_re,val_sup_re))));
      } else if(field == mm_cinterval) {
        ss >> val_inf_re >> val_sup_re >> val_inf_im >> val_sup_im;
        v.push_back(triplet_store<cinterval>(i,j,cinterval(interval(val_inf_re,val_sup_re),interval(-val_inf_im,-val_sup_im))));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(interval(val_inf_re,val_sup_re),interval(-val_inf_im,-val_sup_im))));
      } else if(field == mm_pattern) {
        v.push_back(triplet_store<cinterval>(i,j,cinterval(1.0)));
        if(i!=j) v.push_back(triplet_store<cinterval>(j,i,cinterval(1.0)));
      }

  }
}




template<class TA, class TType>
inline std::istream& sp_m_input(std::istream& is, TA& A) {
  if(ioflags.isset(IOFlags::fullinout)) {

    TA Atmp(A.n,A.m);
    TType tmp;

    for(int i=0 ; i<A.m ; i++) {
      for(int j=0 ; j<A.n ; j++) {
        is >> tmp;
        if(tmp != 0.0) {
          Atmp.x.push_back(tmp);
          Atmp.ind.push_back(j);
          Atmp.p[i+1]++;
        }
      }
    }

    A = transp(Atmp);

  } else if(ioflags.isset(IOFlags::sparseinout)) {

    A.p.clear();
    A.ind.clear();
    A.x.clear();
    
    is >> A.m >> A.n;
    int nnz; 
    is >> nnz;

    int itmp;
    TType xtmp;

    for(int i=0 ; i<A.n+1 ; i++) {
      is >> itmp;
      A.p.push_back(itmp);
    }

    for(int i=0 ; i<nnz ; i++) {
      is >> itmp;
      A.ind.push_back(itmp);
    }

    for(int i=0 ; i<nnz ; i++) {
      is >> xtmp;
      A.x.push_back(xtmp);
    }

  } else {

    //MatrixMarket
    std::string header;
    std::getline(is,header);
    std::stringstream ss(header);
    std::string tmp, s_format, s_field, s_symmetry;
    //MM_FORMAT format;
    MM_FIELD field;
    MM_SYMMETRY symmetry;

    //Read start tag
    ss >> tmp;
    toLower(tmp);
    if(tmp != "%%matrixmarket") return is;

    //Read object tag
    ss >> tmp;
    toLower(tmp);
    if(tmp != "matrix") return is;

    //Read data format tag
    ss >> s_format;
    toLower(s_format);
    if(s_format != "coordinate" && s_format != "array")
      return is;
    

    //Read number field tag
    ss >> s_field;
    toLower(s_field);
    if(s_field == "real")
      field = mm_real;
    else if(s_field == "integer")
      field = mm_integer;
    else if(s_field == "complex")
      field = mm_complex;
    else if(s_field == "pattern") 
      field = mm_pattern;
    else if(s_field == "interval") 
      field = mm_interval;
    else if(s_field == "cinterval") 
      field = mm_cinterval;
    else 
      return is;

    //Read symmetry tag
    ss >> s_symmetry;
    toLower(s_symmetry);
    if(s_symmetry == "general")
      symmetry = mm_general;
    else if(s_symmetry == "symmetric")
      symmetry = mm_symmetric;
    else if(s_symmetry == "skew-symmetric")
      symmetry = mm_skew_symmetric;
    else if(s_symmetry == "hermitian")
      symmetry = mm_hermitian;
    else
      return is;

    //Read past commentary
    std::getline(is,tmp);
    while(tmp[0] == '%')
      std::getline(is,tmp);

    //Read past empty lines
    while(tmp.size() == 0)
      std::getline(is,tmp);

  
    //Read size information
    ss.clear();
    ss.str(tmp);
    ss >> A.m >> A.n;
    int nnz; 
    ss >> nnz;

    A.lb1=1; A.ub1=A.m;
    A.lb2=1; A.ub2=A.n;

    //Clear matrix
    A.p = std::vector<int>((A.n>0) ? A.n+1 : 1, 0);
    A.ind.clear();
    A.x.clear();

    //Read past empty lines
    while(tmp.size() == 0)
      std::getline(is,tmp);

    //Read Data
    std::vector<triplet_store<TType> > work;
    work.reserve(nnz);

    if(s_format == "coordinate") {
      //Sparse data
      while(is.good()) {
	std::getline(is,tmp);
	if(tmp.size() != 0) {
	  mm_add_element(work,tmp,field,symmetry);
	}
      }
    } else {
      //Dense data
      int i=1, j=1;
      while(is.good()) {
	std::getline(is,tmp);
	if(tmp.size() != 0) {
          std::stringstream ss;
          ss << i << j << tmp;	
          tmp = ss.str();	
	  mm_add_element(work,tmp,field,symmetry);
          i++;
          if(i == A.m) {
            i = 1;
            j++;
          }
	}
      }
    }

    sort(work.begin(), work.end());

    int i=0;

    for(int j=0 ; j<A.n ; j++) {        

      while((unsigned int)i < work.size() && work[i].col == j ) {
        A.ind.push_back(work[i].row);
        A.x.push_back(work[i].val);
        i++;
      }

      A.p[j+1] = i;
    }

  }

  return is;
}

//--------------------------------------------------------------------------

template<class TA, class TB, class TElement>
inline TA& slsp_mm_assign(TA& A, const TB& C) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(C) || ColLen(A)!=ColLen(C)) cxscthrow(OP_WITH_WRONG_DIM("operator=(" + nameof(A) + " &, const " + nameof(C) + " &)"));
#endif
    int lb1=A.A.lb1, ub1=A.A.ub1, lb2=A.A.lb2, ub2=A.A.ub2;
    A.A = C;
    A.A.lb1=lb1; A.A.ub1=ub1; A.A.lb2=lb2; A.A.ub2=ub2;

    int start1 = A.A.lb1 - A.M->lb1;
    int end1 = A.A.ub1 - A.M->lb1;
    int start2 = A.A.lb2 - A.M->lb2;
    int end2 = A.A.ub2 - A.M->lb2;
    int tmp1=0, tmp2=0;

    for(int j=start2 ; j<=end2 ; j++) {
      //1. alte loeschen zwischen start1 und end1 
      std::vector<int>::iterator ind_it = A.M->ind.begin()+A.M->p[j];
      TElement x_it  = A.M->x.begin()+A.M->p[j];

      int size = A.M->p[j+1] - tmp1 + tmp2 - A.M->p[j];

      for(int k=0 ; k<size ; k++) {
        if(*ind_it>=start1 && *ind_it<=end1) {
          ind_it = A.M->ind.erase(ind_it);
          x_it = A.M->x.erase(x_it);
          tmp1++;
        } else {
          x_it++;
          ind_it++;
        }
      }

      A.M->p[j+1] -= tmp1;

      ind_it = A.M->ind.begin()+A.M->p[j];
      x_it  = A.M->x.begin()+A.M->p[j];

      //2. neue aus C kopieren zwischen start1 und end1
      for(int k=A.A.p[j-start2+1]-1 ; k>=A.A.p[j-start2] ; k--) {
        if(A.A.ind[k]>=0 && A.A.ind[k]<=end1-start1) {
	  //Sortierung pro Spalte beachten
	  while(ind_it < A.M->ind.begin()+A.M->p[j+1]+tmp2   &&  A.A.ind[k]+start1 > *ind_it) {
	    ind_it++;
	    x_it++;
	  }
          ind_it = A.M->ind.insert(ind_it, A.A.ind[k]+start1);
          x_it = A.M->x.insert(x_it, A.A.x[k]);
          tmp2++;
        }
      }

      A.M->p[j+1] += tmp2;
    }

    for(unsigned int i=end2+2 ; i<A.M->p.size() ; i++) {
      A.M->p[i] += tmp2 - tmp1;
    }
    

    return A;
}

template<class TA, class TB, class TElement, class TType>
inline TA& slf_mm_assign(TA& A, const TB& C) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(C) || ColLen(A)!=ColLen(C)) cxscthrow(OP_WITH_WRONG_DIM("operator=(" + nameof(A) + " &, const " + nameof(C) + " &)"));
#endif
      int lb1=A.A.lb1, ub1=A.A.ub1, lb2=A.A.lb2, ub2=A.A.ub2;
      A.A = C;
      A.A.lb1=lb1; A.A.ub1=ub1; A.A.lb2=lb2; A.A.ub2=ub2;

      int start1 = A.A.lb1 - A.M->lb1;
      int end1 = A.A.ub1 - A.M->lb1;
      int start2 = A.A.lb2 - A.M->lb2;
      int end2 = A.A.ub2 - A.M->lb2;

      for(int j=start2 ; j<=end2 ; j++) {
        //1. alte loeschen zwischen start1 und end1 
        std::vector<int>::iterator ind_it = A.M->ind.begin()+A.M->p[j];
        TElement x_it  = A.M->x.begin()+A.M->p[j];

        int size = A.M->p[j+1]-A.M->p[j];

        for(int k=0 ; k<size ; k++) {
          if(*ind_it>=start1 && *ind_it<=end1) {
            ind_it = A.M->ind.erase(ind_it);
            x_it = A.M->x.erase(x_it);
            for(unsigned int i=j+1 ; i<A.M->p.size() ; i++) {
              A.M->p[i]--;
            }
          } else {
            x_it++;
            ind_it++;
          }
        }

        //2. neue aus C kopieren zwischen start1 und end1
        for(int k=end1 ; k>=start1 ; k--) {
          if(C[k-start1+Lb(C,1)][j-start2+Lb(C,2)] != 0.0) {
	    ind_it = A.M->ind.insert(ind_it, k);
            x_it = A.M->x.insert(x_it, TType(C[k-start1+Lb(C,1)][j-start2+Lb(C,2)]));
            for(unsigned int i=j+1 ; i<A.M->p.size() ; i++) {
              A.M->p[i]++;
            }
          }
        }
      }

      return A;
}

template<class TA, class TB, class TType>
inline TA& spf_mm_assign(TA& A, const TB& B) 
{
      A.m = ColLen(B);
      A.n = RowLen(B);
      A.lb1 = Lb(B,1);
      A.ub1 = Ub(B,1);
      A.lb2 = Lb(B,2);
      A.ub2 = Ub(B,2);
      int nnz = 0;
      A.p = std::vector<int>(A.n+1,0);
      for(unsigned int i=0 ; i<A.p.size() ; i++)
        A.p[i] = 0;

      A.ind.clear();
      A.x.clear();

      for(int j=0 ; j<A.n ; j++) {
        for(int i=0 ; i<A.m ; i++) {
          if(B[i+Lb(B,1)][j+Lb(B,2)] != 0.0) {
             A.ind.push_back(i);
             A.x.push_back(TType(B[i+Lb(B,1)][j+Lb(B,2)]));
             nnz++;
          }
        }
          
        A.p[j+1] = nnz;
      }

      return A;
}

template<class TA, class Ts, class TType>
inline TA& sp_ms_assign(TA& A, const Ts& s) {
      int nnz = 0;      
      for(unsigned int i=0 ; i<A.p.size() ; i++)
        A.p[i] = 0;

      A.ind.clear();
      A.x.clear();

      if(s != 0.0) {
        for(int j=0 ; j<A.n ; j++) {
          for(int i=0 ; i<A.m ; i++) {
             A.ind.push_back(i);
             A.x.push_back(TType(s));
             nnz++;
          }
          A.p[j+1] = nnz; 
        }
          
      }

      return A;
}

template<class TA, class Ts, class TElement, class TType>
inline TA& sl_ms_assign(TA& A, const Ts& s) {
      A.A = s;

      int start1 = A.A.lb1 - A.M->lb1;
      int end1 = A.A.ub1 - A.M->lb1;
      int start2 = A.A.lb2 - A.M->lb2;
      int end2 = A.A.ub2 - A.M->lb2;

      for(int j=start2 ; j<=end2 ; j++) {
        //1. alte loeschen zwischen start1 und end1 
        std::vector<int>::iterator ind_it = A.M->ind.begin()+A.M->p[j];
        TElement x_it  = A.M->x.begin()+A.M->p[j];

        int size = A.M->p[j+1]-A.M->p[j];

        for(int k=0 ; k<size ; k++) {
          if(*ind_it>=start1 && *ind_it<=end1) {
            ind_it = A.M->ind.erase(ind_it);
            x_it = A.M->x.erase(x_it);
            for(unsigned int i=j+1 ; i<A.M->p.size() ; i++) {
              A.M->p[i]--;
            }
          } else {
            x_it++;
            ind_it++;
          }
        }

        //2. neue Werte einfuegen
        if(s != 0.0) {
          for(int k=end1 ; k>=start1 ; k--) {
            ind_it = A.M->ind.insert(ind_it, k);
            x_it = A.M->x.insert(x_it, TType(s));
            for(unsigned int i=j+1 ; i<A.M->p.size() ; i++) {
              A.M->p[i]++;
            }
          }
        }
      }

      return A;
}


//---------------------------------------------------------------------

template<class TA, class Ts>
inline TA& sp_ms_divassign(TA& A, const Ts& r) {
  if(r == 0.0) {
    A = TA();
  } else {
    for(unsigned int i=0 ; i<A.x.size() ; i++) 
      A.x[i] /= r;
  }
  return A;
}

template<class TA, class Ts>
inline TA& sp_ms_multassign(TA& A, const Ts& r) {
  if(r == 0.0) {
    A = TA();
  } else {
    for(unsigned int i=0 ; i<A.x.size() ; i++) 
      A.x[i] *= r;
  }
  return A;
}

template<class TA, class TB, class TDot, class TElement>
inline TA& spsp_mm_multassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=ColLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator*=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int m = ColLen(A);
  //int n = RowLen(A);
  int o = RowLen(B);

  TA C(m,o,A.get_nnz()+B.get_nnz());

  if(opdotprec == 1 && isReal<TElement>()) {
    std::vector<TElement> dot(m,TElement(0.0));
    std::vector<int> w(m,-1);
    int nnz = 0;

    for(int i=0 ; i<o ; i++) {

      for(int k=B.p[i] ; k<B.p[i+1] ; k++) {
        for(int l=A.p[B.ind[k]] ; l<A.p[B.ind[k]+1] ; l++) {
          if(w[A.ind[l]] < i) {
            w[A.ind[l]] = i;
            C.ind.push_back(A.ind[l]);
            dot[A.ind[l]] = A.x[l] * B.x[k];
            nnz++;
          } else {
            dot[A.ind[l]] += A.x[l] * B.x[k];
          }
        }
      }

      sort(C.ind.begin()+C.p[i], C.ind.begin()+nnz);

      for(int j=C.p[i] ; j<nnz ; j++) {
         C.x.push_back(dot[C.ind[j]]);
      }

      C.p[i+1] = nnz;
    }

  } else if(opdotprec == 1 && isComplex<TElement>()) {

    std::vector<TElement> dot(m,TElement(0.0));
    std::vector<int> w(m,-1);
    int nnz = 0;

    for(int i=0 ; i<o ; i++) {

      for(int k=B.p[i] ; k<B.p[i+1] ; k++) {
        for(int l=A.p[B.ind[k]] ; l<A.p[B.ind[k]+1] ; l++) {
          if(w[A.ind[l]] < i) {
            w[A.ind[l]] = i;
            C.ind.push_back(A.ind[l]);
            //dot[A.ind[l]] = A.x[l] * B.x[k];
            fp_mult(A.x[l],B.x[k],dot[A.ind[l]]);
            nnz++;
          } else {
            //dot[A.ind[l]] += A.x[l] * B.x[k];
            fp_multadd(A.x[l],B.x[k],dot[A.ind[l]]);
          }
        }
      }

      sort(C.ind.begin()+C.p[i], C.ind.begin()+nnz);

      for(int j=C.p[i] ; j<nnz ; j++) {
         C.x.push_back(dot[C.ind[j]]);
      }

      C.p[i+1] = nnz;
    }

  } else if(opdotprec != 0 || !SPMULT_SAVE_MEMORY) {

    std::vector<TDot> dot(m,TDot(opdotprec));
    std::vector<int> w(m,-1);
    int nnz = 0;

    for(int i=0 ; i<o ; i++) {

      for(int k=B.p[i] ; k<B.p[i+1] ; k++) {
        for(int l=A.p[B.ind[k]] ; l<A.p[B.ind[k]+1] ; l++) {
          if(w[A.ind[l]] < i) {
            w[A.ind[l]] = i;
            C.ind.push_back(A.ind[l]);
            dot[A.ind[l]].reset();
            dot[A.ind[l]].add_dot(A.x[l],B.x[k]);
            nnz++;
          } else {
            dot[A.ind[l]].add_dot(A.x[l],B.x[k]);
          }
        }
      }

      sort(C.ind.begin()+C.p[i], C.ind.begin()+nnz);

      for(int j=C.p[i] ; j<nnz ; j++) {
         C.x.push_back(dot[C.ind[j]].result());
      }

      C.p[i+1] = nnz;
    }


  } else {
  
    TA At(A.n,A.m);
    At = transp(A);

    TDot dot(opdotprec);

    int nnz = 0;

    for(int i=0 ; i<o ; i++) {
      for(int j=0 ; j<m ; j++) {
        dot.reset();

        for(int k=At.p[j] ; k<At.p[j+1] ; k++) {
          for(int l=B.p[i] ; l<B.p[i+1] && At.ind[k] >= B.ind[l] ; l++) {
            if(At.ind[k] == B.ind[l]) {
              dot.add_dot(At.x[k], B.x[l]);
            }
          }
        }

        TElement result = dot.result();
        if(result != 0.0) {
          C.x.push_back(result);
          C.ind.push_back(j);
          nnz++;
        }
      }

      C.p[i+1] = nnz;
    }

  }

  A = C;

  return A;
}

template<class TA, class TB, class TDot, class TFull>
inline TA& spf_mm_multassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=ColLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator*=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int lb1 = Lb(B,1);
  int lb2 = Lb(B,2);

  int m = ColLen(A);
  //int n = RowLen(A);
  int o = RowLen(B);

  TFull C(m,o);

  C = 0.0;

  TA At(A.n,A.m);
  At = transp(A);

  if(opdotprec == 1 && isReal<TDot>()) {

    for(int i=0 ; i<o ; i++) {
      for(int j=0 ; j<m ; j++) {
        for(int k=At.p[j] ; k<At.p[j+1] ; k++) {
           C[j+1][i+1] += At.x[k] * B[At.ind[k]+lb1][i+lb2];
        }
      }
    }

  } else {
    TDot dot(opdotprec);

    for(int i=0 ; i<o ; i++) {
      for(int j=0 ; j<m ; j++) {
        dot.reset();
        for(int k=At.p[j] ; k<At.p[j+1] ; k++) {
           dot.add_dot(At.x[k], B[At.ind[k]+lb1][i+lb2]);
        }
        C[j+1][i+1] = dot.result();
      }
    }
  }

  A = TA(C);

  return A;
}

template<class TA, class TB, class TDot, class TFull>
inline TA& fsp_mm_multassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=ColLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator*=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int lb1 = Lb(A,1);
  int lb2 = Lb(A,2);

  int m = ColLen(A);
  //int n = RowLen(A);
  int o = RowLen(B);

  TFull C(m,o);

  C = 0.0;

  if(opdotprec == 1 && isReal<TDot>()) {

    for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
        for(int k=B.p[j] ; k<B.p[j+1] ; k++) {
           C[i+1][j+1] = A[i+lb1][B.ind[k]+lb2] * B.x[k];
        }
      }
    }

  } else {

    TDot dot(opdotprec);

    for(int i=0 ; i<m ; i++) {
      for(int j=0 ; j<o ; j++) {
        dot.reset();
        for(int k=B.p[j] ; k<B.p[j+1] ; k++) {
           dot.add_dot(A[i+lb1][B.ind[k]+lb2], B.x[k]);
        }
        C[i+1][j+1] = dot.result();
      }
    }
  }

  A = C;

  return A;
}


//------------------------------------------------------------

template<class TA, class TB, class TFull>
inline TA& spf_mm_addassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator+=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  //int n = RowLen(A);
  
  TFull DA;
  A.full(DA);

  DA += B;
  
  A = TA(DA);

  return A;
}

template<class TA, class TB, class TElement>
inline TA& spsp_mm_addassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator+=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int m = ColLen(A);
  int n = RowLen(A);


  TA C(m, n, A.get_nnz()+B.get_nnz());

  int nnz = 0;


  for(int j=0 ; j<n ; j++) {

    int k = A.p[j];
    int l = B.p[j]; 

    while(k<A.p[j+1] && l<B.p[j+1]) {
      if(A.ind[k] == B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(A.x[k] + B.x[l]);
        k++; l++;
      } else if(A.ind[k] < B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(TElement(A.x[k]));
        k++;
      } else {
        C.ind.push_back(B.ind[l]);
        C.x.push_back(TElement(B.x[l]));
        l++;
      }
      nnz++;
    }

    for( ; k<A.p[j+1] ; k++) {
      C.ind.push_back(A.ind[k]);
      C.x.push_back(TElement(A.x[k]));
      nnz++;
    }

    for( ; l<B.p[j+1] ; l++) {
      C.ind.push_back(B.ind[l]);
      C.x.push_back(TElement(B.x[l]));
      nnz++;
    }

    C.p[j+1] = nnz;

  }

  A = C;

  return A;
}

template<class TA, class TB>
inline TA& fsp_mm_addassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator+=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  int n = RowLen(A);
  int lb1 = Lb(A,1), lb2 = Lb(A,2);
  for(int j=0 ; j<n ; j++) {
    for(int k=B.p[j] ; k<B.p[j+1] ; k++) {
      A[B.ind[k]+lb1][j+lb2] += B.x[k];
    }
  }

  return A;
}


template<class TA, class TB, class TFull>
inline TA& spf_mm_subassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator-=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  //int n = RowLen(A);
  
  TFull DA;
  A.full(DA);

  DA -= B;
  
  A = TA(DA);

  return A;
}

template<class TA, class TB, class TElement>
inline TA& spsp_mm_subassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator-=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int m = ColLen(A);
  int n = RowLen(A);

  TA C(m, n, A.get_nnz()+B.get_nnz());

  int nnz = 0;

  for(int j=0 ; j<n ; j++) {

    int k = A.p[j];
    int l = B.p[j]; 

    while(k<A.p[j+1] && l<B.p[j+1]) {
      if(A.ind[k] == B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(A.x[k] - B.x[l]);
        k++; l++;
      } else if(A.ind[k] < B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(TElement(A.x[k]));
        k++;
      } else {
        C.ind.push_back(B.ind[l]);
        C.x.push_back(TElement(-B.x[l]));
        l++;
      }
      nnz++;
    }

    for( ; k<A.p[j+1] ; k++) {
      C.ind.push_back(A.ind[k]);
      C.x.push_back(TElement(A.x[k]));
      nnz++;
    }

    for( ; l<B.p[j+1] ; l++) {
      C.ind.push_back(B.ind[l]);
      C.x.push_back(TElement(-B.x[l]));
      nnz++;
    }

    C.p[j+1] = nnz;

  }

  A = C;

  return A;
}

template<class TA, class TB>
inline TA& fsp_mm_subassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator-=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  int n = RowLen(A);
  int lb1=Lb(A,1), lb2=Lb(A,2);
  
  for(int j=0 ; j<n ; j++) {
    for(int k=B.p[j] ; k<B.p[j+1] ; k++) {
      A[B.ind[k]+lb1][j+lb2] -= B.x[k];
    }
  }

  return A;
}

//------------------------------------------------------------

template<class TA, class TB, class TFull>
inline TA& spf_mm_hullassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator|=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  //int n = RowLen(A);
  
  TFull DA;
  A.full(DA);

  DA |= B;
  
  A = TA(DA);

  return A;
}

template<class TA, class TB, class TElement>
inline TA& spsp_mm_hullassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator|=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int m = ColLen(A);
  int n = RowLen(A);

  TA C(m, n, A.get_nnz()+B.get_nnz());

  int nnz = 0;

  for(int j=0 ; j<n ; j++) {

    int k = A.p[j];
    int l = B.p[j]; 

    while(k<A.p[j+1] && l<B.p[j+1]) {
      if(A.ind[k] == B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(A.x[k] | B.x[l]);
        k++; l++;
      } else if(A.ind[k] < B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(TElement(A.x[k]));
        k++;
      } else {
        C.ind.push_back(B.ind[l]);
        C.x.push_back(TElement(B.x[l]));
        l++;
      }
      nnz++;
    }

    for( ; k<A.p[j+1] ; k++) {
      C.ind.push_back(A.ind[k]);
      C.x.push_back(TElement(A.x[k]));
      nnz++;
    }

    for( ; l<B.p[j+1] ; l++) {
      C.ind.push_back(B.ind[l]);
      C.x.push_back(TElement(B.x[l]));
      nnz++;
    }

    C.p[j+1] = nnz;

  }

  A = C;

  return A;
}

template<class TA, class TB>
inline TA& fsp_mm_hullassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator|=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  int n = RowLen(A);
  int lb1=Lb(A,1), lb2=Lb(A,2);
  
  for(int j=0 ; j<n ; j++) {
    for(int k=B.p[j] ; k<B.p[j+1] ; k++) {
      A[B.ind[k]+lb1][j+lb2] |= B.x[k];
    }
  }

  return A;
}

//------------------------------------------------------------

template<class TA, class TB, class TFull>
inline TA& spf_mm_intersectassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator&=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  //int n = RowLen(A);
  
  TFull DA;
  A.full(DA);

  DA &= B;
  
  A = TA(DA);

  return A;
}

template<class TA, class TB, class TElement>
inline TA& spsp_mm_intersectassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator&=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  int m = ColLen(A);
  int n = RowLen(A);

  TA C(m, n, A.get_nnz()+B.get_nnz());
  std::vector<TElement> w(m,TElement(0.0));

  int nnz = 0;

  for(int j=0 ; j<n ; j++) {

    int k = A.p[j];
    int l = B.p[j]; 

    while(k<A.p[j+1] && l<B.p[j+1]) {
      if(A.ind[k] == B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(A.x[k] & B.x[l]);
        k++; l++;
      } else if(A.ind[k] < B.ind[l]) {
        C.ind.push_back(A.ind[k]);
        C.x.push_back(TElement(A.x[k]));
        k++;
      } else {
        C.ind.push_back(B.ind[l]);
        C.x.push_back(TElement(B.x[l]));
        l++;
      }
      nnz++;
    }

    for( ; k<A.p[j+1] ; k++) {
      C.ind.push_back(A.ind[k]);
      C.x.push_back(TElement(A.x[k]));
      nnz++;
    }

    for( ; l<B.p[j+1] ; l++) {
      C.ind.push_back(B.ind[l]);
      C.x.push_back(TElement(B.x[l]));
      nnz++;
    }

    C.p[j+1] = nnz;

  }

  A = C;

  return A;
}

template<class TA, class TB>
inline TA& fsp_mm_intersectassign(TA& A, const TB& B) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(RowLen(A)!=RowLen(B) || ColLen(A)!=ColLen(B)) cxscthrow(OP_WITH_WRONG_DIM("operator&=(" + nameof(A) + " &, const " + nameof(B) + " &)"));
#endif
  //int m = ColLen(A);
  int n = RowLen(A);
  int lb1=Lb(A,1), lb2=Lb(A,2);
  
  for(int j=0 ; j<n ; j++) {
    for(int k=B.p[j] ; k<B.p[j+1] ; k++) {
      A[B.ind[k]+lb1][j+lb2] &= B.x[k];
    }
  }

  return A;
}

//------------------------------------------------------------------------

template<class Tx, class Ty>
inline Tx& svsp_vv_assign(Tx& x, const Ty& y) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(x)!=VecLen(y)) cxscthrow(OP_WITH_WRONG_DIM("operator=(" + nameof(x) + " &, const " + nameof(y) + " &)"));
#endif
  if(x.row) {
    for(unsigned int i=0 ; i<y.p.size() ; i++)
       x.dat[x.index][y.p[i]+x.dat.A.lb2] = y.x[i];
  } else {
    for(unsigned int i=0 ; i<y.p.size() ; i++)
       x.dat[y.p[i]+x.dat.A.lb1][x.index] = y.x[i];         
  }
  return x;
}

template<class Tx, class Ty>
inline Tx& svsl_vv_assign(Tx& x, const Ty& y) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(x)!=VecLen(y)) cxscthrow(OP_WITH_WRONG_DIM("operator=(" + nameof(y) + " &, const " + nameof(y) + " &)"));
#endif
  if(x.row) {
    for(unsigned int i=0 ; i<y.p.size() ; i++)
       x.dat[x.index][y.p[i]+x.dat.A.lb2] = y.x[i];
  } else {
    for(unsigned int i=0 ; i<y.p.size() ; i++)
       x.dat[y.p[i]+x.dat.A.lb1][x.index] = y.x[i];         
  }
  return x;
}

template<class Tx, class Ty>
inline Tx& svf_vv_assign(Tx& x, const Ty& y) 
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
  if(VecLen(x)!=VecLen(y)) cxscthrow(OP_WITH_WRONG_DIM("operator=(" + nameof(x) + " &, const " + nameof(y) + " &)"));
#endif
  if(x.row) {
    for(int i=0 ; i<x.dat.A.n ; i++)
       x.dat[x.index][i+x.dat.A.lb2] = y[i+Lb(y)];
  } else {
    for(int i=0 ; i<x.dat.A.m ; i++)
       x.dat[i+x.dat.A.lb1][x.index] = y[i+Lb(y)];         
  }
  return x;
}

template<class Tx, class Ts>
inline Tx& sv_vs_assign(Tx& x, const Ts& s) {
  if(x.row) {
    for(int i=0 ; i<x.dat.A.n ; i++)
       x.dat[x.index][i+x.dat.A.lb2] = s;
  } else {
    for(int i=0 ; i<x.dat.A.m ; i++)
       x.dat[i+x.dat.A.lb1][x.index] = s;         
  }
  return x;
}

//-----------------------------------------------------------------

template<class TA>
inline bool sp_m_not(const TA& M) {
  bool ret = true;
  for(int i=0 ; i<M.get_nnz() ; i++)
    ret = ret && (!M.x[i]);
  return ret;
}

template<class Tx>
inline bool sv_v_not(const Tx& x) {
  bool ret = true;
  if(x.row) {
    for(int i=0 ; i<x.dat.A.n ; i++)
       ret = ret && (!x.dat(x.index,i+x.dat.A.lb2));
  } else {
    for(int i=0 ; i<x.dat.A.m ; i++)
       ret = ret && (!x.dat(i+x.dat.A.lb1,x.index));         
  }
  return ret;
}

//-------------------------------------------------------------------

template <class TA>
inline void sp_m_resize(TA& A) throw() {
    A.p.clear();
    A.ind.clear();
    A.x.clear();
    A.m = A.n = 0;
    A.lb1 = A.lb2 = 1;
    A.ub1 = A.ub2 = 0;
}

template <class TA>
inline void sp_m_resize(TA &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
    throw(WRONG_BOUNDARIES)
#else
    throw()
#endif
{
#if(CXSC_INDEX_CHECK)
    if((m<0)||(n<0)) cxscthrow(WRONG_BOUNDARIES("void Resize("+nameof(A)+" &,const int &, const int &)"));
#endif
    TA tmp(m,n);

    if(m <= A.m && n <= A.n) 
      tmp = A(Lb(A,1), Lb(A,1)+m-1, Lb(A,2), Lb(A,2)+n-1);
    else if(m <= A.m && n >= A.n)
      tmp(1,m,1,A.n) = A(Lb(A,1), Lb(A,1)+m-1, Lb(A,2), Ub(A,2));
    else if(m >= A.m && n <= A.n)
      tmp(1,A.m,1,n) = A(Lb(A,1), Ub(A,1), Lb(A,2), Lb(A,2)+n-1);
    else if(m >= A.m && n >= A.n)    
      tmp(1,A.m,1,A.n) = A;

    A.m = m;
    A.n = n;
    A.lb1 = A.lb2 = 1;
    A.ub1 = m;
    A.ub2 = n;

    A = tmp;
}

template<class TA>
inline void sp_m_resize(TA &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
    throw(WRONG_BOUNDARIES)
#else
    throw()
#endif
{
#if(CXSC_INDEX_CHECK)
    if((m2<m1)||(n2<n1)) cxscthrow(WRONG_BOUNDARIES("void Resize("+nameof(A)+" &,const int &, const int &,const int &,const int &)"));
#endif
    sp_m_resize(A, m2-m1+1, n2-n1+1);
    SetLb(A,ROW,m1); SetUb(A,ROW,m2);
    SetLb(A,COL,n1); SetUb(A,COL,n2);
}

} //namespace cxsc 

#endif
