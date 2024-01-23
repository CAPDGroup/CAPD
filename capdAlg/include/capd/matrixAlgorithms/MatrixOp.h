/////////////////////////////////////////////////////////////////////////////
/// @file MatrixOp.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-19
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details.

#ifndef CAPD_FILE_MATRIXOP_H
#define CAPD_FILE_MATRIXOP_H

#include <capd/vectalg/Matrix.h>

namespace capd
{

namespace matrixAlgorithms
{
/* ------------------------  ------------------------ */
template<typename intType>
inline bool isDivisible(intType a, intType b){
return a % b == intType(0);
}

inline bool isDivisible(double /*a*/, double /*b*/){
return true;
}

/* ------------------------  ------------------------ */
template<typename intType>
inline bool isInvertible(intType a){
return (a == intType(1)) || (a == -intType(1));
}

/* ------------------------  ------------------------ */
template<typename intType>
inline intType inverse(intType a){
return intType(1)/a;
}


template<typename matrix>
capd::vectalg::Matrix<typename matrix::ScalarType, 0, 0> emptyMatrix()
{
return capd::vectalg::Matrix<typename matrix::ScalarType, 0, 0>();
}
/* Elementary row and column operations */

/* ------------------------  ------------------------ */
template<class matrix>
void rowExchange(matrix& A,int i,int j){
if(i==j || A.empty()) return;
MatrixIterator<matrix> it(A),it2(A);
for(
it=A.beginOfRow(i),it2=A.beginOfRow(j);it<A.endOfRow(i);
it.moveToNextColumn(),it2.moveToNextColumn()
   ){
typename matrix::ScalarType s=*it;
*it=*it2;
*it2=s;
}
}

/* ------------------------  ------------------------ */
template<class matrix>
void rowMultiply(matrix& A,int i,typename matrix::ScalarType s){
if (A.empty()) return;

MatrixIterator<matrix> it(A);
for(
it=A.beginOfRow(i);
it<A.endOfRow(i);
it.moveToNextColumn()
){
*it*=s;
}
}

/* ------------------------  ------------------------ */
template<class matrix>
void rowAdd(matrix& A,int i,int j,typename matrix::ScalarType s){
if (A.empty()) return;

MatrixIterator<matrix> it(A),it2(A);
for(
it=A.beginOfRow(i),it2=A.beginOfRow(j);
it<A.endOfRow(i);
it.moveToNextColumn(),it2.moveToNextColumn()
   ){
*it+=s**it2;
}
}

/* ------------------------  ------------------------ */
template<class matrix>
void columnExchange(matrix& A,int i,int j){
if (A.empty()) return;

if(i==j) return;
MatrixIterator<matrix> it(A),it2(A);
for(
it=A.beginOfColumn(i),it2=A.beginOfColumn(j);
it<A.endOfColumn(i);
it.moveToNextRow(),it2.moveToNextRow()
   ){
typename matrix::ScalarType s=*it;
*it=*it2;
*it2=s;
}
}

/* ------------------------  ------------------------ */
template<class matrix>
void columnMultiply(matrix& A,int j,typename matrix::ScalarType s){
if (A.empty()) return;

MatrixIterator<matrix> it2(A);
for(
it2=A.beginOfColumn(j);
it2<A.endOfColumn(j);
it2.moveToNextRow()
){
*it2*=s;
}
}

/* ------------------------  ------------------------ */
template<class matrix>
void columnAdd(matrix& A,int i,int j,typename matrix::ScalarType s){
if (A.empty()) return;

MatrixIterator<matrix> it(A),it2(A);
for(
it=A.beginOfColumn(i),it2=A.beginOfColumn(j);
it<A.endOfColumn(i);
it.moveToNextRow(),it2.moveToNextRow()
   ){
*it2+=s**it;
}
}

/* Elementary row and column operations on matrix and matrices of bases */
/* ------------------------ ------------------------ */
template<class matrix, class sqMatrix>
void rowAdd(matrix& B,sqMatrix& Q,int i,int j,typename matrix::ScalarType q){
rowAdd(B,i,j,q);
columnAdd(Q,i,j,-q);
}

/* ------------------------  ------------------------ */
template<class matrix, class sqMatrix>
void columnExchange(matrix& B,sqMatrix& R,int i,int j){
columnExchange(B,i,j);
columnExchange(R,i,j);
}

/* ------------------------  ------------------------ */
template<class matrix, class sqMatrix>
void columnMultiply(matrix& B,sqMatrix& R,int i,typename matrix::ScalarType q){
columnMultiply(B,i,q);
columnMultiply(R,i,q);
}

/* ------------------------  ------------------------ */
template<class matrix, class sqMatrix>
void columnAdd(matrix& B,sqMatrix& R,int i,int j,typename matrix::ScalarType q){
columnAdd(B,i,j,q);
columnAdd(R,i,j,q);
}

/* ------------------------  ------------------------ */

template<class matrix, class sqMatrix>
void rowExchange(matrix& B,sqMatrix& Q,sqMatrix& Qinv,int i,int j){
rowExchange(B,i,j);
rowExchange(Qinv,i,j);
columnExchange(Q,i,j);
}

/* ------------------------  ------------------------ */
template<class matrix, class sqMatrix>
void rowMultiply(matrix& B,sqMatrix& Q,sqMatrix& Qinv,int i,typename matrix::ScalarType q){
rowMultiply(B,i,q);
rowMultiply(Qinv,i,q);
columnMultiply(Q,i,q);
}

/* ------------------------  ------------------------ */
template<class matrix, class sqMatrix>
void rowAdd(matrix& B,sqMatrix& Q,sqMatrix& Qinv,int i,int j,typename matrix::ScalarType q){
rowAdd(B,i,j,q);
rowAdd(Qinv,i,j,q);
 columnAdd(Q,i,j,-q);
}

  /* ------------------------  ------------------------ */
  template<class matrix, class sqMatrix>
  void columnExchange(matrix& B,sqMatrix& R,sqMatrix& Rinv,int i,int j){
    columnExchange(B,i,j);
    columnExchange(R,i,j);
    rowExchange(Rinv,i,j);
  }

  /* ------------------------  ------------------------ */
  template<class matrix, class sqMatrix>
  void columnMultiply(matrix& B,sqMatrix& R,sqMatrix& Rinv,int i,typename matrix::ScalarType q){
    columnMultiply(B,i,q);
    columnMultiply(R,i,q);
    rowMultiply(Rinv,i,q);
  }

  /* ------------------------  ------------------------ */
  template<class matrix, class sqMatrix>
  void columnAdd(matrix& B,sqMatrix& R,sqMatrix& Rinv,int i,int j,typename matrix::ScalarType q){
    columnAdd(B,i,j,q);
    columnAdd(R,i,j,q);
    rowAdd(Rinv,i,j,-q);
  }

  /* ------------------------  ------------------------ */

  template<class matrix, class sqMatrix>
  void partRowReduce(matrix& B,sqMatrix& Q,sqMatrix& Qinv,int k,int l){
    int m=B.numberOfRows();
    for(int i=k+1;i<=m;i++){
      typename matrix::ScalarType q=B(i,l)/B(k,l);
      rowAdd(B,Q,Qinv,i,k,-q);
    }
  }

  /* ------------------------  ------------------------ */
  template<class matrix, class sqMatrix>
  void partColumnReduce(matrix& B,sqMatrix& R,sqMatrix& Rinv,int k,int l){
    int n=B.numberOfColumns();
    for(int j=l+1;j<=n;j++){
      typename matrix::ScalarType q=B(k,j)/B(k,l);
      columnAdd(B,R,Rinv,l,j,-q);
    }
  }

  // *** Test for nonzero matrices *** /

  /* ------------------------  ------------------------ */
  template<class matrix>
  void smallestNonZero(const matrix& A,typename matrix::ScalarType& s, int& iOpt, int& jOpt){
    typedef typename matrix::ScalarType ScalarType;
    s=ScalarType(0);
    const_MatrixIterator<matrix> it(A),itOpt(A);
    for(int i=1;i<=(int)A.numberOfRows();++i){
      it=A.beginOfRow(i);
      for(int j=1;j<=(int)A.numberOfColumns();++j){
	ScalarType t=*it;
	if(t<ScalarType(0)) t=-t;
	if( s==ScalarType(0) || (s>t && t>ScalarType(0)) ){
	  s=t;itOpt=it;
	}
	it.moveToNextColumn();
      }
    }
    std::pair<int,int> p=itOpt.rowAndColumn();
    iOpt=p.first;
    jOpt=p.second;
  }

  /* ------------------------  ------------------------ */
  template<class matrix>
  bool nonZero(const matrix& A){
    typedef typename matrix::ScalarType ScalarType;
    const_MatrixIterator<matrix> it(A);
    for(int i=1;i<=(int)A.numberOfRows();++i){
      it=A.beginOfRow(i);
      for(int j=1;j<=(int)A.numberOfColumns();++j){
	if(*it!=ScalarType(0)) return true;
	it.moveToNextColumn();
      }
    }
    return false;
  }

  // *** row echelon form *** //

  /* ------------------------  ------------------------ */
  template<class matrix, class sqMatrix>
  void rowPrepare(matrix& B,sqMatrix& Q,sqMatrix& Qinv,int k,int l){
    typedef typename matrix::ScalarType ScalarType;
    int m=B.numberOfRows();
    ScalarType s;
    int i,j;
    ;
    smallestNonZero(MatrixSlice<matrix>(B,k,m,l,l),s,i,j);
    i+=k-1;
    rowExchange(B,Q,Qinv,k,i);
  }

  /* ------------------------  ------------------------ */
  template<class matrix, class sqMatrix>
  void rowReduce(matrix& B,sqMatrix& Q,sqMatrix& Qinv,int k,int l){
    int m=B.numberOfRows();
    while( nonZero(MatrixSlice<matrix>(B,k+1,m,l,l)) ){
      rowPrepare(B,Q,Qinv,k,l);
      partRowReduce(B,Q,Qinv,k,l);
    }
  }

  /* ------------------------  ------------------------ */
  template<class matrix, class sqMatrix>
  void rowEchelon(matrix& B,sqMatrix& Q,sqMatrix& Qinv,int &k){
    int m=B.numberOfRows();
    int n=B.numberOfColumns();
    Q.setToIdentity();
    Qinv.setToIdentity();
    k=0;
    int l=1;
    do{
      while(l<=n && !nonZero(MatrixSlice<matrix>(B,k+1,m,l,l))) l++;
      if(l==n+1) break;
      rowReduce(B,Q,Qinv,++k,l);
    }while(k<m);
  }

  // *** column echelon form *** //
  // *** not tested *** //

  /* ------------------------  ------------------------ */

  template<class matrix, class sqMatrix>
  void columnPrepare(matrix& B,sqMatrix& R,sqMatrix& Rinv,int k,int l){
    typedef typename matrix::ScalarType ScalarType;
    int n=B.numberOfColumns();
    ScalarType s;
    int i,j;
    smallestNonZero(MatrixSlice<matrix>(B,k,k,l,n),s,i,j);
    j+=l-1;
    columnExchange(B,R,Rinv,l,j);
  }

  /* ------------------------  ------------------------ */
  template<class matrix, class sqMatrix>
  void columnReduce(matrix& B,sqMatrix& R,sqMatrix& Rinv,int k,int l){
    int n=B.numberOfColumns();
    while( nonZero(MatrixSlice<matrix>(B,k,k,l+1,n)) ){
      columnPrepare(B,R,Rinv,k,l);
      partColumnReduce(B,R,Rinv,k,l);
    }
  }

  /* ------------------------  ------------------------ */
  template<class matrix, class sqMatrix>
  void columnEchelon(matrix& B,sqMatrix& R,sqMatrix& Rinv,int &l){
    int m=B.numberOfRows();
    int n=B.numberOfColumns();
    R.setToIdentity();
    Rinv.setToIdentity();
    l=0;
    int k=1;
    do{
      while(k<=m && !nonZero(MatrixSlice<matrix>(B,k,k,l+1,n)) ) k++;
      if(k==m+1) break;
      columnReduce(B,R,Rinv,k,++l);
    }while(l<n);
  }

}

}

#endif // CAPD_FILE_MATRIXOP_H
