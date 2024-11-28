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

/* CVS $Id: simatrix.hpp,v 1.20 2014/01/30 17:23:48 cxsc Exp $ */

#ifndef _CXSC_SIMATRIX_HPP_INCLUDED
#define _CXSC_SIMATRIX_HPP_INCLUDED

#include <interval.hpp>
#include <imatrix.hpp>
#include <ivector.hpp>
#include <sivector.hpp>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cidot.hpp>
#include <sparseidot.hpp>
#include <sparsematrix.hpp>
#include <srmatrix.hpp>

namespace cxsc {

//definiert in srmatrix.hpp
//enum STORAGE_TYPE{triplet,compressed_row,compressed_column};

class simatrix_slice;
class simatrix_subv;
class scimatrix;
class scimatrix_slice;
class scimatrix_subv;

inline bool comp_pair_i(std::pair<int,interval> p1, std::pair<int,interval> p2) {
  return p1.first < p2.first;
}

//! A sparse interval matrix
/*!
     Sparse matrices in C-XSC are stored in the Compressed Column Storage (CCS) format. The non zero entries of the matrix are stored in three arrays (STL-vectors) \f$ p \f$, \f$ ind \f$, \f$ x \f$. The array \f$ ind \f$ stores the row indices of the non zero elements, the array \f$ x \f$ the respective values (the \f$ i \f$-th element of \f$ ind \f$ corresponds to the i-th element of \f$ x \f$, \f$ i=0,\ldots,nnz-1 \f$, where \f$ nnz \f$ is the number of non zeros. The entries are sorted by column. The array \f$ p \f$ of size \f$ n+1 \f$ stores the starting indices for the entries of each column of the matrix, so that the entries of the \f$ j \f$-th column of the matrix are stored in the elements with index \f$ p[j] \f$ through \f$ p[j+1]-1 \f$ of the arrays \f$ x \f$ and \f$ ind \f$. The elements of each column are stored as sorted by the row indices, explicitly stored zeros are allowed.
 
     The internal data structure uses 0-based indexing throughout. However, in the interface to the user (using the respective operators) every matrix possesses a lower and an upper
     row index bound \f$ lb_1 \f$, \f$ ub_1 \f$ and a lower and an upper column index bound \f$ lb_2 \f$, \f$ ub_2 \f$ of type int. By default, these indexes are 1-based.
 
     It is possible to directly access the internal data structure through the appropriate member function to allow for easier interfacing with other sparse matrix libraries and writing of more efficient sparse matrix algorithms. In this case, the user has to take care that the data structure remains consistent with the format described above. If the user just works with the operators and functions provided by C-XSC, everything will be handled automatically by the C-XSC library.
     
     All matrix and vector operators which require dot product computations use higher precision dot products provided by the dotprecision classes. The precision to be used for these implicit dot products can be choosen by setting the global variable opdotprec accordingly. A value of 0 means maximum accuracy (the default, always used by all older C-XSC versions), a value of 1 means double accuracy, a value of 2 or higher means k-fold double accuracy. Lower accuracy leads to (significantly) faster computing times, but also to less exact results. For all dot products with an interval result, error bounds are computed to guarantee a correct enclosure. For all other dot products approximations without error bounds are computed.
     
     \sa cxsc::dotprecision
*/
class simatrix {

  private:
    std::vector<int> p;
    std::vector<int> ind;
    std::vector<interval> x;
    int m;
    int n;
    int lb1,ub1,lb2,ub2;

  public:

    //! Returns a reference to the vector containing the column pointers (the \f$ p \f$ array)
    std::vector<int>& column_pointers() {
      return p;
    }

    //! Returns a reference to the vector containing the row indices (the \f$ ind \f$ array)
    std::vector<int>& row_indices() {
      return ind;
    }

    //! Returns a reference to the vector containing the stored values (the \f$ x \f$ array)
    std::vector<interval>& values() {
      return x;
    }

    //! Returns a constant reference to the vector containing the column pointers (the \f$ p \f$ array)    
    const std::vector<int>& column_pointers() const {
      return p;
    }

    //! Returns a constant reference to the vector containing the row indices (the \f$ ind \f$ array)
    const std::vector<int>& row_indices() const {
      return ind;
    }

    //! Returns a constant reference to the vector containing the stored values (the \f$ x \f$ array)
    const std::vector<interval>& values() const {
      return x;
    }

    //! Standard constructor, creates an empty matrix of dimension 0x0
    simatrix() {
      p.push_back(0);
      m = n = 0;
      lb1 = lb2 = ub1 = ub2 = 0;
    }

    //! Creates an empty matrix with r rows and c columns, pre-reserving space for 2*(r+c) elements
    simatrix(const int r, const int c) : m(r),n(c),lb1(1),ub1(r),lb2(1),ub2(c) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve(2*(m+n));
      x.reserve(2*(m+n));

      p[0] = 0;
    }

    //! Creates an empty matrix with r rows and c columns, pre-reserving space for e elements
    simatrix(const int r, const int c, const int e) : m(r),n(c),lb1(1),ub1(r),lb2(1),ub2(c) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve(e);
      x.reserve(e);

      p[0] = 0;
    }

    //! Creates a sparse matrix out of three vectors (arrays) forming a matrix stored in triplet, compressed row or compressed column storage
    /*! 
        Creates a sparse matrix of dimension \f$ m \times n \f$ with the values from the three arrays rows, cols and values. These are interpreted according to the storage type t. For triplet
        rows and cols contains the indices of all non zero values and the array values contains the actual values. For compressed_column, the arrays are interpreted according to the compressed column storage format, where cols contains the column pointers. For compressed_row the arrays are interpreted according to the compressed row storage format, which is equivalent to compressed column storage with the role of the rows and columns interchanged.
     
         In each case the arrays need not be sorted by rows or columns and may contain explicit zero entries.
     */
    simatrix(const int m, const int n, const int nnz, const intvector& rows, const intvector& cols, const ivector& values, const enum STORAGE_TYPE t = triplet) {
      if(t == triplet) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         std::vector<triplet_store<interval> > work;
         work.reserve(nnz);

         for(int k=0 ; k<nnz ; k++) {
           work.push_back(triplet_store<interval>(rows[Lb(rows)+k],cols[Lb(cols)+k],values[Lb(values)+k]));
         }

         sort(work.begin(), work.end());

         int i=0;

         for(int j=0 ; j<n ; j++) {        

	   while((unsigned int)i < work.size() && work[i].col == j ) {
               ind.push_back(work[i].row);
               x.push_back(work[i].val);
               i++;
           }

           p[j+1] = i;
         }
         
      } else if(t == compressed_row) {

         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         for(int i=0 ; i<n+1 ; i++)
           p[i] = rows[Lb(rows)+i];

         std::vector<triplet_store<interval> > work;
         work.reserve(nnz);

         for(int j=0 ; j<n ; j++) {
           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(triplet_store<interval>(j,cols[Lb(cols)+k],values[Lb(values)+k]));
           }
         }

         sort(work.begin(), work.end());

         int i=0;

         for(int j=0 ; j<n ; j++) {        

	   while((unsigned int)i < work.size() && work[i].col == j ) {
               ind.push_back(work[i].row);
               x.push_back(work[i].val);
               i++;
           }

           p[j+1] = i;
         }
    
      } else if(t == compressed_column) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         for(int i=0 ; i<n+1 ; i++)
           p[i] = rows[Lb(rows)+i];

         std::vector<std::pair<int,interval> > work;
         work.reserve(n);

         for(int j=0 ; j<n ; j++) {
           work.clear();

           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(std::make_pair(cols[Lb(cols)+k],values[Lb(values)+k]));
           }

           std::sort(work.begin(),work.end(),comp_pair_i);

           for(unsigned int i=0 ; i<work.size() ; i++) {
             ind.push_back(work[i].first);
             x.push_back(work[i].second);
           }
         }

      }

    }

    //! Creates a sparse matrix out of three arrays forming a matrix stored in triplet, compressed row or compressed column storage
    /*!
     Creates a sparse matrix of dimension \f$ m \times n \f$ with the values from the three arrays rows, cols and values. These are interpreted according to the storage type t. For triplet
     rows and cols contains the indices of all non zero values and the array values contains the actual values. For compressed_column, the arrays are interpreted according to the compressed column storage format, where cols contains the column pointers. For compressed_row the arrays are interpreted according to the compressed row storage format, which is equivalent to compressed column storage with the role of the rows and columns interchanged.
     
     In each case the arrays need not be sorted by rows or columns and may contain explicit zero entries.
     
     */
    simatrix(const int m, const int n, const int nnz, const int* rows, const int* cols, const interval* values, const enum STORAGE_TYPE t = triplet) {
      if(t == triplet) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         std::vector<triplet_store<interval> > work;
         work.reserve(nnz);

         for(int k=0 ; k<nnz ; k++) {
           work.push_back(triplet_store<interval>(rows[k],cols[k],values[k]));
         }

         sort(work.begin(), work.end());

         int i=0;

         for(int j=0 ; j<n ; j++) {        

	   while((unsigned int)i < work.size() && work[i].col == j ) {
               ind.push_back(work[i].row);
               x.push_back(work[i].val);
               i++;
           }

           p[j+1] = i;
         }
         
      } else if(t == compressed_row) {

         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         for(int i=0 ; i<n+1 ; i++)
           p[i] = rows[i];

         std::vector<triplet_store<interval> > work;
         work.reserve(nnz);

         for(int j=0 ; j<n ; j++) {
           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(triplet_store<interval>(j,cols[k],values[k]));
           }
         }

         sort(work.begin(), work.end());

         int i=0;

         for(int j=0 ; j<n ; j++) {        

	   while((unsigned int)i < work.size() && work[i].col == j ) {
               ind.push_back(work[i].row);
               x.push_back(work[i].val);
               i++;
           }

           p[j+1] = i;
         }
    
      } else if(t == compressed_column) {
         this->m = m;
         this->n = n;
         p = std::vector<int>(n+1,0);
         ind.reserve(nnz);
         x.reserve(nnz);
         lb1 = lb2 = 1;
         ub1 = m; ub2 = n;

         for(int i=0 ; i<n+1 ; i++)
           p[i] = rows[i];

         std::vector<std::pair<int,interval> > work;
         work.reserve(n);

         for(int j=0 ; j<n ; j++) {
           work.clear();

           for(int k=p[j] ; k<p[j+1] ; k++) {
             work.push_back(std::make_pair(cols[k],values[k]));
           }

           std::sort(work.begin(),work.end(),comp_pair_i);

           for(unsigned int i=0 ; i<work.size() ; i++) {
             ind.push_back(work[i].first);
             x.push_back(work[i].second);
           }
         }

      }

    }


    //! Creates a sparse interval matrix out of a sparse real matrix A.
    simatrix(const srmatrix& A) : p(A.p), ind(A.ind), m(A.m), n(A.n), lb1(A.lb1), ub1(A.ub1), lb2(A.lb2), ub2(A.ub2) {
      x.reserve(A.get_nnz());
      for(unsigned int i=0 ; i<A.x.size() ; i++)
        x.push_back(interval(A.x[i]));
    }


    //! Creates a sparse matrix out of a dense matrix A. Only the non zero elements of A are stored explicitly.
    simatrix(const rmatrix& A) : m(ColLen(A)),n(RowLen(A)),lb1(Lb(A,1)),ub1(Ub(A,1)),lb2(Lb(A,2)),ub2(Ub(A,2)) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);
      x.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);

      p[0] = 0;
      int nnz = 0;

      for(int j=0 ; j<n ; j++) {
        for(int i=0 ; i<m ; i++) {
          if(A[i+lb1][j+lb2] != 0.0) {
             ind.push_back(i);
             x.push_back(interval(A[i+lb1][j+lb2]));
             nnz++;
          }
        }
          
        p[j+1] = nnz;
      }

    }

    //! Creates a sparse matrix out of a dense matrix A. Only the non zero elements of A are stored explicitly.
    simatrix(const imatrix& A) : m(ColLen(A)),n(RowLen(A)),lb1(Lb(A,1)),ub1(Ub(A,1)),lb2(Lb(A,2)),ub2(Ub(A,2)) {
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);
      x.reserve((m*n*0.1 < 2*m) ? (int)(m*n*0.1) : 2*m);

      p[0] = 0;
      int nnz = 0;

      for(int j=0 ; j<n ; j++) {
        for(int i=0 ; i<m ; i++) {
          if(A[i+lb1][j+lb2] != 0.0) {
             ind.push_back(i);
             x.push_back(interval(A[i+lb1][j+lb2]));
             nnz++;
          }
        }
          
        p[j+1] = nnz;
      }

    }

    //! Constructor for banded matrices
    /*!
       Creates a sparse banded matrix of dimension \f$ ms \times ns \f$, whose bands are defined by the columns of the dense matrix A. The column inde range of A must be set appropriately: The columns of with negative indices are used for bands below the diagonal, the column with index 0 is used for the diagonal, and columns with positive index are use for bands above the diagonal.
     */
    simatrix(const int ms, const int ns, const imatrix& A) : m(ms), n(ns), lb1(1), ub1(ms), lb2(1), ub2(ns)  {
      //Banded matrix constructor
      int nnz = RowLen(A)*ColLen(A);
      p = std::vector<int>((n>0) ? n+1 : 1, 0);
      ind.reserve(nnz);
      x.reserve(nnz);

      std::vector<triplet_store<interval> > work;
      work.reserve(nnz);

      
      for(int i=0 ; i<ColLen(A) ; i++) {
        for(int j=Lb(A,2) ; j<=Ub(A,2) ; j++) {
          if(i+j >=0  &&  i+j < n) {
            work.push_back(triplet_store<interval>(i,i+j,A[i+Lb(A,1)][j]));
          }
        }
      }

      sort(work.begin(), work.end());

      int i=0;

      for(int j=0 ; j<n ; j++) {        

        while((unsigned int)i < work.size() && work[i].col == j ) {
          ind.push_back(work[i].row);
          x.push_back(work[i].val);
          i++;
        }

        p[j+1] = i;
      }

    }

    //! Creates a sparse matrix out of a sparse matrix slice
    simatrix(const srmatrix_slice&);
    //! Creates a sparse matrix out of a sparse matrix slice
    simatrix(const simatrix_slice&);

    //! Creates a full matrix out of the sparse matrix and stores it in A. This should normally be done using the respective constructor of the dense matrix.
    void full(imatrix& A) const {
       A = imatrix(lb1,ub1,lb2,ub2);
       A = 0.0;
       for(int j=0 ; j<n ; j++) {
          for(int k=p[j] ; k<p[j+1] ; k++) {
             A[ind[k]+lb1][j+lb2] = x[k];
          }
       }
    }

    //! Drops explicitly stored zeros from the data structure.
    /*!
        Dropping explicit zero entries can be essential when interfacing with other sparse matrix libraries. Some librares using the CCS storage format do
        not allow explicitly stored zeros.
     */
    void dropzeros() {
      std::vector<int> pnew(n+1,0);
      std::vector<int> indnew;
      std::vector<interval> xnew;
      int nnznew = 0;

      for(int j=0 ; j<n ; j++) {
        for(int k=p[j] ; k<p[j+1] ; k++) {
          if(x[k] != 0.0) {
            xnew.push_back(x[k]);
            indnew.push_back(ind[k]);
            nnznew++;
          }
        }
        pnew[j+1] = nnznew;
      }

      p = pnew;
      ind = indnew;
      x = xnew;
    }


    //! Assigns a real value to all elements of the matrix (resulting in a dense matrix!)
    simatrix& operator=(const real& A) {
      return sp_ms_assign<simatrix,real,interval>(*this,A);
    }

    //! Assigns an interval value to all elements of the matrix (resulting in a dense matrix!)
    simatrix& operator=(const interval& A) {
      return sp_ms_assign<simatrix,interval,interval>(*this,A);
    }

    //! Assigns a dense matrix to the sparse matrix. Only the non zero entries of the dense matrix are used.
    simatrix& operator=(const rmatrix& A) {
      return spf_mm_assign<simatrix,rmatrix,interval>(*this,A);
    }

    //! Assigns a dense matrix to the sparse matrix. Only the non zero entries of the dense matrix are used.
    simatrix& operator=(const imatrix& A) {
      return spf_mm_assign<simatrix,imatrix,interval>(*this,A);
    }

    //! Assigns a dense matrix slice to the sparse matrix. Only the non zero entries of the dense matrix are used.
    simatrix& operator=(const rmatrix_slice& A) {
      return spf_mm_assign<simatrix,rmatrix_slice,interval>(*this,A);
    }

    //! Assigns a dense matrix slice to the sparse matrix. Only the non zero entries of the dense matrix are used.
    simatrix& operator=(const imatrix_slice& A) {
      return spf_mm_assign<simatrix,imatrix_slice,interval>(*this,A);
    }

    //! Assign a sparse real to a sparse interval matrix
    simatrix& operator=(const srmatrix& A) {
      m = A.m;
      n = A.n;
      p = A.p;
      ind = A.ind;
      x.clear();
      x.reserve(A.get_nnz());
      for(unsigned int i=0 ; i<A.x.size() ; i++)
        x.push_back(interval(A.x[i]));
      return *this;
    }

    /* simatrix& operator=(const simatrix& A) {
      p = A.p;
      ind = A.ind;
      x = A.x;
      return *this;
    } */

    //! Assign a sparse matrix slice to a sparse matrix
    simatrix& operator=(const srmatrix_slice&);
    //! Assign a sparse matrix slice to a sparse matrix
    simatrix& operator=(const simatrix_slice&);

    //! Returns a copy of the element in row i and column j
    /*!
       This operator can be used for read access only. The indices i and j must be used according to the current index range of the matrix. A copy of the element (i,j) is returned, or 0 if this element
       is not explicitly stored. For write access to a single element, the []-opeator or the member function element should be used.
     
       Note that due to the underlying data structure the access to single elements of a sparse matrix is much more expensive than to the elements of a dense matrix.
     */
    const interval operator()(int i, int j) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb1 || i>ub1 || j<lb2 || j>ub2)
        cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator()(int, int)"));
#endif
      interval r(0.0);
      for(int k=p[j-lb2] ; k<p[j-lb2+1] && ind[k]<=i-lb1 ; k++) {
        if(ind[k] == i-lb1)  r = x[k];
      }
      return r;
    }

    //! Returns a reference to the element (i,j) of the matrix
    /*!
     This function should only be used for write access. The indices i and j must be used according to the current index range of the matrix. A reference to the element (i,j) is returned. If the element is not stored explicitly, it is created as an explicit 0 entry.
     For read access to a single element, the ()-opeator should be used.
     
     This function is faster than using the []-operator (for example A[i][j]), since it does not need to create a temporary subvector first.
     
     Note that due to the underlying data structure the access to single elements of a sparse matrix is much more expensive than to the elements of a dense matrix.
     */
    interval& element(int i, int j) {
#if(CXSC_INDEX_CHECK)
      if(i<lb1 || i>ub1 || j<lb2 || j>ub2)
        cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator()(int, int)"));
#endif
      int k;
      for(k=p[j-lb2] ; k<p[j-lb2+1] && ind[k]<=i-lb1 ; k++) {
        if(ind[k] == i-lb1)  return x[k];
      }

      //Nicht gefunden, Element muss angelegt werden, da Schreibzugriff moeglich
      std::vector<int>::iterator ind_it = ind.begin() + k;
      std::vector<interval>::iterator x_it  = x.begin() + k;
      ind.insert(ind_it, i-lb1);
      x_it = x.insert(x_it, interval(0.0));
      for(k=j-lb2+1 ; k<(int)p.size() ; k++)
        p[k]++;

      return *x_it;
    }

    //! Returns a column of the matrix as a sparse subvector object
    simatrix_subv operator[](const cxscmatrix_column&);
    //! Returns a row of the matrix as a sparse subvector object
    simatrix_subv operator[](const int);
    //! Returns a column of the matrix as a sparse subvector object
    const simatrix_subv operator[](const cxscmatrix_column&) const;
    //! Returns a row of the matrix as a sparse subvector object
    const simatrix_subv operator[](const int) const;

    //! Returns a slice of the matrix
    simatrix_slice operator()(const int, const int , const int, const int);
    //! Returns a slice of the matrix
    const simatrix_slice operator()(const int, const int , const int, const int) const;

    //! Performs a row and column permutation using two permutation vectors
    simatrix operator()(const intvector& pervec, const intvector& q) {
      simatrix A(m,n,get_nnz());
      intvector per = perminv(pervec);

      int nnz=0;
      for(int k=0 ; k<n ; k++) {
        A.p[k] = nnz;

        std::map<int,interval> work;
        for(int j=p[q[Lb(q)+k]] ; j<p[q[Lb(q)+k]+1] ; j++) 
           work.insert(std::make_pair(per[Lb(per)+ind[j]], x[j]));
        
        for(std::map<int,interval>::iterator it = work.begin() ; it != work.end() ; it++) {
           A.ind.push_back(it->first);
           A.x.push_back(it->second);
        }

        nnz += work.size();
 
      }

      A.p[n] = nnz;

      return A;
    }

    //! Performs a row permutation using a permutation vector
    simatrix operator()(const intvector& pervec) {
      simatrix A(m,n,get_nnz());
      intvector per = perminv(pervec);

      for(int k=0 ; k<n ; k++) {
        A.p[k] = p[k];

        std::map<int,interval> work;
        for(int j=p[k] ; j<p[k+1] ; j++) 
           work.insert(std::make_pair(per[Lb(per)+ind[j]], x[j]));
        
        for(std::map<int,interval>::iterator it = work.begin() ; it != work.end() ; it++) {
           A.ind.push_back(it->first);
           A.x.push_back(it->second);
        }
 
      }

      A.p[n] = p[n];

      return A;
    }

    //! Performs row and column permutations using the two permutation matrices P and Q. Faster than explicitly computing the product.
    simatrix operator()(const intmatrix& P, const intmatrix& Q) {
      intvector p = permvec(P);
      intvector q = perminv(permvec(Q));
      return (*this)(p,q);
    }

    //! Performs a row permutation using the permutation matrix P. Faster than explicitly computing the product.
    simatrix operator()(const intmatrix& P) {
      intvector p = permvec(P);
      return (*this)(p);
    }

    //! Returns the density (the number of non-zeros divided by the number of elements) of the matrix
    real density() const {
      return p[n]/((double)m*n);
    }

    //! Returns the number of non-zero entries (including explicitly stored zeros).
    int get_nnz() const {
      return p[n];
    }

    //! Add B to the sparse matrix and assign the result to it.
    simatrix& operator+=(const rmatrix& B) {
      return spf_mm_addassign<simatrix,rmatrix,imatrix>(*this,B);
    }

    //! Add B to the sparse matrix and assign the result to it.
    simatrix& operator+=(const imatrix& B) {
      return spf_mm_addassign<simatrix,imatrix,imatrix>(*this,B);
    }

    //! Add B to the sparse matrix and assign the result to it.
    simatrix& operator+=(const rmatrix_slice& B) {
      return spf_mm_addassign<simatrix,rmatrix_slice,imatrix>(*this,B);
    }

    //! Add B to the sparse matrix and assign the result to it.
    simatrix& operator+=(const imatrix_slice& B) {
      return spf_mm_addassign<simatrix,imatrix_slice,imatrix>(*this,B);
    }

    //! Add B to the sparse matrix and assign the result to it.
    simatrix& operator+=(const srmatrix& B) {
      return spsp_mm_addassign<simatrix,srmatrix,interval>(*this,B);
    }

    //! Add B to the sparse matrix and assign the result to it.
    simatrix& operator+=(const simatrix& B) {
      return spsp_mm_addassign<simatrix,simatrix,interval>(*this,B);
    }

    //! Subtract B from the sparse matrix and assign the result to it.
    simatrix& operator-=(const rmatrix& B) {
      return spf_mm_subassign<simatrix,rmatrix,imatrix>(*this,B);
    }

    //! Subtract B from the sparse matrix and assign the result to it.
    simatrix& operator-=(const imatrix& B) {
      return spf_mm_subassign<simatrix,imatrix,imatrix>(*this,B);
    }

    //! Subtract B from the sparse matrix and assign the result to it.
    simatrix& operator-=(const rmatrix_slice& B) {
      return spf_mm_subassign<simatrix,rmatrix_slice,imatrix>(*this,B);
    }

    //! Subtract B from the sparse matrix and assign the result to it.
    simatrix& operator-=(const imatrix_slice& B) {
      return spf_mm_subassign<simatrix,imatrix_slice,imatrix>(*this,B);
    }

    //! Subtract B from the sparse matrix and assign the result to it.
    simatrix& operator-=(const srmatrix& B) {
      return spsp_mm_subassign<simatrix,srmatrix,interval>(*this,B);
    }

    //! Subtract B from the sparse matrix and assign the result to it.
    simatrix& operator-=(const simatrix& B) {
      return spsp_mm_subassign<simatrix,simatrix,interval>(*this,B);
    }

    //! Multiply the sparse matrix by B and assign the result to it.
    simatrix& operator*=(const imatrix& B) {
      return spf_mm_multassign<simatrix,imatrix,sparse_idot,imatrix>(*this,B);
    }

    //! Multiply the sparse matrix by B and assign the result to it.
    simatrix& operator*=(const rmatrix& B) {
      return spf_mm_multassign<simatrix,rmatrix,sparse_idot,imatrix>(*this,B);
    }

    //! Multiply the sparse matrix by B and assign the result to it.
    simatrix& operator*=(const rmatrix_slice& B) {
      return spf_mm_multassign<simatrix,rmatrix_slice,sparse_idot,imatrix>(*this,B);
    }

    //! Multiply the sparse matrix by B and assign the result to it.
    simatrix& operator*=(const imatrix_slice& B) {
      return spf_mm_multassign<simatrix,imatrix_slice,sparse_idot,imatrix>(*this,B);
    }

    //! Multiply the sparse matrix by B and assign the result to it.
    simatrix& operator*=(const srmatrix& B) {
      return spsp_mm_multassign<simatrix,srmatrix,sparse_idot,interval>(*this,B);
    }

    //! Multiply the sparse matrix by B and assign the result to it.
    simatrix& operator*=(const simatrix& B) {
      return spsp_mm_multassign<simatrix,simatrix,sparse_idot,interval>(*this,B);
    }

    //! Multiply all elements of the sparse matrix by r and assign the result to it.
    simatrix& operator*=(const real& r) {
      return sp_ms_multassign(*this,r);
    }

    //! Multiply all elements of the sparse matrix by r and assign the result to it.
    simatrix& operator*=(const interval& r) {
      return sp_ms_multassign(*this,r);
    }

    //! Divide all elements of the sparse matrix by r and assign the result to it.
    simatrix& operator/=(const real& r) {
      return sp_ms_divassign(*this,r);
    }

    //! Divide all elements of the sparse matrix by r and assign the result to it.
    simatrix& operator/=(const interval& r) {
      return sp_ms_divassign(*this,r);
    }

    //! Form the convex hull of a sparse matrix and B and assign the result to it.
    simatrix& operator|=(const rmatrix& B) {
      return spf_mm_hullassign<simatrix,rmatrix,imatrix>(*this,B);
    }

    //! Form the convex hull of a sparse matrix and B and assign the result to it.
    simatrix& operator|=(const imatrix& B) {
      return spf_mm_hullassign<simatrix,imatrix,imatrix>(*this,B);
    }

    //! Form the convex hull of a sparse matrix and B and assign the result to it.
    simatrix& operator|=(const rmatrix_slice& B) {
      return spf_mm_hullassign<simatrix,rmatrix_slice,imatrix>(*this,B);
    }

    //! Form the convex hull of a sparse matrix and B and assign the result to it.
    simatrix& operator|=(const imatrix_slice& B) {
      return spf_mm_hullassign<simatrix,imatrix_slice,imatrix>(*this,B);
    }

    //! Form the convex hull of a sparse matrix and B and assign the result to it.
    simatrix& operator|=(const srmatrix& B) {
      return spsp_mm_hullassign<simatrix,srmatrix,interval>(*this,B);
    }

    //! Form the convex hull of a sparse matrix and B and assign the result to it.
    simatrix& operator|=(const simatrix& B) {
      return spsp_mm_hullassign<simatrix,simatrix,interval>(*this,B);
    }

    //! Form the intersection of a sparse matrix and B and assign the result to it.
    simatrix& operator&=(const imatrix& B) {
      return spf_mm_intersectassign<simatrix,imatrix,imatrix>(*this,B);
    }

    //! Form the convex hull of a sparse matrix and B and assign the result to it.
    simatrix& operator&=(const imatrix_slice& B) {
      return spf_mm_intersectassign<simatrix,imatrix_slice,imatrix>(*this,B);
    }

    //! Form the convex hull of a sparse matrix and B and assign the result to it.
    simatrix& operator&=(const simatrix& B) {
      return spsp_mm_intersectassign<simatrix,simatrix,interval>(*this,B);
    }

    friend void SetLb(simatrix&, const int, const int);
    friend void SetUb(simatrix&, const int, const int);    
    friend int Lb(const simatrix&, int);
    friend int Ub(const simatrix&, int);
    friend int RowLen(const simatrix&);
    friend int ColLen(const simatrix&);
    friend srmatrix Inf(const simatrix&);
    friend srmatrix Sup(const simatrix&);
    friend simatrix Re(const scimatrix&);
    friend simatrix Im(const scimatrix&);
    friend simatrix abs(const simatrix&);
    friend srmatrix mid(const simatrix&);
    friend srmatrix diam(const simatrix&);
    friend simatrix abs(const scimatrix&);
    friend srmatrix absmin(const simatrix&);
    friend srmatrix absmax(const simatrix&);

    friend simatrix transp(const simatrix&);
    friend simatrix Id(const simatrix&);
    friend srmatrix CompMat(const simatrix&);

    friend std::istream& operator>>(std::istream&, simatrix_slice&);
    friend std::istream& operator>>(std::istream&, simatrix_subv&);

    friend class srmatrix_slice;
    friend class srmatrix_subv;
    friend class srvector;
    friend class simatrix_slice;
    friend class simatrix_subv;
    friend class sivector;
    friend class scimatrix;
    friend class scimatrix_slice;
    friend class scimatrix_subv;
    friend class scivector;
    friend class rmatrix;
    friend class imatrix;
    friend class cimatrix;

#include "matrix_friend_declarations.inl"
};

inline imatrix::imatrix(const srmatrix& A) {
  dat = new interval[A.m*A.n];
  lb1 = A.lb1; lb2 = A.lb2; ub1 = A.ub1; ub2 = A.ub2;
  xsize = A.n;
  ysize = A.m;
  *this = 0.0;
  for(int j=0 ; j<A.n ; j++) {
     for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
        dat[A.ind[k]*A.n+j] = A.x[k];
     }
  }
}

inline imatrix::imatrix(const simatrix& A) {
  dat = new interval[A.m*A.n];
  lb1 = A.lb1; lb2 = A.lb2; ub1 = A.ub1; ub2 = A.ub2;
  xsize = A.n;
  ysize = A.m;
  *this = 0.0;
  for(int j=0 ; j<A.n ; j++) {
     for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
        dat[A.ind[k]*A.n+j] = A.x[k];
     }
  }
}

//! Return a sparse unity matrix of the same dimension as A
inline simatrix Id(const simatrix& A) {
  simatrix I(A.m, A.n, (A.m>A.n) ? A.m : A.n);
  I.lb1 = A.lb1; I.lb2 = A.lb2;
  I.ub1 = A.ub1; I.ub2 = A.ub2;

  if(A.m < A.n) {
    for(int i=0 ; i<A.m ; i++) {
      I.p[i+1] = I.p[i] + 1;
      I.ind.push_back(i);
      I.x.push_back(interval(1.0));
    }
  } else {
    for(int i=0 ; i<A.n ; i++) {
      I.p[i+1] = I.p[i] + 1;
      I.ind.push_back(i);
      I.x.push_back(interval(1.0));
    }
  }

  return I;
}

//! Returns the transpose of A
inline simatrix transp(const simatrix& A) {
  simatrix B(A.n, A.m, A.get_nnz());
     
  //NIchtnullen pro Zeile bestimmen
  std::vector<int> w(A.m,0);
  for(unsigned int i=0 ; i<A.ind.size() ; i++) 
    w[A.ind[i]]++;

  //Spalten"pointer" setzen
  B.p.resize(A.m+1);
  B.p[0] = 0;
  for(unsigned int i=1 ; i<B.p.size() ; i++)
    B.p[i] = w[i-1] + B.p[i-1];

  //w vorbereiten
  w.insert(w.begin(), 0); 
  for(unsigned int i=1 ; i<w.size() ; i++) {
    w[i] += w[i-1];
  }

  //neuer zeilenindex und wert wird gesetzt
  int q;
  B.ind.resize(A.get_nnz());
  B.x.resize(A.get_nnz());
  for(int j=0 ; j<A.n ; j++) {
    for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
      q = w[A.ind[k]]++;
      B.ind[q] = j;
      B.x[q] = A.x[k];
    }
  }

  return B;
}

//! Sets the lower index bound of the rows (i==ROW) or columns (i==COL) to j.
/*!
    If i==ROW, the lower index bound for the rows of A ist set to j. if i==COL, the lower index
    bound of the columns is set to j. The upper index bound is automatically set according to the number of rows or columns of the matrix.
 */
inline void SetLb(simatrix& A, const int i, const int j) {
  if(i==1) {
    A.lb1 = j;
    A.ub1 = j + A.m - 1;
  } else if(i==2) {
    A.lb2 = j;
    A.ub2 = j + A.n - 1;
  }
}

//! Sets the upper index bound of the rows (i==ROW) or columns (i==COL) to j.
/*!
If i==ROW, the upper index bound for the rows of A ist set to j. if i==COL, the upper index
bound of the columns is set to j. The lower index bound is automatically set according to the number of rows or columns of the matrix.
*/
inline void SetUb(simatrix& A, const int i, const int j) {
  if(i==1) {
    A.ub1 = j;
    A.lb1 = j - A.m + 1;
  } else if(i==2) {
    A.ub2 = j;
    A.lb2 = j - A.n + 1;
  }
}

//! Returns the lower index bound for the rows or columns of A
/*! 
   If i==ROW, the lower index bound for the rows is returned, if i==COL, the lower index bound for the columns is returned.
 */
inline int Lb(const simatrix& A, int i) {
  if(i==1) 
    return A.lb1;
  else if(i==2)
    return A.lb2;
  else
    return 1;
}

//! Returns the upper index bound for the rows or columns of A
/*!
  If i==COL, the upper index bound for the rows is returned, if i==COL, the upper index bound for the columns is returned.
*/
inline int Ub(const simatrix& A, int i) {
  if(i==1) 
    return A.ub1;
  else if(i==2)
    return A.ub2;
  else
    return 1;
}

//! Returns the number of columns of the matrix
inline int RowLen(const simatrix& A) {
  return A.n;
}

//! Returns the number of rows of the matrix
inline int ColLen(const simatrix& A) {
  return A.m;
}

//! Resizes the matrix to a \f$ 0 \times 0 \f$ matrix
inline void Resize(simatrix& A) {
  sp_m_resize(A);
}

//! Resizes the matrix to a \f$ m \times n \f$ matrix, preserving as many of the old entries as possible.
inline void Resize(simatrix& A, const int m, const int n) {
  sp_m_resize(A,m,n);
}

//! Resizes the matrix to u1-l1+1 rows and u2-l2+1 columns, preserving as many of the old entries as possible and setting the index range accordingly.
inline void Resize(simatrix& A, const int l1, const int u1, const int l2, const int u2) {
  sp_m_resize(A,l1,u1,l2,u2);
}

//! Returns the Infimum of the matrix A
inline srmatrix Inf(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(Inf(A.x[i]));

  res.dropzeros();

  return res; 
}

//! Returns the Supremum of the matrix A
inline srmatrix Sup(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(Sup(A.x[i]));

  res.dropzeros();

  return res; 
}

//! Returns the componentwise absolute value as the interval hull of \f$ \{a_{ij} \mid a_{ij} \in [a_{ij}] \} \f$
inline simatrix abs(const simatrix& A) {
  simatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(abs(A.x[i]));

  res.dropzeros();

  return res; 
}

//! Returns the componentwise minimum absolute value 
inline srmatrix absmin(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(AbsMin(A.x[i]));

  res.dropzeros();

  return res; 
}

//! Returns the componentwise maximum absolute value 
inline srmatrix absmax(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(AbsMax(A.x[i]));

  res.dropzeros();

  return res; 
}

//! Returns Ostroswkis comparison matrix for A
inline srmatrix CompMat(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;
  res.p[A.n] = A.p[A.n];

  for(int j=0 ; j<res.n ; j++) {
    for(int k=A.p[j] ; k<A.p[j+1] ; k++) {
      if(A.ind[k] == j)
        res.x.push_back(AbsMin(A.x[k]));
      else
        res.x.push_back(-AbsMax(A.x[k]));
    }
  }

  res.dropzeros();

  return res; 
}

//! Returns the midpoint matrix for A
inline srmatrix mid(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(mid(A.x[i]));

  res.dropzeros();

  return res; 
}

//! Returns the componentwise diameter of A
inline srmatrix diam(const simatrix& A) {
  srmatrix res(A.m,A.n,A.get_nnz());
  res.lb1 = A.lb1;
  res.lb2 = A.lb2;
  res.ub1 = A.ub1;
  res.ub2 = A.ub2;
  res.p   = A.p;
  res.ind = A.ind;

  for(int i=0 ; i<res.get_nnz() ; i++)
    res.x.push_back(diam(A.x[i]));

  res.dropzeros();

  return res; 
}


//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const imatrix& A, const srmatrix& B) {
  return fsp_mm_mult<imatrix,srmatrix,imatrix,sparse_idot>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const rmatrix& A, const simatrix& B) {
  return fsp_mm_mult<rmatrix,simatrix,imatrix,sparse_idot>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const imatrix& A, const simatrix& B) {
  return fsp_mm_mult<imatrix,simatrix,imatrix,sparse_idot>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const simatrix& A, const rmatrix& B) {
  return spf_mm_mult<simatrix,rmatrix,imatrix,sparse_idot>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const srmatrix& A, const imatrix& B) {
  return spf_mm_mult<srmatrix,imatrix,imatrix,sparse_idot>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const simatrix& A, const imatrix& B) {
  return spf_mm_mult<simatrix,imatrix,imatrix,sparse_idot>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_mult<imatrix_slice,srmatrix,imatrix,sparse_idot>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_mult<rmatrix_slice,simatrix,imatrix,sparse_idot>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_mult<imatrix_slice,simatrix,imatrix,sparse_idot>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_mult<simatrix,rmatrix_slice,imatrix,sparse_idot>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_mult<srmatrix,imatrix_slice,imatrix,sparse_idot>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_mult<simatrix,imatrix_slice,imatrix,sparse_idot>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline simatrix operator*(const simatrix& A, const srmatrix& B) {
  return spsp_mm_mult<simatrix,srmatrix,simatrix,sparse_idot,interval>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline simatrix operator*(const srmatrix& A, const simatrix& B) {
  return spsp_mm_mult<srmatrix,simatrix,simatrix,sparse_idot,interval>(A,B);
}

//! Returns the product of the matrices A and B.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline simatrix operator*(const simatrix& A, const simatrix& B) {
  return spsp_mm_mult<simatrix,simatrix,simatrix,sparse_idot,interval>(A,B);
}

//! Divides every element of A by r and returns the result
inline simatrix operator/(const simatrix& A, const real& r) {
  return sp_ms_div<simatrix,real,simatrix>(A,r);
}

//! Divides every element of A by r and returns the result
inline simatrix operator/(const simatrix& A, const interval& r) {
  return sp_ms_div<simatrix,interval,simatrix>(A,r);
}

//! Divides every element of A by r and returns the result
inline simatrix operator/(const srmatrix& A, const interval& r) {
  return sp_ms_div<srmatrix,interval,simatrix>(A,r);
}

//! Multiplies every element of A by r and returns the result
inline simatrix operator*(const simatrix& A, const real& r) {
  return sp_ms_mult<simatrix,real,simatrix>(A,r);
}

//! Multiplies every element of A by r and returns the result
inline simatrix operator*(const simatrix& A, const interval& r) {
  return sp_ms_mult<simatrix,interval,simatrix>(A,r);
}

//! Multiplies every element of A by r and returns the result
inline simatrix operator*(const srmatrix& A, const interval& r) {
  return sp_ms_mult<srmatrix,interval,simatrix>(A,r);
}

//! Multiplies every element of A by r and returns the result
inline simatrix operator*(const real& r, const simatrix& A) {
  return sp_sm_mult<real,simatrix,simatrix>(r,A);
}

//! Multiplies every element of A by r and returns the result
inline simatrix operator*(const interval& r, const simatrix& A) {
  return sp_sm_mult<interval,simatrix,simatrix>(r,A);
}

//! Multiplies every element of A by r and returns the result
inline simatrix operator*(const interval& r, const srmatrix& A) {
  return sp_sm_mult<interval,srmatrix,simatrix>(r,A);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const simatrix& A, const rvector& v) {
  return spf_mv_mult<simatrix,rvector,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const srmatrix& A, const ivector& v) {
  return spf_mv_mult<srmatrix,ivector,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const simatrix& A, const ivector& v) {
  return spf_mv_mult<simatrix,ivector,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const simatrix& A, const rvector_slice& v) {
  return spf_mv_mult<simatrix,rvector_slice,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const srmatrix& A, const ivector_slice& v) {
  return spf_mv_mult<srmatrix,ivector_slice,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const simatrix& A, const ivector_slice& v) {
  return spf_mv_mult<simatrix,ivector_slice,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline sivector operator*(const simatrix& A, const srvector& v) {
  return spsp_mv_mult<simatrix,srvector,sivector,sparse_idot,interval>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline sivector operator*(const srmatrix& A, const sivector& v) {
  return spsp_mv_mult<srmatrix,sivector,sivector,sparse_idot,interval>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline sivector operator*(const simatrix& A, const sivector& v) {
  return spsp_mv_mult<simatrix,sivector,sivector,sparse_idot,interval>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline sivector operator*(const simatrix& A, const srvector_slice& v) {
  return spsl_mv_mult<simatrix,srvector_slice,sivector,sparse_idot,interval>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline sivector operator*(const srmatrix& A, const sivector_slice& v) {
  return spsl_mv_mult<srmatrix,sivector_slice,sivector,sparse_idot,interval>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline sivector operator*(const simatrix& A, const sivector_slice& v) {
  return spsl_mv_mult<simatrix,sivector_slice,sivector,sparse_idot,interval>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const imatrix& A, const srvector& v) {
  return fsp_mv_mult<imatrix,srvector,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const rmatrix& A, const sivector& v) {
  return fsp_mv_mult<rmatrix,sivector,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const imatrix& A, const sivector& v) {
  return fsp_mv_mult<imatrix,sivector,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const imatrix_slice& A, const srvector& v) {
  return fsp_mv_mult<imatrix_slice,srvector,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const rmatrix_slice& A, const sivector& v) {
  return fsp_mv_mult<rmatrix_slice,sivector,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const imatrix_slice& A, const sivector& v) {
  return fsp_mv_mult<imatrix_slice,sivector,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const imatrix& A, const srvector_slice& v) {
  return fsl_mv_mult<imatrix,srvector_slice,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const rmatrix& A, const sivector_slice& v) {
  return fsl_mv_mult<rmatrix,sivector_slice,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const imatrix& A, const sivector_slice& v) {
  return fsl_mv_mult<imatrix,sivector_slice,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const imatrix_slice& A, const srvector_slice& v) {
  return fsl_mv_mult<imatrix_slice,srvector_slice,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const rmatrix_slice& A, const sivector_slice& v) {
  return fsl_mv_mult<rmatrix_slice,sivector_slice,ivector,sparse_idot>(A,v);
}

//! Returns the product of the matrix A and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const imatrix_slice& A, const sivector_slice& v) {
  return fsl_mv_mult<imatrix_slice,sivector_slice,ivector,sparse_idot>(A,v);
}

//! Returns the elementwise sum of the matrices A and B.
inline imatrix operator+(const imatrix& A, const srmatrix& B) {
  return fsp_mm_add<imatrix,srmatrix,imatrix>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline imatrix operator+(const rmatrix& A, const simatrix& B) {
  return fsp_mm_add<rmatrix,simatrix,imatrix>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline imatrix operator+(const imatrix& A, const simatrix& B) {
  return fsp_mm_add<imatrix,simatrix,imatrix>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline imatrix operator+(const simatrix& A, const rmatrix& B) {
  return spf_mm_add<simatrix,rmatrix,imatrix>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline imatrix operator+(const srmatrix& A, const imatrix& B) {
  return spf_mm_add<srmatrix,imatrix,imatrix>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline imatrix operator+(const simatrix& A, const imatrix& B) {
  return spf_mm_add<simatrix,imatrix,imatrix>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline imatrix operator+(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_add<imatrix_slice,srmatrix,imatrix>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline imatrix operator+(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_add<rmatrix_slice,simatrix,imatrix>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline imatrix operator+(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_add<imatrix_slice,simatrix,imatrix>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline imatrix operator+(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_add<simatrix,rmatrix_slice,imatrix>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline imatrix operator+(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_add<srmatrix,imatrix_slice,imatrix>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline imatrix operator+(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_add<simatrix,imatrix_slice,imatrix>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline simatrix operator+(const simatrix& A, const srmatrix& B) {
  return spsp_mm_add<simatrix,srmatrix,simatrix,interval>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline simatrix operator+(const srmatrix& A, const simatrix& B) {
  return spsp_mm_add<srmatrix,simatrix,simatrix,interval>(A,B);
}

//! Returns the elementwise sum of the matrices A and B.
inline simatrix operator+(const simatrix& A, const simatrix& B) {
  return spsp_mm_add<simatrix,simatrix,simatrix,interval>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline imatrix operator-(const imatrix& A, const srmatrix& B) {
  return fsp_mm_sub<imatrix,srmatrix,imatrix>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline imatrix operator-(const rmatrix& A, const simatrix& B) {
  return fsp_mm_sub<rmatrix,simatrix,imatrix>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline imatrix operator-(const imatrix& A, const simatrix& B) {
  return fsp_mm_sub<imatrix,simatrix,imatrix>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline imatrix operator-(const simatrix& A, const rmatrix& B) {
  return spf_mm_sub<simatrix,rmatrix,imatrix>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline imatrix operator-(const srmatrix& A, const imatrix& B) {
  return spf_mm_sub<srmatrix,imatrix,imatrix>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline imatrix operator-(const simatrix& A, const imatrix& B) {
  return spf_mm_sub<simatrix,imatrix,imatrix>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline imatrix operator-(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_sub<imatrix_slice,srmatrix,imatrix>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline imatrix operator-(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_sub<rmatrix_slice,simatrix,imatrix>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline imatrix operator-(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_sub<imatrix_slice,simatrix,imatrix>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline imatrix operator-(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_sub<simatrix,rmatrix_slice,imatrix>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline imatrix operator-(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_sub<srmatrix,imatrix_slice,imatrix>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline imatrix operator-(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_sub<simatrix,imatrix_slice,imatrix>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline simatrix operator-(const simatrix& A, const srmatrix& B) {
  return spsp_mm_sub<simatrix,srmatrix,simatrix,interval>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline simatrix operator-(const srmatrix& A, const simatrix& B) {
  return spsp_mm_sub<srmatrix,simatrix,simatrix,interval>(A,B);
}

//! Returns the elementwise difference of the matrices A and B.
inline simatrix operator-(const simatrix& A, const simatrix& B) {
  return spsp_mm_sub<simatrix,simatrix,simatrix,interval>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const imatrix& A, const srmatrix& B) {
  return fsp_mm_hull<imatrix,srmatrix,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const rmatrix& A, const simatrix& B) {
  return fsp_mm_hull<rmatrix,simatrix,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const imatrix& A, const simatrix& B) {
  return fsp_mm_hull<imatrix,simatrix,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const simatrix& A, const rmatrix& B) {
  return spf_mm_hull<simatrix,rmatrix,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const srmatrix& A, const imatrix& B) {
  return spf_mm_hull<srmatrix,imatrix,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const simatrix& A, const imatrix& B) {
  return spf_mm_hull<simatrix,imatrix,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_hull<imatrix_slice,srmatrix,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_hull<rmatrix_slice,simatrix,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_hull<imatrix_slice,simatrix,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_hull<simatrix,rmatrix_slice,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_hull<srmatrix,imatrix_slice,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_hull<simatrix,imatrix_slice,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline simatrix operator|(const simatrix& A, const srmatrix& B) {
  return spsp_mm_hull<simatrix,srmatrix,simatrix,interval>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline simatrix operator|(const srmatrix& A, const simatrix& B) {
  return spsp_mm_hull<srmatrix,simatrix,simatrix,interval>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline simatrix operator|(const simatrix& A, const simatrix& B) {
  return spsp_mm_hull<simatrix,simatrix,simatrix,interval>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const rmatrix& A, const srmatrix& B) {
  return fsp_mm_hull<rmatrix,srmatrix,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const srmatrix& A, const rmatrix& B) {
  return spf_mm_hull<srmatrix,rmatrix,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const rmatrix_slice& A, const srmatrix& B) {
  return fsp_mm_hull<rmatrix_slice,srmatrix,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline imatrix operator|(const srmatrix& A, const rmatrix_slice& B) {
  return spf_mm_hull<srmatrix,rmatrix_slice,imatrix>(A,B);
}

//! Returns the elementwise convex hull of the matrices A and B.
inline simatrix operator|(const srmatrix& A, const srmatrix& B) {
  return spsp_mm_hull<srmatrix,srmatrix,simatrix,interval>(A,B);
}

//! Returns the elementwise intersection of the matrices A and B.
inline imatrix operator&(const imatrix& A, const simatrix& B) {
  return fsp_mm_intersect<imatrix,simatrix,imatrix>(A,B);
}

//! Returns the elementwise intersection of the matrices A and B.
inline imatrix operator&(const simatrix& A, const imatrix& B) {
  return spf_mm_intersect<simatrix,imatrix,imatrix>(A,B);
}

//! Returns the elementwise intersection of the matrices A and B.
inline imatrix operator&(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_intersect<imatrix_slice,simatrix,imatrix>(A,B);
}

//! Returns the elementwise intersection of the matrices A and B.
inline imatrix operator&(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_intersect<simatrix,imatrix_slice,imatrix>(A,B);
}

//! Returns the elementwise intersection of the matrices A and B.
inline simatrix operator&(const simatrix& A, const simatrix& B) {
  return spsp_mm_intersect<simatrix,simatrix,simatrix,interval>(A,B);
}

//! Unary component-wise negation of M
inline simatrix operator-(const simatrix& M) {
  return sp_m_negative<simatrix,simatrix>(M);
}

//! Unary component-wise operator +
inline simatrix& operator+(simatrix& A) {
  return A;
}

inline imatrix& imatrix::operator=(const srmatrix& B) {
  *this = rmatrix(B);
  return *this;
}

inline imatrix& imatrix::operator=(const simatrix& B) {
  *this = imatrix(B);
  return *this;
}

inline imatrix& imatrix::operator+=(const srmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline imatrix& imatrix::operator+=(const simatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator+=(const srmatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator+=(const simatrix& B) {
  return fsp_mm_addassign(*this,B);
}

inline imatrix& imatrix::operator-=(const srmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline imatrix& imatrix::operator-=(const simatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator-=(const srmatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator-=(const simatrix& B) {
  return fsp_mm_subassign(*this,B);
}

inline imatrix& imatrix::operator*=(const srmatrix& B) {
  return fsp_mm_multassign<imatrix,srmatrix,sparse_idot,imatrix>(*this,B);
}

inline imatrix& imatrix::operator*=(const simatrix& B) {
  return fsp_mm_multassign<imatrix,simatrix,sparse_idot,imatrix>(*this,B);
}

inline imatrix_slice& imatrix_slice::operator*=(const srmatrix& B) {
  return fsp_mm_multassign<imatrix_slice,srmatrix,sparse_idot,imatrix>(*this,B);
}

inline imatrix_slice& imatrix_slice::operator*=(const simatrix& B) {
  return fsp_mm_multassign<imatrix_slice,simatrix,sparse_idot,imatrix>(*this,B);
}

inline imatrix& imatrix::operator|=(const srmatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline imatrix& imatrix::operator|=(const simatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator|=(const srmatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator|=(const simatrix& B) {
  return fsp_mm_hullassign(*this,B);
}

inline imatrix& imatrix::operator&=(const srmatrix& B) {
  return fsp_mm_intersectassign(*this,B);
}

inline imatrix& imatrix::operator&=(const simatrix& B) {
  return fsp_mm_intersectassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator&=(const srmatrix& B) {
  return fsp_mm_intersectassign(*this,B);
}

inline imatrix_slice& imatrix_slice::operator&=(const simatrix& B) {
  return fsp_mm_intersectassign(*this,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const simatrix& A, const srmatrix& B) {
  return spsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const srmatrix& A, const simatrix& B) {
  return spsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const simatrix& A, const simatrix& B) {
  return spsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const simatrix& A, const rmatrix& B) {
  return spf_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const srmatrix& A, const imatrix& B) {
  return spf_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const simatrix& A, const imatrix& B) {
  return spf_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const imatrix& A, const srmatrix& B) {
  return fsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const rmatrix& A, const simatrix& B) {
  return fsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const imatrix& A, const simatrix& B) {
  return fsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns true iff all elements of A and B are identical.
inline bool operator==(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const simatrix& A, const srmatrix& B) {
  return !spsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const srmatrix& A, const simatrix& B) {
  return !spsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const simatrix& A, const simatrix& B) {
  return !spsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const simatrix& A, const rmatrix& B) {
  return !spf_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const srmatrix& A, const imatrix& B) {
  return !spf_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const simatrix& A, const imatrix& B) {
  return !spf_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const imatrix& A, const srmatrix& B) {
  return !fsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const rmatrix& A, const simatrix& B) {
  return !fsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const imatrix& A, const simatrix& B) {
  return !fsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const imatrix_slice& A, const srmatrix& B) {
  return !fsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const rmatrix_slice& A, const simatrix& B) {
  return !fsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const imatrix_slice& A, const simatrix& B) {
  return !fsp_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const simatrix& A, const rmatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const srmatrix& A, const imatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

//! Element-wise comparison of A and B. Returns false iff all elements of A and B are identical.
inline bool operator!=(const simatrix& A, const imatrix_slice& B) {
  return !spf_mm_comp(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<(const srmatrix& A, const simatrix& B) {
  return spsp_mm_less<srmatrix,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<(const simatrix& A, const simatrix& B) {
  return spsp_mm_less<simatrix,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<(const srmatrix& A, const imatrix& B) {
  return spf_mm_less<srmatrix,imatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<(const simatrix& A, const imatrix& B) {
  return spf_mm_less<simatrix,imatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<(const rmatrix& A, const simatrix& B) {
  return fsp_mm_less<rmatrix,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<(const imatrix& A, const simatrix& B) {
  return fsp_mm_less<imatrix,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_less<rmatrix_slice,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_less<imatrix_slice,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_less<srmatrix,imatrix_slice,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_less<simatrix,imatrix_slice,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<=(const srmatrix& A, const simatrix& B) {
  return spsp_mm_leq<srmatrix,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<=(const simatrix& A, const simatrix& B) {
  return spsp_mm_leq<simatrix,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<=(const srmatrix& A, const imatrix& B) {
  return spf_mm_leq<srmatrix,imatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<=(const simatrix& A, const imatrix& B) {
  return spf_mm_leq<simatrix,imatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<=(const rmatrix& A, const simatrix& B) {
  return fsp_mm_leq<rmatrix,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<=(const imatrix& A, const simatrix& B) {
  return fsp_mm_leq<imatrix,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<=(const rmatrix_slice& A, const simatrix& B) {
  return fsp_mm_leq<rmatrix_slice,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<=(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_leq<imatrix_slice,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<=(const srmatrix& A, const imatrix_slice& B) {
  return spf_mm_leq<srmatrix,imatrix_slice,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator<=(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_leq<simatrix,imatrix_slice,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>(const simatrix& A, const srmatrix& B) {
  return spsp_mm_greater<simatrix,srmatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>(const simatrix& A, const simatrix& B) {
  return spsp_mm_greater<simatrix,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>(const simatrix& A, const rmatrix& B) {
  return spf_mm_greater<simatrix,rmatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>(const simatrix& A, const imatrix& B) {
  return spf_mm_greater<simatrix,imatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>(const imatrix& A, const srmatrix& B) {
  return fsp_mm_greater<imatrix,srmatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>(const imatrix& A, const simatrix& B) {
  return fsp_mm_greater<imatrix,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_greater<imatrix_slice,srmatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_greater<imatrix_slice,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_greater<simatrix,rmatrix_slice,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_greater<simatrix,imatrix_slice,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>=(const simatrix& A, const srmatrix& B) {
  return spsp_mm_geq<simatrix,srmatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>=(const simatrix& A, const simatrix& B) {
  return spsp_mm_geq<simatrix,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>=(const simatrix& A, const rmatrix& B) {
  return spf_mm_geq<simatrix,rmatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>=(const simatrix& A, const imatrix& B) {
  return spf_mm_geq<simatrix,imatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>=(const imatrix& A, const srmatrix& B) {
  return fsp_mm_geq<imatrix,srmatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>=(const imatrix& A, const simatrix& B) {
  return fsp_mm_geq<imatrix,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>=(const imatrix_slice& A, const srmatrix& B) {
  return fsp_mm_geq<imatrix_slice,srmatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>=(const imatrix_slice& A, const simatrix& B) {
  return fsp_mm_geq<imatrix_slice,simatrix,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>=(const simatrix& A, const rmatrix_slice& B) {
  return spf_mm_geq<simatrix,rmatrix_slice,interval>(A,B);
}

//! Element-wise comparison of A and B.
inline bool operator>=(const simatrix& A, const imatrix_slice& B) {
  return spf_mm_geq<simatrix,imatrix_slice,interval>(A,B);
}

//! Element-wise logical negation of A. Return true if all elements of A are equal to zero
inline bool operator!(const simatrix& A) {
  return sp_m_not(A);
}

//! Standard output operator for sparse matrices
/**
 *  The output format is set by global flags, default is dense output.
 *  Use cout << SparseInOut; for sparse output or cout << MatrixMarketInOut; for
 *  output in matrix market format.
 */
inline std::ostream& operator<<(std::ostream& os, const simatrix& A) {
  return sp_m_output<simatrix,interval>(os,A);
}

//! Standard input operator for sparse matrices
/**
 *  The input format is set by global flags, default is dense input.
 *  Use cout << SparseInOut; for sparse input or cout << MatrixMarketInOut; for
 *  input in matrix market format.
 */
inline std::istream& operator>>(std::istream& is, simatrix& A) {
  return sp_m_input<simatrix,interval>(is,A);
}

//! A slice of a sparse real interval matrix
/*!
    Represents a slice of a sparse real interval matrix. This helper class provides read and write access to such a slice using the standard operators. It should normally not be necessary
    for the user to explicitly work with this data type, which is why the constructors are private.
 */
class simatrix_slice {
  public:
    simatrix  A;
    simatrix* M; //Originalmatrix

  private:
    simatrix_slice(simatrix& Mat, int sl1l, int sl1u, int sl2l, int sl2u) {    
        A.lb1 = sl1l;
        A.lb2 = sl2l;
        A.ub1 = sl1u;
        A.ub2 = sl2u;
        A.m   = sl1u-sl1l+1;
        A.n   = sl2u-sl2l+1;
 
        //Kopieren der Werte aus A
        A.p = std::vector<int>(A.n+1, 0);
        A.ind.reserve(A.m + A.n);
        A.x.reserve(A.m + A.n);

        for(int i=0 ; i<A.n ; i++) {
           A.p[i+1] = A.p[i];
           for(int j=Mat.p[sl2l-Mat.lb2+i] ; j<Mat.p[sl2l-Mat.lb2+i+1] ; j++) {
              if(Mat.ind[j] >= sl1l-Mat.lb1  &&  Mat.ind[j] <= sl1u-Mat.lb1) {
                A.ind.push_back(Mat.ind[j]-(sl1l-Mat.lb1));
                A.x.push_back(Mat.x[j]);
                A.p[i+1]++;
              }
           }
        }

        //Zeiger auf A fuer Datenmanipulationen
        M = &Mat;
    }

    simatrix_slice(const simatrix& Mat, int sl1l, int sl1u, int sl2l, int sl2u) {    
        A.lb1 = sl1l;
        A.lb2 = sl2l;
        A.ub1 = sl1u;
        A.ub2 = sl2u;
        A.m   = sl1u-sl1l+1;
        A.n   = sl2u-sl2l+1;
 
        //Kopieren der Werte aus A
        A.p = std::vector<int>(A.n+1, 0);
        A.ind.reserve(A.m + A.n);
        A.x.reserve(A.m + A.n);

        for(int i=0 ; i<A.n ; i++) {
           A.p[i+1] = A.p[i];
           for(int j=Mat.p[sl2l-Mat.lb2+i] ; j<Mat.p[sl2l-Mat.lb2+i+1] ; j++) {
              if(Mat.ind[j] >= sl1l-Mat.lb1  &&  Mat.ind[j] <= sl1u-Mat.lb1) {
                A.ind.push_back(Mat.ind[j]-(sl1l-Mat.lb1));
                A.x.push_back(Mat.x[j]);
                A.p[i+1]++;
              }
           }
        }

        //Zeiger auf A fuer Datenmanipulationen
        M = const_cast<simatrix*>(&Mat); //Vorgehen noetig um Schreibweise A[i][j] bei auslesen von const A zu ermoeglichen
    }


  public:
    //! Assing C to all elements of the slice
    simatrix_slice& operator=(const real& C) {
      return sl_ms_assign<simatrix_slice, real, std::vector<interval>::iterator, interval>(*this,C);
    }

    //! Assing C to all elements of the slice
    simatrix_slice& operator=(const interval& C) {
      return sl_ms_assign<simatrix_slice, interval, std::vector<interval>::iterator, interval>(*this,C);
    }

    //! Assing C to the slice
    simatrix_slice& operator=(const srmatrix& C) {
      return slsp_mm_assign<simatrix_slice, srmatrix, std::vector<interval>::iterator>(*this,C);
    }

    //! Assing C to the slice
    simatrix_slice& operator=(const simatrix& C) {
      return slsp_mm_assign<simatrix_slice, simatrix, std::vector<interval>::iterator>(*this,C);
    }

    //! Assing C to the slice
    simatrix_slice& operator=(const rmatrix& C) {
      return slf_mm_assign<simatrix_slice, rmatrix, std::vector<interval>::iterator, interval>(*this,C);
    }

    //! Assing C to the slice
    simatrix_slice& operator=(const imatrix& C) {
      return slf_mm_assign<simatrix_slice, imatrix, std::vector<interval>::iterator, interval>(*this,C);
    }

    //! Assing C to the slice
    simatrix_slice& operator=(const rmatrix_slice& C) {
      return slf_mm_assign<simatrix_slice, rmatrix_slice, std::vector<interval>::iterator, interval>(*this,C);
    }

    //! Assing C to the slice
    simatrix_slice& operator=(const imatrix_slice& C) {
      return slf_mm_assign<simatrix_slice, imatrix_slice, std::vector<interval>::iterator, interval>(*this,C);
    }

    //! Assing C to the slice
    simatrix_slice& operator=(const srmatrix_slice& C) {
      *this = C.A;
      return *this;
    }

    //! Assing C to the slice
    simatrix_slice& operator=(const simatrix_slice& C) {
      *this = C.A;
      return *this;
    }

    //! Assigns the product of the sparse slice and M to the slice
    simatrix_slice& operator*=(const srmatrix_slice& M) {
      *this = A*M.A;
      return *this;
    }

    //! Assigns the product of the sparse slice and M to the slice
    simatrix_slice& operator*=(const simatrix_slice& M) {
      *this = A*M.A;
      return *this;
    }

    //! Assigns the product of the sparse slice and M to the slice
    simatrix_slice& operator*=(const srmatrix& M) {
      *this = A*M;
      return *this;
    }

    //! Assigns the product of the sparse slice and M to the slice
    simatrix_slice& operator*=(const simatrix& M) {
      *this = A*M;
      return *this;
    }

    //! Assigns the product of the sparse slice and M to the slice
    simatrix_slice& operator*=(const rmatrix& M) {
      *this = A*M;
      return *this;
    }

    //! Assigns the product of the sparse slice and M to the slice
    simatrix_slice& operator*=(const imatrix& M) {
      *this = A*M;
      return *this;
    }

    //! Assigns the product of the sparse slice and M to the slice
    simatrix_slice& operator*=(const rmatrix_slice& M) {
      *this = A*M;
      return *this;
    }

    //! Assigns the product of the sparse slice and M to the slice
    simatrix_slice& operator*=(const imatrix_slice& M) {
      *this = A*M;
      return *this;
    }

    //! Assigns the component wise product of the sparse slice and r to the slice
    simatrix_slice& operator*=(const real& r) {
      *this = A*r;
      return *this;
    }

    //! Assigns the component wise product of the sparse slice and r to the slice
    simatrix_slice& operator*=(const interval& r) {
      *this = A*r;
      return *this;
    }

    //! Assigns the component wise division of the sparse slice and M to the slice
    simatrix_slice& operator/=(const real& r) {
      *this = A/r;
      return *this;
    }

    //! Assigns the component wise division of the sparse slice and M to the slice
    simatrix_slice& operator/=(const interval& r) {
      *this = A/r;
      return *this;
    }

    //! Assigns the element wise sum of the sparse slice and M to the slice
    simatrix_slice& operator+=(const srmatrix_slice& M) {
      *this = A+M.A;
      return *this;
    } 

    //! Assigns the element wise sum of the sparse slice and M to the slice
    simatrix_slice& operator+=(const simatrix_slice& M) {
      *this = A+M.A;
      return *this;
    } 

    //! Assigns the element wise sum of the sparse slice and M to the slice
    simatrix_slice& operator+=(const srmatrix& M) {
      *this = A+M;
      return *this;
    } 

    //! Assigns the element wise sum of the sparse slice and M to the slice
    simatrix_slice& operator+=(const simatrix& M) {
      *this = A+M;
      return *this;
    } 

    //! Assigns the element wise sum of the sparse slice and M to the slice
    simatrix_slice& operator+=(const rmatrix& M) {
      *this = A+M;
      return *this;
    } 

    //! Assigns the element wise sum of the sparse slice and M to the slice
    simatrix_slice& operator+=(const imatrix& M) {
      *this = A+M;
      return *this;
    } 

    //! Assigns the element wise sum of the sparse slice and M to the slice
    simatrix_slice& operator+=(const rmatrix_slice& M) {
      *this = A+M;
      return *this;
    } 

    //! Assigns the element wise sum of the sparse slice and M to the slice
    simatrix_slice& operator+=(const imatrix_slice& M) {
      *this = A+M;
      return *this;
    } 

    //! Assigns the element wise difference of the sparse slice and M to the slice
    simatrix_slice& operator-=(const srmatrix_slice& M) {
      *this = A-M.A;
      return *this;
    } 

    //! Assigns the element wise difference of the sparse slice and M to the slice
    simatrix_slice& operator-=(const simatrix_slice& M) {
      *this = A-M.A;
      return *this;
    } 

    //! Assigns the element wise difference of the sparse slice and M to the slice
    simatrix_slice& operator-=(const srmatrix& M) {
      *this = A-M;
      return *this;
    } 

    //! Assigns the element wise difference of the sparse slice and M to the slice
    simatrix_slice& operator-=(const simatrix& M) {
      *this = A-M;
      return *this;
    } 

    //! Assigns the element wise difference of the sparse slice and M to the slice
    simatrix_slice& operator-=(const rmatrix& M) {
      *this = A-M;
      return *this;
    } 

    //! Assigns the element wise difference of the sparse slice and M to the slice
    simatrix_slice& operator-=(const imatrix& M) {
      *this = A-M;
      return *this;
    } 

    //! Assigns the element wise difference of the sparse slice and M to the slice
    simatrix_slice& operator-=(const rmatrix_slice& M) {
      *this = A-M;
      return *this;
    }

    //! Assigns the element wise difference of the sparse slice and M to the slice
    simatrix_slice& operator-=(const imatrix_slice& M) {
      *this = A-M;
      return *this;
    }

    //! Assigns the element wise convex hull of the sparse slice and M to the slice
    simatrix_slice& operator|=(const srmatrix_slice& M) {
      *this = A|M.A;
      return *this;
    } 

    //! Assigns the element wise convex hull of the sparse slice and M to the slice
    simatrix_slice& operator|=(const simatrix_slice& M) {
      *this = A|M.A;
      return *this;
    } 

    //! Assigns the element wise convex hull of the sparse slice and M to the slice
    simatrix_slice& operator|=(const srmatrix& M) {
      *this = A|M;
      return *this;
    } 

    //! Assigns the element wise convex hull of the sparse slice and M to the slice
    simatrix_slice& operator|=(const simatrix& M) {
      *this = A|M;
      return *this;
    } 

    //! Assigns the element wise convex hull of the sparse slice and M to the slice
    simatrix_slice& operator|=(const rmatrix& M) {
      *this = A|M;
      return *this;
    } 

    //! Assigns the element wise convex hull of the sparse slice and M to the slice
    simatrix_slice& operator|=(const imatrix& M) {
      *this = A|M;
      return *this;
    } 

    //! Assigns the element wise convex hull of the sparse slice and M to the slice
    simatrix_slice& operator|=(const rmatrix_slice& M) {
      *this = A|M;
      return *this;
    } 

    //! Assigns the element wise convex hull of the sparse slice and M to the slice
    simatrix_slice& operator|=(const imatrix_slice& M) {
      *this = A|M;
      return *this;
    } 

    //! Returns a copy of the element (i,j) of the matrix
    /*!
        This operators can only be usd for read access. Note that accessing single elements of a sparse matrix is more expensive than for dense matrices and 
        should in general be avoided.
     */
    const interval operator()(const int i, const int j) const {
#if(CXSC_INDEX_CHECK)
      if(i<A.lb1 || i>A.ub1 || j<A.lb2 || j>A.ub2)
        cxscthrow(ELEMENT_NOT_IN_VEC("simatrix_slice::operator()(int, int)"));
#endif
      interval r = A(i,j);
      return r;
    }

    //! Returns a reference to the element (i,j) of the matrix
    /*!
        Returns a reference to the (i,j)-th element. If the element is not explicitly stored, it is added as an explicit zero entry to the data structure.
        Using this function is faster than using A[i][j], since no temporary subvecto object must be created.
     */
    interval& element(const int i, const int j) {
#if(CXSC_INDEX_CHECK)
      if(i<A.lb1 || i>A.ub1 || j<A.lb2 || j>A.ub2)
        cxscthrow(ELEMENT_NOT_IN_VEC("simatrix_slice::element(int, int)"));
#endif
      return M->element(i,j);
    }

    //! Returns a row of the matrix
    simatrix_subv operator[](const int);
    //! Returns a column of the matrix
    simatrix_subv operator[](const cxscmatrix_column&);
    //! Returns a row of the matrix
    const simatrix_subv operator[](const int) const;
    //! Returns a column of the matrix
    const simatrix_subv operator[](const cxscmatrix_column&) const;

    friend std::ostream& operator<<(std::ostream&, const simatrix_slice&);

    friend int Lb(const simatrix_slice&, const int);
    friend int Ub(const simatrix_slice&, const int);
    friend srmatrix Inf(const simatrix_slice&);
    friend srmatrix Sup(const simatrix_slice&);
    friend simatrix abs(const simatrix_slice&);
    friend srmatrix mid(const simatrix_slice&);
    friend srmatrix diam(const simatrix_slice&);
    friend int RowLen(const simatrix_slice&);
    friend int ColLen(const simatrix_slice&);

    friend class srmatrix;
    friend class srmatrix_subv;
    friend class srvector;
    friend class simatrix;
    friend class simatrix_subv;
    friend class sivector;
    friend class scimatrix;
    friend class scimatrix_subv;
    friend class scimatrix_slice;
    friend class scivector;
    friend class imatrix;
    friend class cimatrix;


#include "matrix_friend_declarations.inl"    
};

inline imatrix::imatrix(const srmatrix_slice& A) {
  dat = new interval[A.A.m*A.A.n];
  lb1 = A.A.lb1; lb2 = A.A.lb2; ub1 = A.A.ub1; ub2 = A.A.ub2;
  xsize = A.A.n;
  ysize = A.A.m;
  *this = 0.0;
  for(int j=0 ; j<A.A.n ; j++) {
     for(int k=A.A.p[j] ; k<A.A.p[j+1] ; k++) {
        dat[A.A.ind[k]*A.A.n+j] = A.A.x[k];
     }
  }
}

inline imatrix::imatrix(const simatrix_slice& A) {
  dat = new interval[A.A.m*A.A.n];
  lb1 = A.A.lb1; lb2 = A.A.lb2; ub1 = A.A.ub1; ub2 = A.A.ub2;
  xsize = A.A.n;
  ysize = A.A.m;
  *this = 0.0;
  for(int j=0 ; j<A.A.n ; j++) {
     for(int k=A.A.p[j] ; k<A.A.p[j+1] ; k++) {
        dat[A.A.ind[k]*A.A.n+j] = A.A.x[k];
     }
  }
}

//! Returns the number columns of the matrix slice
inline int RowLen(const simatrix_slice& S) {
  return RowLen(S.A);
}

//! Returns the number of rows of the matrix slice
inline int ColLen(const simatrix_slice& S) {
  return ColLen(S.A);
}

inline simatrix_slice simatrix::operator()(const int i, const int j, const int k, const int l) {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || j>ub1 || k<lb2 || l>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator()(int, int, int, int)"));
#endif
  return simatrix_slice(*this, i, j, k, l);
}

inline const simatrix_slice simatrix::operator()(const int i, const int j, const int k, const int l) const {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || j>ub1 || k<lb2 || l>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator()(int, int, int, int) const"));
#endif
  return simatrix_slice(*this, i, j, k, l);
}

//! Returns the lower index bound of the rows (if i==ROW) or columns (if i==COL) of the slice
inline int Lb(const simatrix_slice& S, const int i) {
  return Lb(S.A, i);
}

//! Returns the upper index bound of the rows (if i==ROW) or columns (if i==COL) of the slice
inline int Ub(const simatrix_slice& S, const int i) {
  return Ub(S.A, i);
}

//! Returns the infimum of the slice S
inline srmatrix Inf(const simatrix_slice& S) {
  return Inf(S.A);
}

//! Returns the supremum of the slice S
inline srmatrix Sup(const simatrix_slice& S) {
  return Sup(S.A);
}

//! Returns the elementwise absolute value of S
inline simatrix abs(const simatrix_slice& S) {
  return abs(S.A);
}

//! Returns the elementwise midpoint of S
inline srmatrix mid(const simatrix_slice& S) {
  return mid(S.A);
}

//! Returns the elementwise diameter of S
inline srmatrix diam(const simatrix_slice& S) {
  return diam(S.A);
}

inline simatrix::simatrix(const srmatrix_slice& S) {
  m = S.A.m;
  n = S.A.n;
  lb1 = S.A.lb1;
  ub1 = S.A.ub1;
  lb2 = S.A.lb2;
  ub2 = S.A.ub2;
  *this = S.A;
}

inline simatrix::simatrix(const simatrix_slice& S) {
  m = S.A.m;
  n = S.A.n;
  lb1 = S.A.lb1;
  ub1 = S.A.ub1;
  lb2 = S.A.lb2;
  ub2 = S.A.ub2;
  *this = S.A;
}

inline simatrix& simatrix::operator=(const srmatrix_slice& S) {
  *this = S.A;
  return *this;
}

inline simatrix& simatrix::operator=(const simatrix_slice& S) {
  *this = S.A;
  return *this;
}

inline imatrix& imatrix::operator=(const srmatrix_slice& M) {
  *this = rmatrix(M);
  return *this;
}

inline imatrix& imatrix::operator=(const simatrix_slice& M) {
  *this = imatrix(M);
  return *this;
}

inline imatrix& imatrix::operator+=(const srmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline imatrix& imatrix::operator+=(const simatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator+=(const srmatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator+=(const simatrix_slice& M) {
  *this += M.A;
  return *this;
}

inline imatrix& imatrix::operator-=(const srmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline imatrix& imatrix::operator-=(const simatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator-=(const srmatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator-=(const simatrix_slice& M) {
  *this -= M.A;
  return *this;
}

inline imatrix& imatrix::operator*=(const srmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline imatrix& imatrix::operator*=(const simatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator*=(const srmatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator*=(const simatrix_slice& M) {
  *this *= M.A;
  return *this;
}

inline imatrix& imatrix::operator|=(const srmatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline imatrix& imatrix::operator|=(const simatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator|=(const srmatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator|=(const simatrix_slice& M) {
  *this |= M.A;
  return *this;
}

inline imatrix& imatrix::operator&=(const srmatrix_slice& M) {
  *this &= M.A;
  return *this;
}

inline imatrix& imatrix::operator&=(const simatrix_slice& M) {
  *this &= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator&=(const srmatrix_slice& M) {
  *this &= M.A;
  return *this;
}

inline imatrix_slice& imatrix_slice::operator&=(const simatrix_slice& M) {
  *this &= M.A;
  return *this;
}

//! Unary negation operator for matrix slices
inline simatrix operator-(const simatrix_slice& M) {
  return sp_m_negative<simatrix,simatrix>(M.A);
}

//! Unary operator+ for matrix slices
inline simatrix operator+(const simatrix_slice& M) {
  return M.A;
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline simatrix operator*(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_mult<simatrix,srmatrix,simatrix,sparse_idot,interval>(M1.A,M2.A);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline simatrix operator*(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_mult<srmatrix,simatrix,simatrix,sparse_idot,interval>(M1.A,M2.A);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline simatrix operator*(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_mult<simatrix,simatrix,simatrix,sparse_idot,interval>(M1.A,M2.A);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline simatrix operator*(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_mult<simatrix,srmatrix,simatrix,sparse_idot,interval>(M1.A,M2);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline simatrix operator*(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_mult<srmatrix,simatrix,simatrix,sparse_idot,interval>(M1.A,M2);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline simatrix operator*(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_mult<simatrix,simatrix,simatrix,sparse_idot,interval>(M1.A,M2);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline simatrix operator*(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_mult<simatrix,srmatrix,simatrix,sparse_idot,interval>(M1,M2.A);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline simatrix operator*(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_mult<srmatrix,simatrix,simatrix,sparse_idot,interval>(M1,M2.A);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline simatrix operator*(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_mult<simatrix,simatrix,simatrix,sparse_idot,interval>(M1,M2.A);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_mult<simatrix,rmatrix,imatrix,sparse_idot>(M1.A,M2);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_mult<srmatrix,imatrix,imatrix,sparse_idot>(M1.A,M2);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_mult<simatrix,imatrix,imatrix,sparse_idot>(M1.A,M2);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_mult<imatrix,srmatrix,imatrix,sparse_idot>(M1,M2.A);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_mult<rmatrix,simatrix,imatrix,sparse_idot>(M1,M2.A);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_mult<imatrix,simatrix,imatrix,sparse_idot>(M1,M2.A);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_mult<simatrix,rmatrix_slice,imatrix,sparse_idot>(M1.A,M2);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_mult<srmatrix,imatrix_slice,imatrix,sparse_idot>(M1.A,M2);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_mult<simatrix,imatrix_slice,imatrix,sparse_idot>(M1.A,M2);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_mult<imatrix,srmatrix,imatrix,sparse_idot>(M1,M2.A);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_mult<rmatrix,simatrix,imatrix,sparse_idot>(M1,M2.A);
}

//! Returns the product of the matrices M1 and M2.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline imatrix operator*(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_mult<imatrix,simatrix,imatrix,sparse_idot>(M1,M2.A);
}

//! Returns the product of the matrix M and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline sivector operator*(const simatrix_slice& M, const srvector& v) {
  return spsp_mv_mult<simatrix,srvector,sivector,sparse_idot,interval>(M.A,v);
}

//! Returns the product of the matrix M and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline sivector operator*(const srmatrix_slice& M, const sivector& v) {
  return spsp_mv_mult<srmatrix,sivector,sivector,sparse_idot,interval>(M.A,v);
}

//! Returns the product of the matrix M and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline sivector operator*(const simatrix_slice& M, const sivector& v) {
  return spsp_mv_mult<simatrix,sivector,sivector,sparse_idot,interval>(M.A,v);
}

//! Returns the product of the matrix M and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline sivector operator*(const simatrix_slice& M, const srvector_slice& v) {
  return spsl_mv_mult<simatrix,srvector_slice,sivector,sparse_idot,interval>(M.A,v);
}

//! Returns the product of the matrix M and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline sivector operator*(const srmatrix_slice& M, const sivector_slice& v) {
  return spsl_mv_mult<srmatrix,sivector_slice,sivector,sparse_idot,interval>(M.A,v);
}

//! Returns the product of the matrix M and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline sivector operator*(const simatrix_slice& M, const sivector_slice& v) {
  return spsl_mv_mult<simatrix,sivector_slice,sivector,sparse_idot,interval>(M.A,v);
}

//! Returns the product of the matrix M and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const simatrix_slice& M, const rvector& v) {
  return spf_mv_mult<simatrix,rvector,ivector,sparse_idot>(M.A,v);
}

//! Returns the product of the matrix M and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const srmatrix_slice& M, const ivector& v) {
  return spf_mv_mult<srmatrix,ivector,ivector,sparse_idot>(M.A,v);
}

//! Returns the product of the matrix M and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const simatrix_slice& M, const ivector& v) {
  return spf_mv_mult<simatrix,ivector,ivector,sparse_idot>(M.A,v);
}

//! Returns the product of the matrix M and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const simatrix_slice& M, const rvector_slice& v) {
  return spf_mv_mult<simatrix,rvector_slice,ivector,sparse_idot>(M.A,v);
}

//! Returns the product of the matrix M and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const srmatrix_slice& M, const ivector_slice& v) {
  return spf_mv_mult<srmatrix,ivector_slice,ivector,sparse_idot>(M.A,v);
}

//! Returns the product of the matrix M and the vector v.
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline ivector operator*(const simatrix_slice& M, const ivector_slice& v) {
  return spf_mv_mult<simatrix,ivector_slice,ivector,sparse_idot>(M.A,v);
}

//! Returns the element wise division of the matrix M and r.
inline simatrix operator/(const simatrix_slice& M, const real& r) {
  return sp_ms_div<simatrix,real,simatrix>(M.A,r);
}

//! Returns the element wise division of the matrix M and r.
inline simatrix operator/(const simatrix_slice& M, const interval& r) {
  return sp_ms_div<simatrix,interval,simatrix>(M.A,r);
}

//! Returns the element wise division of the matrix M and r.
inline simatrix operator/(const srmatrix_slice& M, const interval& r) {
  return sp_ms_div<srmatrix,interval,simatrix>(M.A,r);
}

//! Returns the element wise product of the matrix M and r.
inline simatrix operator*(const simatrix_slice& M, const real& r) {
  return sp_ms_mult<simatrix,real,simatrix>(M.A,r);
}

//! Returns the element wise product of the matrix M and r.
inline simatrix operator*(const simatrix_slice& M, const interval& r) {
  return sp_ms_mult<simatrix,interval,simatrix>(M.A,r);
}

//! Returns the element wise product of the matrix M and r.
inline simatrix operator*(const srmatrix_slice& M, const interval& r) {
  return sp_ms_mult<srmatrix,interval,simatrix>(M.A,r);
}

//! Returns the element wise product of the matrix M and r.
inline simatrix operator*(const real& r, const simatrix_slice& M) {
  return sp_sm_mult<real,simatrix,simatrix>(r,M.A);
}

//! Returns the element wise product of the matrix M and r.
inline simatrix operator*(const interval& r, const srmatrix_slice& M) {
  return sp_sm_mult<interval,srmatrix,simatrix>(r,M.A);
}

//! Returns the element wise product of the matrix M and r.
inline simatrix operator*(const interval& r, const simatrix_slice& M) {
  return sp_sm_mult<interval,simatrix,simatrix>(r,M.A);
}

//! Returns the element-wise sum of M1 and M2
inline simatrix operator+(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_add<simatrix,srmatrix,simatrix,interval>(M1.A,M2.A);
}

//! Returns the element-wise sum of M1 and M2
inline simatrix operator+(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_add<srmatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

//! Returns the element-wise sum of M1 and M2
inline simatrix operator+(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_add<simatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

//! Returns the element-wise sum of M1 and M2
inline simatrix operator+(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_add<simatrix,srmatrix,simatrix,interval>(M1.A,M2);
}

//! Returns the element-wise sum of M1 and M2
inline simatrix operator+(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_add<srmatrix,simatrix,simatrix,interval>(M1.A,M2);
}

//! Returns the element-wise sum of M1 and M2
inline simatrix operator+(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_add<simatrix,simatrix,simatrix,interval>(M1.A,M2);
}

//! Returns the element-wise sum of M1 and M2
inline simatrix operator+(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_add<simatrix,srmatrix,simatrix,interval>(M1,M2.A);
}

//! Returns the element-wise sum of M1 and M2
inline simatrix operator+(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_add<srmatrix,simatrix,simatrix,interval>(M1,M2.A);
}

//! Returns the element-wise sum of M1 and M2
inline simatrix operator+(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_add<simatrix,simatrix,simatrix,interval>(M1,M2.A);
}

//! Returns the element-wise sum of M1 and M2
inline imatrix operator+(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_add<simatrix,rmatrix,imatrix>(M1.A,M2);
}

//! Returns the element-wise sum of M1 and M2
inline imatrix operator+(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_add<srmatrix,imatrix,imatrix>(M1.A,M2);
}

//! Returns the element-wise sum of M1 and M2
inline imatrix operator+(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_add<simatrix,imatrix,imatrix>(M1.A,M2);
}

//! Returns the element-wise sum of M1 and M2
inline imatrix operator+(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_add<imatrix,srmatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise sum of M1 and M2
inline imatrix operator+(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_add<rmatrix,simatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise sum of M1 and M2
inline imatrix operator+(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_add<imatrix,simatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise sum of M1 and M2
inline imatrix operator+(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_add<simatrix,rmatrix_slice,imatrix>(M1.A,M2);
}

//! Returns the element-wise sum of M1 and M2
inline imatrix operator+(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_add<srmatrix,imatrix_slice,imatrix>(M1.A,M2);
}

//! Returns the element-wise sum of M1 and M2
inline imatrix operator+(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_add<simatrix,imatrix_slice,imatrix>(M1.A,M2);
}

//! Returns the element-wise sum of M1 and M2
inline imatrix operator+(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_add<imatrix_slice,srmatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise sum of M1 and M2
inline imatrix operator+(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_add<rmatrix_slice,simatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise sum of M1 and M2
inline imatrix operator+(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_add<imatrix_slice,simatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise difference of M1 and M2
inline simatrix operator-(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_sub<simatrix,srmatrix,simatrix,interval>(M1.A,M2.A);
}

//! Returns the element-wise difference of M1 and M2
inline simatrix operator-(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_sub<srmatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

//! Returns the element-wise difference of M1 and M2
inline simatrix operator-(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_sub<simatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

//! Returns the element-wise difference of M1 and M2
inline simatrix operator-(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_sub<simatrix,srmatrix,simatrix,interval>(M1.A,M2);
}

//! Returns the element-wise difference of M1 and M2
inline simatrix operator-(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_sub<srmatrix,simatrix,simatrix,interval>(M1.A,M2);
}

//! Returns the element-wise difference of M1 and M2
inline simatrix operator-(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_sub<simatrix,simatrix,simatrix,interval>(M1.A,M2);
}

//! Returns the element-wise difference of M1 and M2
inline simatrix operator-(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_sub<simatrix,srmatrix,simatrix,interval>(M1,M2.A);
}

//! Returns the element-wise difference of M1 and M2
inline simatrix operator-(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_sub<srmatrix,simatrix,simatrix,interval>(M1,M2.A);
}

//! Returns the element-wise difference of M1 and M2
inline simatrix operator-(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_sub<simatrix,simatrix,simatrix,interval>(M1,M2.A);
}

//! Returns the element-wise difference of M1 and M2
inline imatrix operator-(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_sub<simatrix,rmatrix,imatrix>(M1.A,M2);
}

//! Returns the element-wise difference of M1 and M2
inline imatrix operator-(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_sub<srmatrix,imatrix,imatrix>(M1.A,M2);
}

//! Returns the element-wise difference of M1 and M2
inline imatrix operator-(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_sub<simatrix,imatrix,imatrix>(M1.A,M2);
}

//! Returns the element-wise difference of M1 and M2
inline imatrix operator-(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_sub<imatrix,srmatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise difference of M1 and M2
inline imatrix operator-(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_sub<rmatrix,simatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise difference of M1 and M2
inline imatrix operator-(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_sub<imatrix,simatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise difference of M1 and M2
inline imatrix operator-(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_sub<simatrix,rmatrix_slice,imatrix>(M1.A,M2);
}

//! Returns the element-wise difference of M1 and M2
inline imatrix operator-(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_sub<srmatrix,imatrix_slice,imatrix>(M1.A,M2);
}

//! Returns the element-wise difference of M1 and M2
inline imatrix operator-(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_sub<simatrix,imatrix_slice,imatrix>(M1.A,M2);
}

//! Returns the element-wise difference of M1 and M2
inline imatrix operator-(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_sub<imatrix_slice,srmatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise difference of M1 and M2
inline imatrix operator-(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_sub<rmatrix_slice,simatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise difference of M1 and M2
inline imatrix operator-(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_sub<imatrix_slice,simatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline simatrix operator|(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_hull<simatrix,srmatrix,simatrix,interval>(M1.A,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline simatrix operator|(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_hull<srmatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline simatrix operator|(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_hull<simatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline simatrix operator|(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_hull<simatrix,srmatrix,simatrix,interval>(M1.A,M2);
}

//! Returns the element-wise convex hull of M1 and M2
inline simatrix operator|(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_hull<srmatrix,simatrix,simatrix,interval>(M1.A,M2);
}

//! Returns the element-wise convex hull of M1 and M2
inline simatrix operator|(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_hull<simatrix,simatrix,simatrix,interval>(M1.A,M2);
}

//! Returns the element-wise convex hull of M1 and M2
inline simatrix operator|(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_hull<simatrix,srmatrix,simatrix,interval>(M1,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline simatrix operator|(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_hull<srmatrix,simatrix,simatrix,interval>(M1,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline simatrix operator|(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_hull<simatrix,simatrix,simatrix,interval>(M1,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_hull<simatrix,rmatrix,imatrix>(M1.A,M2);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_hull<srmatrix,imatrix,imatrix>(M1.A,M2);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_hull<simatrix,imatrix,imatrix>(M1.A,M2);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_hull<imatrix,srmatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_hull<rmatrix,simatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_hull<imatrix,simatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_hull<simatrix,rmatrix_slice,imatrix>(M1.A,M2);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_hull<srmatrix,imatrix_slice,imatrix>(M1.A,M2);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_hull<simatrix,imatrix_slice,imatrix>(M1.A,M2);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_hull<imatrix_slice,srmatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_hull<rmatrix_slice,simatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_hull<imatrix_slice,simatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline simatrix operator|(const srmatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_hull<srmatrix,srmatrix,simatrix,interval>(M1.A,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline simatrix operator|(const srmatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_hull<srmatrix,srmatrix,simatrix,interval>(M1.A,M2);
}

//! Returns the element-wise convex hull of M1 and M2
inline simatrix operator|(const srmatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_hull<srmatrix,srmatrix,simatrix,interval>(M1,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const srmatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_hull<srmatrix,rmatrix,imatrix>(M1.A,M2);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const rmatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_hull<rmatrix,srmatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const srmatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_hull<srmatrix,rmatrix_slice,imatrix>(M1.A,M2);
}

//! Returns the element-wise convex hull of M1 and M2
inline imatrix operator|(const rmatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_hull<rmatrix_slice,srmatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise intersection of M1 and M2
inline simatrix operator&(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_intersect<simatrix,simatrix,simatrix,interval>(M1.A,M2.A);
}

//! Returns the element-wise intersection of M1 and M2
inline simatrix operator&(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_intersect<simatrix,simatrix,simatrix,interval>(M1.A,M2);
}

//! Returns the element-wise intersection of M1 and M2
inline simatrix operator&(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_intersect<simatrix,simatrix,simatrix,interval>(M1,M2.A);
}

//! Returns the element-wise intersection of M1 and M2
inline imatrix operator&(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_intersect<simatrix,imatrix,imatrix>(M1.A,M2);
}

//! Returns the element-wise intersection of M1 and M2
inline imatrix operator&(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_intersect<imatrix,simatrix,imatrix>(M1,M2.A);
}

//! Returns the element-wise intersection of M1 and M2
inline imatrix operator&(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_intersect<simatrix,imatrix_slice,imatrix>(M1.A,M2);
}

//! Returns the element-wise intersection of M1 and M2
inline imatrix operator&(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_intersect<imatrix_slice,simatrix,imatrix>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_comp(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator==(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const simatrix_slice& M1, const simatrix_slice& M2) {
  return !spsp_mm_comp(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const simatrix_slice& M1, const srmatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const srmatrix_slice& M1, const simatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const simatrix_slice& M1, const simatrix& M2) {
  return !spsp_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const simatrix& M1, const srmatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const srmatrix& M1, const simatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const simatrix& M1, const simatrix_slice& M2) {
  return !spsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const simatrix_slice& M1, const rmatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const srmatrix_slice& M1, const imatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const simatrix_slice& M1, const imatrix& M2) {
  return !spf_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const imatrix& M1, const srmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const rmatrix& M1, const simatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const imatrix& M1, const simatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const imatrix_slice& M1, const simatrix_slice& M2) {
  return !fsp_mm_comp(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator!=(const simatrix_slice& M1, const imatrix_slice& M2) {
  return !spf_mm_comp(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_less<srmatrix,simatrix,interval>(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_less<simatrix,simatrix,interval>(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_less<srmatrix,simatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_less<simatrix,simatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_less<srmatrix,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_less<simatrix,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_less<srmatrix,imatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_less<simatrix,imatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_less<rmatrix,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_less<imatrix,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_less<rmatrix_slice,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_less<imatrix_slice,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_less<srmatrix,imatrix_slice,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator<(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_less<simatrix,imatrix_slice,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const srmatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_leq<srmatrix,simatrix,interval>(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_leq<simatrix,simatrix,interval>(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const srmatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_leq<srmatrix,simatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_leq<simatrix,simatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const srmatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_leq<srmatrix,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_leq<simatrix,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const srmatrix_slice& M1, const imatrix& M2) {
  return spf_mm_leq<srmatrix,imatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_leq<simatrix,imatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const rmatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_leq<rmatrix,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_leq<imatrix,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const rmatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_leq<rmatrix_slice,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_leq<imatrix_slice,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const srmatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_leq<srmatrix,imatrix_slice,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator<=(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_leq<simatrix,imatrix_slice,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_greater<simatrix,srmatrix,interval>(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_greater<simatrix,simatrix,interval>(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_greater<simatrix,srmatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_greater<simatrix,simatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_greater<simatrix,srmatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_greater<simatrix,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_greater<simatrix,rmatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_greater<simatrix,imatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_greater<imatrix,srmatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_greater<imatrix,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_greater<imatrix,srmatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_greater<imatrix_slice,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_greater<simatrix,rmatrix_slice,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator>(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_greater<simatrix,imatrix_slice,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const simatrix_slice& M1, const srmatrix_slice& M2) {
  return spsp_mm_geq<simatrix,srmatrix,interval>(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const simatrix_slice& M1, const simatrix_slice& M2) {
  return spsp_mm_geq<simatrix,simatrix,interval>(M1.A,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const simatrix_slice& M1, const srmatrix& M2) {
  return spsp_mm_geq<simatrix,srmatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const simatrix_slice& M1, const simatrix& M2) {
  return spsp_mm_geq<simatrix,simatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const simatrix& M1, const srmatrix_slice& M2) {
  return spsp_mm_geq<simatrix,srmatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const simatrix& M1, const simatrix_slice& M2) {
  return spsp_mm_geq<simatrix,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const simatrix_slice& M1, const rmatrix& M2) {
  return spf_mm_geq<simatrix,rmatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const simatrix_slice& M1, const imatrix& M2) {
  return spf_mm_geq<simatrix,imatrix,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const imatrix& M1, const srmatrix_slice& M2) {
  return fsp_mm_geq<imatrix,srmatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const imatrix& M1, const simatrix_slice& M2) {
  return fsp_mm_geq<imatrix,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const imatrix_slice& M1, const srmatrix_slice& M2) {
  return fsp_mm_geq<imatrix,srmatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const imatrix_slice& M1, const simatrix_slice& M2) {
  return fsp_mm_geq<imatrix_slice,simatrix,interval>(M1,M2.A);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const simatrix_slice& M1, const rmatrix_slice& M2) {
  return spf_mm_geq<simatrix,rmatrix_slice,interval>(M1.A,M2);
}

//! Componentwise comparison of M1 and M2
inline bool operator>=(const simatrix_slice& M1, const imatrix_slice& M2) {
  return spf_mm_geq<simatrix,imatrix_slice,interval>(M1.A,M2);
}

//! Logical negation of M
inline bool operator!(const simatrix_slice& M) {
  return sp_m_not(M.A);
}

//! Standard output operator for sparse matrix slice
/**
 *  The output format is set by global flags, default is dense output.
 *  Use cout << SparseInOut; for sparse output or cout << MatrixMarketInOut; for
 *  output in matrix market format.
 */
inline std::ostream& operator<<(std::ostream& os, const simatrix_slice& M) {
  return sp_m_output<simatrix,interval>(os, M.A);
}

//! Standard input operator for sparse matrix slice
/**
 *  The input format is set by global flags, default is dense input.
 *  Use cout << SparseInOut; for sparse input or cout << MatrixMarketInOut; for
 *  input in matrix market format.
 */
inline std::istream& operator>>(std::istream& is, simatrix_slice& M) {
  simatrix tmp(M.A.m, M.A.n);
  sp_m_input<simatrix,interval>(is, tmp);
  M = tmp;
  return is;
}


//! Represents a row or column vector of a sparse matrix
/*!
    This is a helper class created by the [] operator to represent a row or column of a sparse matrix. This helper class provides read 
    and write access to the subvector using the standard operators. It is normally not necessary for the user to use this class explicitly, which
    is why the constructors are private.
 */
class simatrix_subv {
  private:
    simatrix_slice dat;
    bool row;
    int index;

    simatrix_subv(simatrix& A, bool r, int i, int j, int k, int l) : dat(A,i,j,k,l), row(r) {
       if(row) index=i; else index=k;
    }

    simatrix_subv(const simatrix& A, bool r, int i, int j, int k, int l) : dat(A,i,j,k,l), row(r) {
       if(row) index=i; else index=k;
    }

  public:
    //! Returns a reference to the i-th element of the subvector.
    /*!
       A refernce to the i-th element is returned. If this element is not explicitly stored, it is added as an
       explicit zero entry to the data structure.
    */
    interval& operator[](const int i) {
      if(row) {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb2 || i>dat.A.ub2)
          cxscthrow(ELEMENT_NOT_IN_VEC("simatrix_subv::operator[](int)"));
#endif
        return dat.element(index,i);
      } else {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb1 || i>dat.A.ub1)
          cxscthrow(ELEMENT_NOT_IN_VEC("simatrix_subv::operator[](int)"));
#endif
        return dat.element(i,index);
      }
    }

    //! Returns a copy of the i-th element of the subvector.
    /*!
       A copy to the i-th element is returned. If this element is not explicitly stored, 0 is returned
     */
    const interval operator[](const int i) const{
      if(row) {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb2 || i>dat.A.ub2)
          cxscthrow(ELEMENT_NOT_IN_VEC("simatrix_subv::operator[](int)"));
#endif
        return dat(index,i);
      } else {
#if(CXSC_INDEX_CHECK)
        if(i<dat.A.lb1 || i>dat.A.ub1)
          cxscthrow(ELEMENT_NOT_IN_VEC("simatrix_subv::operator[](int)"));
#endif
        return dat(i,index);
      }
    }

    //! Assigns v to all elements of the subvector
    simatrix_subv& operator=(const real& v) {
      return sv_vs_assign(*this,v);
    }

    //! Assigns v to all elements of the subvector
    simatrix_subv& operator=(const interval& v) {
      return sv_vs_assign(*this,v);
    }

    //! Assigns a vector to a subvector
    simatrix_subv& operator=(const srvector& v) {
      return svsp_vv_assign(*this,v);
    }

    //! Assigns a vector to a subvector
    simatrix_subv& operator=(const sivector& v) {
      return svsp_vv_assign(*this,v);
    }

    //! Assigns a vector to a subvector
    simatrix_subv& operator=(const srvector_slice& v) {
      return svsl_vv_assign(*this,v);
    }

    //! Assigns a vector to a subvector
    simatrix_subv& operator=(const sivector_slice& v) {
      return svsl_vv_assign(*this,v);
    }

    //! Assigns a vector to a subvector
    simatrix_subv& operator=(const rvector& v) {
      return svf_vv_assign(*this,v);
    }

    //! Assigns a vector to a subvector
    simatrix_subv& operator=(const ivector& v) {
      return svf_vv_assign(*this,v);
    }

    //! Assigns a vector to a subvector
    simatrix_subv& operator=(const rvector_slice& v) {
      return svf_vv_assign(*this,v);
    }

    //! Assigns a vector to a subvector
    simatrix_subv& operator=(const ivector_slice& v) {
      return svf_vv_assign(*this,v);
    }

    //! Assigns a vector to a subvector
    simatrix_subv& operator=(const srmatrix_subv& v) {
      return svsp_vv_assign(*this,srvector(v));
    }

    //! Assigns a vector to a subvector
    simatrix_subv& operator=(const simatrix_subv& v) {
      return svsp_vv_assign(*this,sivector(v));
    }

    //! Assign the componentwise product of the subvector with a scalar to the subvector
    simatrix_subv& operator*=(const real&);
    //! Assign the componentwise product of the subvector with a scalar to the subvector
    simatrix_subv& operator*=(const interval&);
    //! Assign the componentwise division of the subvector with a scalar to the subvector
    simatrix_subv& operator/=(const real&);
    //! Assign the componentwise division of the subvector with a scalar to the subvector
    simatrix_subv& operator/=(const interval&);
    //! Assign the sum of the subvector with a vector to the subvector
    simatrix_subv& operator+=(const srvector&);
    //! Assign the sum of the subvector with a vector to the subvector
    simatrix_subv& operator+=(const srvector_slice&);
    //! Assign the sum of the subvector with a vector to the subvector
    simatrix_subv& operator+=(const rvector&);
    //! Assign the sum of the subvector with a vector to the subvector
    simatrix_subv& operator+=(const rvector_slice&);
    //! Assign the difference of the subvector with a vector to the subvector
    simatrix_subv& operator-=(const srvector&);
    //! Assign the difference of the subvector with a vector to the subvector
    simatrix_subv& operator-=(const srvector_slice&);
    //! Assign the difference of the subvector with a vector to the subvector
    simatrix_subv& operator-=(const rvector&);
    //! Assign the difference of the subvector with a vector to the subvector
    simatrix_subv& operator-=(const rvector_slice&);
    //! Assign the sum of the subvector with a vector to the subvector
    simatrix_subv& operator+=(const sivector&);
    //! Assign the sum of the subvector with a vector to the subvector
    simatrix_subv& operator+=(const sivector_slice&);
    //! Assign the sum of the subvector with a vector to the subvector
    simatrix_subv& operator+=(const ivector&);
    //! Assign the sum of the subvector with a vector to the subvector
    simatrix_subv& operator+=(const ivector_slice&);
    //! Assign the difference of the subvector with a vector to the subvector
    simatrix_subv& operator-=(const sivector&);
    //! Assign the difference of the subvector with a vector to the subvector
    simatrix_subv& operator-=(const sivector_slice&);
    //! Assign the difference of the subvector with a vector to the subvector
    simatrix_subv& operator-=(const ivector&);
    //! Assign the difference of the subvector with a vector to the subvector
    simatrix_subv& operator-=(const ivector_slice&);
    //! Assign the convex hull of the subvector and a vector to the subvector
    simatrix_subv& operator|=(const srvector&);
    //! Assign the convex hull of the subvector and a vector to the subvector
    simatrix_subv& operator|=(const srvector_slice&);
    //! Assign the convex hull of the subvector and a vector to the subvector
    simatrix_subv& operator|=(const rvector&);
    //! Assign the convex hull of the subvector and a vector to the subvector
    simatrix_subv& operator|=(const rvector_slice&);
    //! Assign the convex hull of the subvector and a vector to the subvector
    simatrix_subv& operator|=(const sivector&);
    //! Assign the convex hull of the subvector and a vector to the subvector
    simatrix_subv& operator|=(const sivector_slice&);
    //! Assign the convex hull of the subvector and a vector to the subvector
    simatrix_subv& operator|=(const ivector&);
    //! Assign the convex hull of the subvector and a vector to the subvector
    simatrix_subv& operator|=(const ivector_slice&);

    friend sivector operator-(const simatrix_subv&);

    friend std::istream& operator>>(std::istream&, simatrix_subv&);

    friend int Lb(const simatrix_subv&);
    friend int Ub(const simatrix_subv&);
    friend int VecLen(const simatrix_subv&);
    friend srvector Inf(const simatrix_subv&);
    friend srvector Sup(const simatrix_subv&);

    friend class srvector;
    friend class srmatrix;
    friend class srmatrix_slice;
    friend class sivector;
    friend class simatrix;
    friend class simatrix_slice;
    friend class scivector;
    friend class scimatrix;
    friend class scimatrix_slice;

#include "vector_friend_declarations.inl"
};

//! Returns the lower index bound of the subvector
inline int Lb(const simatrix_subv& S) {
  if(S.row)
    return Lb(S.dat, 2);
  else
    return Lb(S.dat, 1);
}

//! Returns the upper index bound of the subvector
inline int Ub(const simatrix_subv& S) {
  if(S.row)
    return Ub(S.dat, 2);
  else
    return Ub(S.dat, 1);
}

//! Returns the length of the subvector
inline int VecLen(const simatrix_subv& S) {
  return Ub(S)-Lb(S)+1;
}

//! Returns the infimum of the subvector
inline srvector Inf(const simatrix_subv& S) {
  return Inf(sivector(S));
}

//! Returns the supremum of the subvector
inline srvector Sup(const simatrix_subv& S) {
  return Sup(sivector(S));
}

//! Returns the midpoint of the subvector
inline srvector mid(const simatrix_subv& S) {
  return mid(sivector(S));
}

//! Returns the diameter of the subvector
inline srvector diam(const simatrix_subv& S) {
  return diam(sivector(S));
}

//! Returns the componentwise absolute value of the subvector
inline sivector abs(const simatrix_subv& S) {
  return abs(sivector(S));
}

//! Standard output operator for subvectors
inline std::ostream& operator<<(std::ostream& os, const simatrix_subv& v) {
  os << sivector(v);
  return os;
}

//! Standard input operator for subvectors
inline std::istream& operator>>(std::istream& is, simatrix_subv& v) {
  int n = 0;
  if(v.row) n=v.dat.A.n; else n=v.dat.A.m;
  sivector tmp(n);
  is >> tmp;
  v = tmp;
  return is;
}

inline simatrix_subv simatrix::operator[](const cxscmatrix_column& c) {
#if(CXSC_INDEX_CHECK)
  if(c.col()<lb2 || c.col()>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const cxscmatrix_column&)"));
#endif
  return simatrix_subv(*this, false, lb1, ub1, c.col(), c.col());
}

inline simatrix_subv simatrix::operator[](const int i) {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || i>ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const int)"));
#endif
  return simatrix_subv(*this, true, i, i, lb2, ub2);
}

inline const simatrix_subv simatrix::operator[](const cxscmatrix_column& c) const {
#if(CXSC_INDEX_CHECK)
  if(c.col()<lb2 || c.col()>ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const cxscmatrix_column&)"));
#endif
  return simatrix_subv(*this, false, lb1, ub1, c.col(), c.col());
}

inline const simatrix_subv simatrix::operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
  if(i<lb1 || i>ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const int)"));
#endif
  return simatrix_subv(*this, true, i, i, lb2, ub2);
}

inline simatrix_subv simatrix_slice::operator[](const int i) {
#if(CXSC_INDEX_CHECK)
  if(i<A.lb1 || i>A.ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const int"));
#endif
  return simatrix_subv(*M, true, i, i, A.lb2, A.ub2);
}

inline simatrix_subv simatrix_slice::operator[](const cxscmatrix_column& c) {
#if(CXSC_INDEX_CHECK)
  if(c.col()<A.lb2 || c.col()>A.ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const cxscmatrix_column&)"));
#endif
  return simatrix_subv(*M, false, A.lb1, A.ub1, c.col(), c.col());
}

inline const simatrix_subv simatrix_slice::operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
  if(i<A.lb1 || i>A.ub1)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const int"));
#endif
  return simatrix_subv(*M, true, i, i, A.lb2, A.ub2);
}

inline const simatrix_subv simatrix_slice::operator[](const cxscmatrix_column& c) const {
#if(CXSC_INDEX_CHECK)
  if(c.col()<A.lb2 || c.col()>A.ub2)
    cxscthrow(ROW_OR_COL_NOT_IN_MAT("simatrix::operator[](const cxscmatrix_column&)"));
#endif
  return simatrix_subv(*M, false, A.lb1, A.ub1, c.col(), c.col());
}

inline sivector::sivector(const simatrix_subv& A) {
  int nnz = A.dat.A.get_nnz();
  p.reserve(nnz);
  x.reserve(nnz);

  if(A.row) {
    lb = A.dat.A.lb2;
    ub = A.dat.A.ub2;
    n = ub-lb+1; 

    for(int j=0 ; j<n ; j++) {
      for(int k=A.dat.A.p[j] ; k<A.dat.A.p[j+1] ; k++) {
        p.push_back(j);
        x.push_back(A.dat.A.x[k]);
      }
    }

  } else {
    lb = A.dat.A.lb1;
    ub = A.dat.A.ub1;
    n = ub-lb+1; 

    for(unsigned int k=0 ; k<A.dat.A.ind.size() ; k++) {
        p.push_back(A.dat.A.ind[k]);
        x.push_back(A.dat.A.x[k]);
    }
  }
}

//! Unary negation operator
inline sivector operator-(const simatrix_subv& v) {
 sivector s(v);
 return -s;
}

//! Computes the componentwise division of v1 and v2
inline sivector operator/(const simatrix_subv& v1, const real& v2) {
  return sivector(v1) / v2;
}

//! Computes the componentwise division of v1 and v2
inline sivector operator/(const simatrix_subv& v1, const interval& v2) {
  return sivector(v1) / v2;
}

//! Computes the componentwise division of v1 and v2
inline sivector operator/(const srmatrix_subv& v1, const interval& v2) {
  return srvector(v1) / v2;
}

//! Computes the componentwise product of v1 and v2
inline sivector operator*(const simatrix_subv& v1, const real& v2) {
  return sivector(v1) * v2;
}

//! Computes the componentwise product of v1 and v2
inline sivector operator*(const simatrix_subv& v1, const interval& v2) {
  return sivector(v1) * v2;
}

//! Computes the componentwise product of v1 and v2
inline sivector operator*(const srmatrix_subv& v1, const interval& v2) {
  return srvector(v1) * v2;
}

//! Computes the componentwise product of v1 and v2
inline sivector operator*(const real& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

//! Computes the componentwise product of v1 and v2
inline sivector operator*(const interval& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

//! Computes the componentwise product of v1 and v2
inline sivector operator*(const interval& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) * v2;
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) * v2;
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) * v2;
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) * v2;
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) * v2;
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) * v2;
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) * v2;
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) * v2;
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) * v2;
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) * v2;
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) * v2;
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) * v2;
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srvector& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const ivector& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const rvector& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const ivector& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 * srvector(v2);
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

//! Returns the dot product of v1 and v2
/*!
 * Note that the precision used for the computation is set by the global variable opdotprec.
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 * sivector(v2);
}

//! Returns the sum of v1 and v2
inline sivector operator+(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) + v2;
}

//! Returns the sum of v1 and v2
inline sivector operator+(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) + v2;
}

//! Returns the sum of v1 and v2
inline sivector operator+(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) + v2;
}

//! Returns the sum of v1 and v2
inline sivector operator+(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) + v2;
}

//! Returns the sum of v1 and v2
inline sivector operator+(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) + v2;
}

//! Returns the sum of v1 and v2
inline sivector operator+(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) + v2;
}

//! Returns the sum of v1 and v2
inline ivector operator+(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) + v2;
}

//! Returns the sum of v1 and v2
inline ivector operator+(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) + v2;
}

//! Returns the sum of v1 and v2
inline ivector operator+(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) + v2;
}

//! Returns the sum of v1 and v2
inline ivector operator+(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) + v2;
}

//! Returns the sum of v1 and v2
inline ivector operator+(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) + v2;
}

//! Returns the sum of v1 and v2
inline ivector operator+(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) + v2;
}

//! Returns the sum of v1 and v2
inline sivector operator+(const sivector& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

//! Returns the sum of v1 and v2
inline sivector operator+(const srvector& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

//! Returns the sum of v1 and v2
inline sivector operator+(const sivector& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

//! Returns the sum of v1 and v2
inline sivector operator+(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

//! Returns the sum of v1 and v2
inline sivector operator+(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

//! Returns the sum of v1 and v2
inline sivector operator+(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

//! Returns the sum of v1 and v2
inline ivector operator+(const ivector& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

//! Returns the sum of v1 and v2
inline ivector operator+(const rvector& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

//! Returns the sum of v1 and v2
inline ivector operator+(const ivector& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

//! Returns the sum of v1 and v2
inline ivector operator+(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 + srvector(v2);
}

//! Returns the sum of v1 and v2
inline ivector operator+(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

//! Returns the sum of v1 and v2
inline ivector operator+(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 + sivector(v2);
}

//! Returns the difference of v1 and v2
inline sivector operator-(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) - v2;
}

//! Returns the difference of v1 and v2
inline sivector operator-(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) - v2;
}

//! Returns the difference of v1 and v2
inline sivector operator-(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) - v2;
}

//! Returns the difference of v1 and v2
inline sivector operator-(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) - v2;
}

//! Returns the difference of v1 and v2
inline sivector operator-(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) - v2;
}

//! Returns the difference of v1 and v2
inline sivector operator-(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) - v2;
}

//! Returns the difference of v1 and v2
inline ivector operator-(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) - v2;
}

//! Returns the difference of v1 and v2
inline ivector operator-(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) - v2;
}

//! Returns the difference of v1 and v2
inline ivector operator-(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) - v2;
}

//! Returns the difference of v1 and v2
inline ivector operator-(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) - v2;
}

//! Returns the difference of v1 and v2
inline ivector operator-(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) - v2;
}

//! Returns the difference of v1 and v2
inline ivector operator-(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) - v2;
}

//! Returns the difference of v1 and v2
inline sivector operator-(const sivector& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

//! Returns the difference of v1 and v2
inline sivector operator-(const srvector& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

//! Returns the difference of v1 and v2
inline sivector operator-(const sivector& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

//! Returns the difference of v1 and v2
inline sivector operator-(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

//! Returns the difference of v1 and v2
inline sivector operator-(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

//! Returns the difference of v1 and v2
inline sivector operator-(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

//! Returns the difference of v1 and v2
inline ivector operator-(const ivector& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

//! Returns the difference of v1 and v2
inline ivector operator-(const rvector& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

//! Returns the difference of v1 and v2
inline ivector operator-(const ivector& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

//! Returns the difference of v1 and v2
inline ivector operator-(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 - srvector(v2);
}

//! Returns the difference of v1 and v2
inline ivector operator-(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

//! Returns the difference of v1 and v2
inline ivector operator-(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 - sivector(v2);
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const sivector& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const srvector& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const sivector& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const ivector& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const rvector& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const ivector& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 | sivector(v2);
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const srmatrix_subv& v1, const srvector& v2) {
  return srvector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const srmatrix_subv& v1, const srvector_slice& v2) {
  return srvector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const srmatrix_subv& v1, const rvector& v2) {
  return srvector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const srmatrix_subv& v1, const rvector_slice& v2) {
  return srvector(v1) | v2;
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const srvector& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

//! Returns the convex hull of v1 and v2
inline sivector operator|(const srvector_slice& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const rvector& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

//! Returns the convex hull of v1 and v2
inline ivector operator|(const rvector_slice& v1, const srmatrix_subv& v2) {
  return v1 | srvector(v2);
}

inline simatrix_subv& simatrix_subv::operator*=(const real& v) {
  *this = *this * v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator*=(const interval& v) {
  *this = *this * v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator/=(const real& v) {
  *this = *this / v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator/=(const interval& v) {
  *this = *this / v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const srvector& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const srvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const rvector& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const rvector_slice& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const srvector& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const srvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const rvector& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const rvector_slice& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const sivector& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const sivector_slice& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const ivector& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator+=(const ivector_slice& v) {
  *this = *this + v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const sivector& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const sivector_slice& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const ivector& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator-=(const ivector_slice& v) {
  *this = *this - v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const srvector& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const srvector_slice& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const rvector& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const rvector_slice& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const sivector& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const sivector_slice& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const ivector& v) {
  *this = *this | v;
  return *this;
}

inline simatrix_subv& simatrix_subv::operator|=(const ivector_slice& v) {
  *this = *this | v;
  return *this;
}

inline imatrix_subv& imatrix_subv::operator+=(const srmatrix_subv& v) {
  *this += rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator+=(const simatrix_subv& v) {
  *this += ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator+=(const srvector& v) {
  *this += rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator+=(const sivector& v) {
  *this += ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator+=(const srvector_slice& v) {
  *this += rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator+=(const sivector_slice& v) {
  *this += ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator-=(const srmatrix_subv& v) {
  *this -= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator-=(const simatrix_subv& v) {
  *this -= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator-=(const srvector& v) {
  *this -= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator-=(const sivector& v) {
  *this -= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator-=(const srvector_slice& v) {
  *this -= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator-=(const sivector_slice& v) {
  *this -= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator|=(const srmatrix_subv& v) {
  *this |= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator|=(const simatrix_subv& v) {
  *this |= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator|=(const srvector& v) {
  *this |= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator|=(const sivector& v) {
  *this |= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator|=(const srvector_slice& v) {
  *this |= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator|=(const sivector_slice& v) {
  *this |= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator&=(const srmatrix_subv& v) {
  *this &= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator&=(const simatrix_subv& v) {
  *this &= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator&=(const srvector& v) {
  *this &= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator&=(const sivector& v) {
  *this &= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator&=(const srvector_slice& v) {
  *this &= rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator&=(const sivector_slice& v) {
  *this &= ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator=(const srvector& v) {
  *this = rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator=(const sivector& v) {
  *this = ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator=(const srvector_slice& v) {
  *this = rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator=(const sivector_slice& v) {
  *this = ivector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator=(const srmatrix_subv& v) {
  *this = rvector(v);
  return *this;
}

inline imatrix_subv& imatrix_subv::operator=(const simatrix_subv& v) {
  *this = ivector(v);
  return *this;
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) == v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) == v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) == v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) == v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) == v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) == v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) == v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) == v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) == v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) == v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) == v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) == v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const sivector& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const srvector& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const sivector& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const ivector& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const rvector& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const ivector& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 == srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator==(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 == sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) != v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) != v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) != v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) != v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) != v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) != v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) != v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) != v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) != v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) != v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) != v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) != v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const sivector& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const srvector& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const sivector& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const ivector& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const rvector& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const ivector& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 != srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator!=(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 != sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) < v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) < v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) < v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) < v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) < v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) < v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) < v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) < v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const srvector& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const sivector& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const rvector& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const ivector& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 < sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const srmatrix_subv& v1, const sivector& v2) {
  return srvector(v1) <= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) <= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const srmatrix_subv& v1, const sivector_slice& v2) {
  return srvector(v1) <= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) <= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const srmatrix_subv& v1, const ivector& v2) {
  return srvector(v1) <= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) <= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const srmatrix_subv& v1, const ivector_slice& v2) {
  return srvector(v1) <= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) <= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const srvector& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const sivector& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const srvector_slice& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const rvector& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const ivector& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const rvector_slice& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator<=(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 <= sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) > v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) > v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) > v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) > v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) > v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) > v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) > v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) > v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const sivector& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const sivector& v1, const simatrix_subv& v2) {
  return v1 > sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 > sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const ivector& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const ivector& v1, const simatrix_subv& v2) {
  return v1 > sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 > srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 > sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const simatrix_subv& v1, const srvector& v2) {
  return sivector(v1) >= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const simatrix_subv& v1, const sivector& v2) {
  return sivector(v1) >= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const simatrix_subv& v1, const srvector_slice& v2) {
  return sivector(v1) >= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const simatrix_subv& v1, const sivector_slice& v2) {
  return sivector(v1) >= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const simatrix_subv& v1, const rvector& v2) {
  return sivector(v1) >= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const simatrix_subv& v1, const ivector& v2) {
  return sivector(v1) >= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const simatrix_subv& v1, const rvector_slice& v2) {
  return sivector(v1) >= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const simatrix_subv& v1, const ivector_slice& v2) {
  return sivector(v1) >= v2;
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const sivector& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const sivector& v1, const simatrix_subv& v2) {
  return v1 >= sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const sivector_slice& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const sivector_slice& v1, const simatrix_subv& v2) {
  return v1 >= sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const ivector& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const ivector& v1, const simatrix_subv& v2) {
  return v1 >= sivector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const ivector_slice& v1, const srmatrix_subv& v2) {
  return v1 >= srvector(v2);
}

//! Componentwise comparison of v1 and v2
inline bool operator>=(const ivector_slice& v1, const simatrix_subv& v2) {
  return v1 >= sivector(v2);
}

//! Logical negation operator
inline bool operator!(const simatrix_subv& x) {
  return sv_v_not(x);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const simatrix_subv& v2) {
  spsp_vv_accu<idotprecision,sivector,sivector,sparse_idot>(dot, sivector(v1), sivector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const srmatrix_subv& v2) {
  spsp_vv_accu<idotprecision,sivector,srvector,sparse_idot>(dot, sivector(v1), srvector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const simatrix_subv& v2) {
  spsp_vv_accu<idotprecision,srvector,sivector,sparse_idot>(dot, srvector(v1), sivector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const sivector& v2) {
  spsp_vv_accu<idotprecision,sivector,sivector,sparse_idot>(dot, sivector(v1), v2);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const srvector& v2) {
  spsp_vv_accu<idotprecision,sivector,srvector,sparse_idot>(dot, sivector(v1), v2);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const sivector& v2) {
  spsp_vv_accu<idotprecision,srvector,sivector,sparse_idot>(dot, srvector(v1), v2);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const sivector_slice& v2) {
  spsl_vv_accu<idotprecision,sivector,sivector_slice,sparse_idot>(dot, sivector(v1), v2);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const srvector_slice& v2) {
  spsl_vv_accu<idotprecision,sivector,srvector_slice,sparse_idot>(dot, sivector(v1), v2);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const sivector_slice& v2) {
  spsl_vv_accu<idotprecision,srvector,sivector_slice,sparse_idot>(dot, srvector(v1), v2);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const ivector& v2) {
  spf_vv_accu<idotprecision,sivector,ivector,sparse_idot>(dot, sivector(v1), v2);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const rvector& v2) {
  spf_vv_accu<idotprecision,sivector,rvector,sparse_idot>(dot, sivector(v1), v2);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const ivector& v2) {
  spf_vv_accu<idotprecision,srvector,ivector,sparse_idot>(dot, srvector(v1), v2);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const ivector_slice& v2) {
  spf_vv_accu<idotprecision,sivector,ivector_slice,sparse_idot>(dot, sivector(v1), v2);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const simatrix_subv& v1, const rvector_slice& v2) {
  spf_vv_accu<idotprecision,sivector,rvector_slice,sparse_idot>(dot, sivector(v1), v2);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srmatrix_subv& v1, const ivector_slice& v2) {
  spf_vv_accu<idotprecision,srvector,ivector_slice,sparse_idot>(dot, srvector(v1), v2);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector& v1, const simatrix_subv& v2) {
  spsp_vv_accu<idotprecision,sivector,sivector,sparse_idot>(dot, v1, sivector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector& v1, const srmatrix_subv& v2) {
  spsp_vv_accu<idotprecision,sivector,srvector,sparse_idot>(dot, v1, srvector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector& v1, const simatrix_subv& v2) {
  spsp_vv_accu<idotprecision,srvector,sivector,sparse_idot>(dot, v1, sivector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector_slice& v1, const simatrix_subv& v2) {
  slsp_vv_accu<idotprecision,sivector_slice,sivector,sparse_idot>(dot, v1, sivector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector_slice& v1, const srmatrix_subv& v2) {
  slsp_vv_accu<idotprecision,sivector_slice,srvector,sparse_idot>(dot, v1, srvector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector_slice& v1, const simatrix_subv& v2) {
  slsp_vv_accu<idotprecision,srvector_slice,sivector,sparse_idot>(dot, v1, sivector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const ivector& v1, const simatrix_subv& v2) {
  fsp_vv_accu<idotprecision,ivector,sivector,sparse_idot>(dot, v1, sivector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const ivector& v1, const srmatrix_subv& v2) {
  fsp_vv_accu<idotprecision,ivector,srvector,sparse_idot>(dot, v1, srvector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const rvector& v1, const simatrix_subv& v2) {
  fsp_vv_accu<idotprecision,rvector,sivector,sparse_idot>(dot, v1, sivector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const ivector_slice& v1, const simatrix_subv& v2) {
  fsp_vv_accu<idotprecision,ivector_slice,sivector,sparse_idot>(dot, v1, sivector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const ivector_slice& v1, const srmatrix_subv& v2) {
  fsp_vv_accu<idotprecision,ivector_slice,srvector,sparse_idot>(dot, v1, srvector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const rvector_slice& v1, const simatrix_subv& v2) {
  fsp_vv_accu<idotprecision,rvector_slice,sivector,sparse_idot>(dot, v1, sivector(v2));
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const srmatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),srvector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const sivector& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const srvector& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const sivector& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const sivector_slice& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const srvector_slice& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const sivector_slice& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const ivector& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const rvector& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const ivector& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const ivector_slice& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const simatrix_subv& v1, const rvector_slice& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,sivector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srmatrix_subv& v1, const ivector_slice& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,srvector(v1),v2);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector& v1, const srmatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector_slice& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector_slice& v1, const srmatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector_slice& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const ivector& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const ivector& v1, const srmatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const ivector_slice& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const ivector_slice& v1, const srmatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,srvector(v2));
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector_slice& v1, const simatrix_subv& v2) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,v1,sivector(v2));
  SetRe(dot, Re(dot) + tmp);
}

}  //namespace cxsc;

#include "sparsematrix.inl"

#endif 
 
