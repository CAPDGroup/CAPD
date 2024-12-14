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

/* CVS $Id: srvector.hpp,v 1.16 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_SRVECTOR_HPP_INCLUDED
#define _CXSC_SRVECTOR_HPP_INCLUDED

#include <real.hpp>
#include <intvector.hpp>
#include <rvector.hpp>
#include <intmatrix.hpp>
#include <vector>
#include <map>
#include <iostream>
#include <except.hpp>
#include <cidot.hpp>
#include <sparsedot.hpp>
#include <sparsevector.hpp>

namespace cxsc {

class srvector_slice;
class srmatrix;
class srmatrix_slice;
class srmatrix_subv;

class scvector;
class sivector;
class sivector_slice;
class scivector;

//! A sparse real vector
/*!
This data type represents a sparse real vector. Only the non zero elements are stored explicitly with their value and
the respective index. All operators are overloaded to take advantage of the sparsity.
*/
class srvector {
  private:
    std::vector<int> p;
    std::vector<real> x;
    int lb;
    int ub;
    int n; 

  public:
    //! Default constructor, creates an empty vector of size 0
    srvector() : lb(0), ub(-1) , n(0) {
    }

    //! Constructor for creating an empty vector of size s
    explicit srvector(const int s) : lb(1), ub(s), n(s) {
	p.reserve((int)(s*0.1));
	x.reserve((int)(s*0.1));
    }

    //! Constructor for creating an empty vector of size s and reserving memory for b elements
    srvector(const int s, const int b) : lb(1), ub(s), n(s) {
	p.reserve(b);
	x.reserve(b);
    }

    //! Constructor for creating a sparse vector our of a dense vector v. Only the non-zero elements of v are stored explicitly.
    srvector(const rvector& v) : lb(Lb(v)), ub(Ub(v)), n(VecLen(v)) {
        for(int i=lb ; i<=ub ; i++) {
          if(v[i] != 0.0) {
            p.push_back(i-lb);
            x.push_back(v[i]);
          }
        }
    }

    //! Creates a sparse vector of dimension n with nnz non zeros, who are defined by the elements of index and values
    srvector(const int n, const int nnz, const intvector& index, const rvector& values) : lb(1), ub(n) {
      this->n = n;
      for(int i=0 ; i<nnz ; i++) {
        if(values[Lb(values)+i] != 0.0) {
          p.push_back(index[Lb(index)+i]);
          x.push_back(values[Lb(values)+i]);
        }
      }
      
    }

    //! Creates a sparse vector of dimension n with nnz non zeros, who are defined by the elements of index and values
    srvector(const int n, const int nnz, const int *index, const real *values) : lb(1), ub(n) {
      this->n = n;
      for(int i=0 ; i<nnz ; i++) {
        if(values[i] != 0.0) {
          p.push_back(index[i]);
          x.push_back(values[i]);
        }
      } 
    }

    //! Creates a sparse vector out of a sparse vector slice
    srvector(const srvector_slice&);
    //! Creates a sparse vector out of a row or column of a sparse matrix
    srvector(const srmatrix_subv& A);

    //! Returns a reference to the STL-vector storing the row indices of the non zero elements.
    /*! This function is provided to allow easy interfacing to other software interfaces and
     *  efficient implementation of sparse algorithms. The user has to take care that the 
     *  data remains consistent (the i-th element of the STL-vector storing the indices
     *  refers to the i-th element of the STL-vector storing the values of the elements). 
     *  Note that the stored indices are always 0-based, independent of the starting index set in C-XSC. */
    std::vector<int>& row_indices() {
      return p;
    }

    //! Returns a reference to the STL-vector storing the values of the non zero elements.
    /*! This function is provided to allow easy interfacing to other software interfaces and
     *  efficient implementation of sparse algorithms. The user has to take care that the 
     *  data remains consistent (the i-th element of the STL-vector storing the indices
     *  refers to the i-th element of the STL-vector storing the values of the elements). */
    std::vector<real>& values() {
      return x;
    }

    //! Returns a reference to the STL-vector storing the row indices of the non zero elements.
    /*! This function is provided to allow easy interfacing to other software interfaces and
     *  efficient implementation of sparse algorithms. The user has to take care that the 
     *  data remains consistent (the i-th element of the STL-vector storing the indices
     *  refers to the i-th element of the STL-vector storing the values of the elements). 
     *  Note that the stored indices are always 0-based, independent of the starting index set in C-XSC. */
    const std::vector<int>& row_indices() const {
      return p;
    }

    //! Returns a reference to the STL-vector storing the values of the non zero elements.
    /*! This function is provided to allow easy interfacing to other software interfaces and
     *  efficient implementation of sparse algorithms. The user has to take care that the 
     *  data remains consistent (the i-th element of the STL-vector storing the indices
     *  refers to the i-th element of the STL-vector storing the values of the elements). */
    const std::vector<real>& values() const {
      return x;
    }

    //! Returns the number of non zero elements of this vector (note that this includes explicitly stored zeros)
    int get_nnz() const {
      return x.size();
    }

    //! Returns the density of the vector (the number of non zero elements divided by the number of elements)
    real density() const {
      return (double)x.size()/n;
    }

    //! Erases explicitly stored zeros from the data structure
    void dropzeros() {
      for(int i=0 ; i<get_nnz() ; i++) {
        if(x[i] == 0.0) {
           x.erase(x.begin()+i);
           p.erase(p.begin()+i);
        }
      }
    }

    //! Assigns v to all elements of the vector (resulting in a dense vector!)
    srvector& operator=(const real& v) {
      return sp_vs_assign<srvector,real,real>(*this,v);
    }

    //! Assign the dense vector v to a sparse vector. Only the non zero elements of v are used.
    srvector& operator=(const rvector& v) {
      return spf_vv_assign<srvector,rvector,real>(*this,v);
    }

    //! Assign the dense vector slice v to a sparse vector. Only the non zero elements of v are used.
    srvector& operator=(const rvector_slice& v) {
      return spf_vv_assign<srvector,rvector_slice,real>(*this,v);
    }

    //! Assign the sparse vector slice v to a sparse vector.
    srvector& operator=(const srvector_slice&);

    //! Returns a reference to the i-th (according to the currently used indexing) element of the vector.
    /*! If the i-th element is explicitly stored, a reference to the value is returned. If is not explicitly stored,
     *  it will be added to the data structure as a zero element. The returned reference then points to this new 
     *  element. Hence ths []-operator should only be used for write access to the elements of a sparse vector. 
     *  Use the ()-operator for read access. */
    real& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector::operator[](const int)"));
#endif
      int k;

      for(k=0 ; k<get_nnz() && p[k]<=i-lb ; k++) {
        if(p[k] == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, 0.0);

      return x[k];
    }

    //! Returns a copy of the i-th (according to the currently used indexing) element of the vector.
    /*! If the i-th element is explicitly stored, a copy to the value is returned. If is not explicitly stored,
     *  zero will be returned. This is the const-Version of this operator, added for convenience. It is suggested to use thei
     * ()-operator for read access. 
     */
    real operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector::operator[](const int)"));
#endif
      return (*this)(i);
    }

    //! Returns a copy of the i-th (according to the currently used indexing) element of the vector.
    /*! If the i-th element is explicitly stored, a copy of this value is returned. Otherwise, 0.0
     *  will be returned. Unlike with the []-operator, the data structure remains unchanged either way.
     *  Thus this operator should always be used for read-only access to the elements of a sparse vector.
     */
    const real operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector::operator()(const int)"));
#endif
      real r = 0.0;

      for(int k=0 ; k<get_nnz() && p[k]<=i-lb ; k++) {
        if(p[k] == i-lb) 
          r = x[k];
      }

      return r; 
    }

    //! Returns a slice of the vector from the i-th to the j-th (according to the currently used indexing) element.
    /*! This operator can be used for read and write access to a slice of the sparse vector.
     */
    srvector_slice operator()(const int i, const int j);

    //! Returns a vector whose elemnts are rearranged according to a given permutation vector
    /*! For a permutation vector p with p[Lb(p)+i]=j, the j-th element of the given sparse vector
     *  will be the i-th element of the returned permuted vector.
     */
    srvector operator()(const intvector& per) {
      srvector v(n,get_nnz());
      intvector pinv = perminv(per);

      std::map<int,real> work;
      for(int i=0 ; i<get_nnz() ; i++)
         work.insert(std::make_pair(pinv[Lb(pinv)+p[i]], x[i]));
 
      for(std::map<int,real>::iterator it=work.begin() ; it!=work.end() ; it++) {
         v.p.push_back(it->first);
         v.x.push_back(it->second);
      }

      return v;
    }

    //! Returns a vector whose elemnts are rearranged according to a given permutation matrix
    /*! 
     * For a given sparse vector x and permutation matrix P, this operator return the result of 
     * P*x. This operator is more efficient than explicitly computing P*x, since P is transformed
     * into a permutation vector and the permutation ist performed directly, without explicitly
     * computing P*x.
     */
    srvector operator()(const intmatrix& P) {
      intvector p = permvec(P);
      return (*this)(p);
    }

    //! Operator for multiplication with a scalar, result is assigned to the vector
    srvector& operator*=(const real& s) {
      return sp_vs_multassign(*this,s);
    }

    //! Operator for division of each element of the vector with a scalar, result is assigned to the vector
    srvector& operator/=(const real& s) {
      return sp_vs_divassign(*this,s);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    srvector& operator+=(const rvector& v) 
    {
      return spf_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    srvector& operator+=(const rvector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    srvector& operator+=(const srvector& v) {
      return spsp_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    srvector& operator+=(const srvector_slice&);

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    srvector& operator-=(const rvector& v) {
      return spf_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    srvector& operator-=(const rvector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    srvector& operator-=(const srvector& v) {
      return spsp_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    srvector& operator-=(const srvector_slice&);

    friend int Lb(const srvector&);
    friend int Ub(const srvector&);
    friend void SetLb(srvector&, const int);
    friend void SetUb(srvector&, const int);

    friend int VecLen(const srvector&);
    friend srvector Re(const scvector&);
    friend srvector Im(const scvector&);
    friend srvector Inf(const sivector&);
    friend srvector Sup(const sivector&);
    friend srvector InfRe(const scivector&);
    friend srvector SupRe(const scivector&);
    friend srvector InfIm(const scivector&);
    friend srvector SupIm(const scivector&);
    friend srvector mid(const sivector&);
    friend srvector diam(const sivector&);
    friend srvector absmin(const sivector&);
    friend srvector absmax(const sivector&);
    friend srvector abs(const srvector&);
    friend srvector mid(const sivector_slice&);
    friend srvector diam(const sivector_slice&);


    friend class srvector_slice;
    friend class scvector_slice;
    friend class scvector;
    friend class sivector_slice;
    friend class sivector;
    friend class scivector_slice;
    friend class scivector;
    friend class srmatrix_subv;
    friend class rvector;
    friend class rvector_slice;
    friend class ivector;
    friend class ivector_slice;
    friend class cvector;
    friend class cvector_slice;
    friend class civector;
    friend class civector_slice;


#include "vector_friend_declarations.inl"
};

inline rvector::rvector(const srvector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new real[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline rvector& rvector::operator=(const srvector& v) {
  return fsp_vv_assign<rvector,srvector,real>(*this,v);
}

inline rvector& rvector::operator=(const srvector_slice& v) {
  return fsl_vv_assign<rvector,srvector_slice,real>(*this,v);
}


//! Sets the lower index bound of the vector v to i
/** 
 * After setting the lower index bound to i, the indexing of the vector is i-based.
 */
inline void SetLb(srvector& v, const int i) {
  v.lb = i;
  v.ub = v.lb + v.n - 1;
}

//! Sets the upper index bound of the vector v to i
/** 
 * After setting the upper index bound to i, the indexing of the vector of dimension n is (i-n+1)-based.
 */
inline void SetUb(srvector& v, const int j) {
  v.ub = j;
  v.lb = v.ub - v.n + 1;
}

//! Returns the lower index bound of the vector v
inline int Lb(const srvector& v) {
  return v.lb;
}

//! Returns the upper index bound of the vector v
inline int Ub(const srvector& v) {
  return v.ub;
}

//! Returns the length of the vector (the dimension)
inline int VecLen(const srvector& v) {
  return v.n;
}

//! Resizes the vector to length 0 (all elements are deleted)
inline void Resize(srvector& v) {
  sp_v_resize(v);
}

//! Resizes the vector to length n.
/**
 * All elements of the vector that can still be stored after the resizing are copied into the resized vector.
 */
inline void Resize(srvector& v, const int n) {
  sp_v_resize(v,n);
}

//! Resizes the vector to length u-l+1.
/**
 * The new vector has lower index bound l and upper index bound u. 
 * All elements of the vector that can still be stored after the resizing are copied into the resized vector.
 */
inline void Resize(srvector& v, const int l, const int u) {
  sp_v_resize(v,l,u);
}

//! Unary operator, returns -v
inline srvector operator-(const srvector& v) {
  return sp_v_negative(v);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline real operator*(const srvector& v1, const rvector& v2) {
  return spf_vv_mult<srvector,rvector,real,sparse_dot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline real operator*(const rvector& v1, const srvector& v2) {
  return fsp_vv_mult<rvector,srvector,real,sparse_dot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline real operator*(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_mult<srvector,rvector_slice,real,sparse_dot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline real operator*(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_mult<rvector_slice,srvector,real,sparse_dot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline real operator*(const srvector& v1, const srvector& v2) {
  return spsp_vv_mult<srvector,srvector,real,sparse_dot>(v1,v2);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline srvector operator*(const srvector& v, const real& s) {
  return sp_vs_mult<srvector,real,srvector>(v,s);
}

//! Divides all elements of v by the scalar s and returns the result as a new vector
inline srvector operator/(const srvector& v, const real& s) {
  return sp_vs_div<srvector,real,srvector>(v,s);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline srvector operator*(const real& s, const srvector& v) {
  return sp_sv_mult<real,srvector,srvector>(s,v);
}

//! Element-wise addition of the vectors v1 and v2
inline rvector operator+(const rvector& v1, const srvector& v2) {
  return fsp_vv_add<rvector,srvector,rvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline rvector operator+(const srvector& v1, const rvector& v2) {
  return spf_vv_add<srvector,rvector,rvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline rvector operator+(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_add<rvector_slice,srvector,rvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline rvector operator+(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_add<srvector,rvector_slice,rvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline srvector operator+(const srvector& v1, const srvector& v2) {
  return spsp_vv_add<srvector,srvector,srvector,real>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline rvector operator-(const rvector& v1, const srvector& v2) {
  return fsp_vv_sub<rvector,srvector,rvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline rvector operator-(const srvector& v1, const rvector& v2) {
  return spf_vv_sub<srvector,rvector,rvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline rvector operator-(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_sub<rvector_slice,srvector,rvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline rvector operator-(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_sub<srvector,rvector_slice,rvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline srvector operator-(const srvector& v1, const srvector& v2) {
  return spsp_vv_sub<srvector,srvector,srvector,real>(v1,v2);
}

inline rvector& rvector::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline rvector_slice& rvector_slice::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}
 
inline rvector& rvector::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline rvector_slice& rvector_slice::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const srvector& v1, const srvector& v2) {
  return spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const srvector& v1, const rvector& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const rvector& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector& v1, const srvector& v2) {
  return !spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector& v1, const rvector& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const rvector& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector& v1, const rvector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const rvector_slice& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector& v1, const srvector& v2) {
  return spsp_vv_less<srvector,srvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector& v1, const rvector& v2) {
  return spf_vv_less<srvector,rvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const rvector& v1, const srvector& v2) {
  return fsp_vv_less<rvector,srvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_less<srvector,rvector_slice,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_less<rvector_slice,srvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector& v1, const srvector& v2) {
  return spsp_vv_leq<srvector,srvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector& v1, const rvector& v2) {
  return spf_vv_leq<srvector,rvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const rvector& v1, const srvector& v2) {
  return fsp_vv_leq<rvector,srvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_leq<srvector,rvector_slice,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_leq<rvector_slice,srvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const srvector& v1, const srvector& v2) {
  return spsp_vv_greater<srvector,srvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const srvector& v1, const rvector& v2) {
  return spf_vv_greater<srvector,rvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const rvector& v1, const srvector& v2) {
  return fsp_vv_greater<rvector,srvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_greater<srvector,rvector_slice,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_greater<rvector_slice,srvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const srvector& v1, const srvector& v2) {
  return spsp_vv_geq<srvector,srvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const srvector& v1, const rvector& v2) {
  return spf_vv_geq<srvector,rvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const rvector& v1, const srvector& v2) {
  return fsp_vv_geq<rvector,srvector,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_geq<srvector,rvector_slice,real>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_geq<rvector_slice,srvector,real>(v1,v2);
}

//! Unary logical negation of x
/** 
 * Returns true only if all elements of x are not equal to zero.
 */
inline bool operator!(const srvector& x) {
  return sp_v_not(x);
}

//! Output operator for sparse vector v.
/**
 *  The output format is set by global flags, default is dense output.
 *  Use cout << SparseInOut; for sparse output or cout << MatrixMarketInOut; for
 *  output in matrix market format.
 */
inline std::ostream& operator<<(std::ostream& os, const srvector& v) {
  return sp_v_output<srvector,real>(os,v);
}

//! Input operator for sparse vector v.
/**
 *  The input format is set by global flags, default is dense input.
 *  Use cout << SparseInOut; for sparse input or cout << MatrixMarketInOut; for
 *  input in matrix market format.
 */
inline std::istream& operator>>(std::istream& is, srvector& v) {
  return sp_v_input<srvector,real>(is,v);
}

//! Helper class for slices of sparse vectors.
/**
 * This class stores a reference to a sparse vector and operates on a slice of it. This class
 * is used internally by C-XSC, it should normally not be necessary for the user to use it explicitly.
 * 
 */
class srvector_slice {
  private:
    std::vector<int>& p;
    std::vector<real>& x;
    srvector& orig;
    int start,end;
    int lb;
    int ub;
    int n;
    int nnz;
    int offset;

    //! Constructs a reference to the slice from index l to index u of the vector v
    /**
     * The current indexing of v is used for the index parameters l and u. All read an write
     * operations with the slice directly reference the data of the vector v.
     */
    srvector_slice(srvector& v, int l, int u) : p(v.p), x(v.x), orig(v), lb(l), ub(u), n(u-l+1)  {
      int i;

      for(i=0 ; i<v.get_nnz() && p[i]<lb-v.lb ; i++);

      start = i;

      for(i=start ; i<v.get_nnz() && p[i]<=ub-v.lb ; i++);

      end = i-1;

      nnz = end-start+1;
      offset = lb-v.lb;
    }

  public:
    //! Returns the number of non zero elements of this vector slice (note that this includes explicitly stored zeros)
    int get_nnz() const {
      return nnz;
    }
    //! Returns the density of the vector slice (the number of non zero elements divided by the number of elements)
    real density() const {
      return (double)nnz/n;
    }

    //! Returns a reference to the i-th (according to the currently used indexing) element of the vector slice.
    /*! If the i-th element is explicitly stored, a reference to the value is returned. If is not explicitly stored,
     *  it will be added to the data structure as a zero element. The returned reference then points to this new 
     *  element. Hence the []-operator should only be used for write access to the elements of a sparse vector. 
     *  Use the ()-operator for read access. */
    real& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector_slice::operator[](const int)"));
#endif
      int k;

      for(k=start ; k<end+1 && p[k]-start<=i-lb ; k++) {
        if(p[k]-offset == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, 0.0);
      end++;

      return x[k];
    }

    //! Returns a copy of the i-th (according to the currently used indexing) element of the vector slice.
    /*! If the i-th element is explicitly stored, a copy to the value is returned. If is not explicitly stored,
     *  zero will be returned. This is the const-version of this operator, added for convenience. It is suggested to use the
     * ()-operator for read access. 
     */
    real operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector_slice::operator[](const int)"));
#endif
      return (*this)(i);
    }

    //! Returns a copy of the i-th (according to the currently used indexing) element of the vector slice.
    /*! If the i-th element is explicitly stored, a copy of this value is returned. Otherwise, 0.0
     *  will be returned. Unlike with the []-operator, the data structure remains unchanged either way.
     *  Thus this operator should always be used for read-only access to the elements of a sparse vector slice.
     */
    const real operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector_slice::operator()(const int)"));
#endif
      real r = 0.0;

      for(int k=start ; k<end && p[k]-start<=i-lb ; k++) {
        if(p[k]-start == i-lb) 
          r = x[k];
      }

      return r; 
    }

    //! Assigns v to all elements of the vector slice
    srvector_slice& operator=(const real& v) {
      return sl_vs_assign<srvector_slice,real,real,std::vector<real>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    srvector_slice& operator=(const srvector_slice& v) {
      return slsl_vv_assign<srvector_slice,srvector_slice,real,std::vector<real>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    srvector_slice& operator=(const srvector& v) {
      return slsp_vv_assign<srvector_slice,srvector,real,std::vector<real>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    srvector_slice& operator=(const rvector& v) {
      return slf_vv_assign<srvector_slice,rvector,real,std::vector<real>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    srvector_slice& operator=(const rvector_slice& v) {
      return slf_vv_assign<srvector_slice,rvector,real,std::vector<real>::iterator>(*this,v);
    }

    //! Operator for multiplication with a scalar, result is assigned to the vector slice
    srvector_slice& operator*=(const real& s) {
      return sl_vs_multassign(*this,s);
    }

    //! Operator for division of each element of the vector slice with a scalar, result is assigned to the vector slice
    srvector_slice& operator/=(const real& s) {
      return sl_vs_divassign(*this,s);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    srvector_slice& operator+=(const rvector& v) {
      return slf_vv_addassign<srvector_slice,rvector,real>(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    srvector_slice& operator+=(const rvector_slice& v) {
      return slf_vv_addassign<srvector_slice,rvector_slice,real>(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    srvector_slice& operator+=(const srvector& v) {
      return slsp_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    srvector_slice& operator+=(const srvector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    srvector_slice& operator-=(const rvector& v) {
      return slf_vv_subassign<srvector_slice,rvector,real>(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    srvector_slice& operator-=(const rvector_slice& v) {
      return slf_vv_subassign<srvector_slice,rvector_slice,real>(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    srvector_slice& operator-=(const srvector& v) {
      return slsp_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    srvector_slice& operator-=(const srvector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }


    friend int Lb(const srvector_slice&);
    friend int Ub(const srvector_slice&);
    friend int VecLen(const srvector_slice&);

    friend srvector operator*(const srmatrix&, const srvector_slice&); //ok
    friend srvector operator*(const srmatrix_slice&, const srvector_slice&); //ok

    friend class srvector;
    friend class scvector;
    friend class sivector;
    friend class scivector;
    friend class srmatrix_subv;
    friend class rvector;
    friend class rvector_slice;
    friend class ivector;
    friend class ivector_slice;
    friend class cvector;
    friend class cvector_slice;
    friend class civector;
    friend class civector_slice;

#include "vector_friend_declarations.inl"
};

inline rvector::rvector(const srvector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new real[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline rvector_slice& rvector_slice::operator=(const srvector& v) {
  *this = rvector(v);
  return *this;
}

inline rvector_slice& rvector_slice::operator=(const srvector_slice& v) {
  *this = rvector(v);
  return *this;
}

inline srvector::srvector(const srvector_slice& s) : lb(s.lb), ub(s.ub), n(s.n)  {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(s.x[i]);
  }

}

inline srvector& srvector::operator=(const srvector_slice& v) {
  return spsl_vv_assign<srvector,srvector_slice,real>(*this,v);
}

inline srvector_slice srvector::operator()(const int i, const int j) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || j>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector::operator()(const int,const int)"));
#endif
  return srvector_slice(*this,i,j);
}

//! Returns the vector -v
inline srvector operator-(const srvector_slice& v) {
  return sl_v_negative<srvector_slice,srvector>(v);
}

//! Returns the lower index bound of the vector slice v
inline int Lb(const srvector_slice& v) {
  return v.lb;
}

//! Returns the upper index bound of the vector slice v
inline int Ub(const srvector_slice& v) {
  return v.ub;
}

//! Returns the length of the vector slice
inline int VecLen(const srvector_slice& v) {
  return v.n;
}

//! Returns the vector whose elements are the respective absolute values of the elements of v
inline srvector abs(const srvector& v) {
  srvector ret(v);
  std::vector<real>& x = ret.values();
  for(unsigned int i=0 ; i<x.size() ; i++)
    x[i] = abs(x[i]);
  return ret;
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline real operator*(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_mult<srvector_slice,rvector,real,sparse_dot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */inline real operator*(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_mult<rvector,srvector_slice,real,sparse_dot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */inline real operator*(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_mult<srvector_slice,rvector_slice,real,sparse_dot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */inline real operator*(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_mult<rvector_slice,srvector_slice,real,sparse_dot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */inline real operator*(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_mult<srvector,srvector_slice,real,sparse_dot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */inline real operator*(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_mult<srvector_slice,srvector,real,sparse_dot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */inline real operator*(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_mult<srvector_slice,srvector_slice,real,sparse_dot>(v1,v2);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline srvector operator*(const srvector_slice& v, const real& s) {
  return sp_vs_mult<srvector_slice,real,srvector>(v,s);
}

//! Divides all elements of v by the scalar s and returns the result as a new vector
inline srvector operator/(const srvector_slice& v, const real& s) {
  return sp_vs_div<srvector_slice,real,srvector>(v,s);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline srvector operator*(const real& s, const srvector_slice& v) {
  return sp_sv_mult<real,srvector_slice,srvector>(s,v);
}

//! Element-wise addition of v1 and v2
inline rvector operator+(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_add<rvector,srvector_slice,rvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline rvector operator+(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_add<srvector_slice,rvector,rvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline rvector operator+(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_add<rvector_slice,srvector_slice,rvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline rvector operator+(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_add<srvector_slice,rvector_slice,rvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline srvector operator+(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_add<srvector_slice,srvector_slice,srvector,real>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline srvector operator+(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_add<srvector,srvector_slice,srvector,real>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline srvector operator+(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_add<srvector_slice,srvector,srvector,real>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline rvector operator-(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_sub<rvector,srvector_slice,rvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline rvector operator-(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_sub<srvector_slice,rvector,rvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline rvector operator-(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_sub<rvector_slice,srvector_slice,rvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline rvector operator-(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_sub<srvector_slice,rvector_slice,rvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline srvector operator-(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_sub<srvector_slice,srvector_slice,srvector,real>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline srvector operator-(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_sub<srvector,srvector_slice,srvector,real>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline srvector operator-(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_sub<srvector_slice,srvector,srvector,real>(v1,v2);
}

inline rvector& rvector::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline rvector_slice& rvector_slice::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline srvector& srvector::operator+=(const srvector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline rvector& rvector::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline rvector_slice& rvector_slice::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline srvector& srvector::operator-=(const srvector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector_slice& v1, const srvector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector_slice& v1, const rvector& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const rvector& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector_slice& v1, const srvector& v2) {
  return !slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector& v1, const srvector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector_slice& v1, const rvector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const rvector_slice& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_less<srvector_slice,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_less<srvector_slice,srvector,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_less<srvector,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_less<srvector_slice,rvector,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_less<rvector,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_less<srvector_slice,rvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_less<rvector_slice,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_leq<srvector_slice,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_leq<srvector_slice,srvector,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_leq<srvector,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_leq<srvector_slice,rvector,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_leq<rvector,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_leq<srvector_slice,rvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_leq<rvector_slice,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_greater<srvector_slice,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_greater<srvector_slice,srvector,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_greater<srvector,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_greater<srvector_slice,rvector,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_greater<rvector,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_greater<srvector_slice,rvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_greater<rvector_slice,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_geq<srvector_slice,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_geq<srvector_slice,srvector,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_geq<srvector,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_geq<srvector_slice,rvector,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_geq<rvector,srvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_geq<srvector_slice,rvector_slice,real>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_geq<rvector_slice,srvector_slice,real>(v1,v2);
}

//! Unary logical negation of x
/** 
 * Returns true only if all elements of x are not equal to zero.
 */
inline bool operator!(const srvector_slice& x) {
  return sl_v_not(x);
}

//! Output operator for sparse vector slice v.
/**
 *  The output format is set by global flags, default is dense output.
 *  Use cout << SparseInOut; for sparse output or cout << MatrixMarketInOut; for
 *  output in matrix market format.
 */
inline std::ostream& operator<<(std::ostream& os, const srvector_slice& v) {
  return sl_v_output<srvector_slice, real>(os,v);
}

//! Input operator for sparse vector slice v.
/**
 *  The input format is set by global flags, default is dense input.
 *  Use cout << SparseInOut; for sparse input or cout << MatrixMarketInOut; for
 *  input in matrix market format.
 */
inline std::istream& operator>>(std::istream& is, srvector_slice& v) {
  return sl_v_input<srvector_slice, real>(is,v);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(dotprecision& dot, const srvector& x, const srvector& y) {
  spsp_vv_accu<dotprecision,srvector,srvector,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(dotprecision& dot, const srvector& x, const rvector& y) {
  spf_vv_accu<dotprecision,srvector,rvector,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(dotprecision& dot, const srvector& x, const rvector_slice& y) {
  spf_vv_accu<dotprecision,srvector,rvector_slice,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(dotprecision& dot, const rvector& x, const srvector& y) {
  fsp_vv_accu<dotprecision,rvector,srvector,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(dotprecision& dot, const rvector_slice& x, const srvector& y) {
  fsp_vv_accu<dotprecision,rvector_slice,srvector,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(dotprecision& dot, const srvector_slice& x, const rvector& y) {
  slf_vv_accu<dotprecision,srvector_slice,rvector,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(dotprecision& dot, const srvector_slice& x, const rvector_slice& y) {
  slf_vv_accu<dotprecision,srvector_slice,rvector_slice,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(dotprecision& dot, const rvector& x, const srvector_slice& y) {
  fsl_vv_accu<dotprecision,rvector,srvector_slice,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(dotprecision& dot, const rvector_slice& x, const srvector_slice& y) {
  fsl_vv_accu<dotprecision,rvector_slice,srvector_slice,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(dotprecision& dot, const srvector_slice& x, const srvector_slice& y) {
  slsl_vv_accu<dotprecision,srvector_slice,srvector_slice,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(dotprecision& dot, const srvector& x, const srvector_slice& y) {
  spsl_vv_accu<dotprecision,srvector,srvector_slice,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(dotprecision& dot, const srvector_slice& x, const srvector& y) {
  slsp_vv_accu<dotprecision,srvector_slice,srvector,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(dotprecision& dot, const srvector& x, const srvector& y) {
  spsp_vv_accuapprox<dotprecision,srvector,srvector,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(dotprecision& dot, const srvector& x, const rvector& y) {
  spf_vv_accuapprox<dotprecision,srvector,rvector,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(dotprecision& dot, const srvector& x, const rvector_slice& y) {
  spf_vv_accuapprox<dotprecision,srvector,rvector_slice,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(dotprecision& dot, const rvector& x, const srvector& y) {
  fsp_vv_accuapprox<dotprecision,rvector,srvector,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(dotprecision& dot, const rvector_slice& x, const srvector& y) {
  fsp_vv_accuapprox<dotprecision,rvector_slice,srvector,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(dotprecision& dot, const srvector_slice& x, const rvector& y) {
  slf_vv_accuapprox<dotprecision,srvector_slice,rvector,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(dotprecision& dot, const srvector_slice& x, const rvector_slice& y) {
  slf_vv_accuapprox<dotprecision,srvector_slice,rvector_slice,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(dotprecision& dot, const rvector& x, const srvector_slice& y) {
  fsl_vv_accuapprox<dotprecision,rvector,srvector_slice,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(dotprecision& dot, const rvector_slice& x, const srvector_slice& y) {
  fsl_vv_accuapprox<dotprecision,rvector_slice,srvector_slice,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(dotprecision& dot, const srvector_slice& x, const srvector_slice& y) {
  slsl_vv_accuapprox<dotprecision,srvector_slice,srvector_slice,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(dotprecision& dot, const srvector& x, const srvector_slice& y) {
  spsl_vv_accuapprox<dotprecision,srvector,srvector_slice,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(dotprecision& dot, const srvector_slice& x, const srvector& y) {
  slsp_vv_accuapprox<dotprecision,srvector_slice,srvector,sparse_dot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const rvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const rvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector_slice& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector_slice& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const rvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const rvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const rvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const rvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector_slice& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector_slice& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const rvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const rvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const rvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const rvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const rvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const rvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate_approx(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector_slice& x, const rvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector_slice& x, const rvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector_slice& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector& x, const srvector_slice& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector_slice& x, const srvector& y) {
  dotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}


} //namespace cxsc

#include "sparsevector.inl"

#endif
