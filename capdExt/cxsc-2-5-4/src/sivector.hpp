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

/* CVS $Id: sivector.hpp,v 1.15 2014/01/30 17:23:48 cxsc Exp $ */

#ifndef _CXSC_SIVECTOR_HPP_INCLUDED
#define _CXSC_SIVECTOR_HPP_INCLUDED

#include <interval.hpp>
#include <ivector.hpp>
#include <vector>
#include <iostream>
#include <cidot.hpp>
#include <srvector.hpp>
#include <sparseidot.hpp>
#include <sparsevector.hpp>

namespace cxsc {

class srvector_slice;
class srmatrix;
class srmatrix_slice;
class srmatrix_subv;
class sivector_slice;
class simatrix;
class simatrix_slice;
class simatrix_subv;
class scivector;
class scivector_slice;
class scimatrix;
class scimatrix_slice;
class scimatrix_subv;

//! A sparse interval vector
/*!
This data type represents a sparse interval vector. Only the non zero elements are stored explicitly with their value and
the respective index. All operators are overloaded to take advantage of the sparsity.
*/
class sivector {
  private:
    std::vector<int> p;
    std::vector<interval> x;
    int lb;
    int ub;
    int n; 

  public:
    //! Default constructor, creates an empty vector of size 0
    sivector() : lb(0), ub(-1) , n(0) {
    }

    //! Constructor for creating an empty vector of size s
    explicit sivector(const int s) : lb(1), ub(s), n(s) {
	p.reserve((int)(s*0.1));
	x.reserve((int)(s*0.1));
    }

    //! Constructor for creating an empty vector of size s and reserving memory for b elements
    sivector(const int s, const int b) : lb(1), ub(s), n(s) {
	p.reserve(b);
	x.reserve(b);
    }

    //! Constructor for creating a sparse vector our of a dense vector v. Only the non-zero elements of v are stored explicitly.
    sivector(const ivector& v) : lb(Lb(v)), ub(Ub(v)), n(VecLen(v)) {
        for(int i=lb ; i<=ub ; i++) {
          if(v[i] != 0.0) {
            p.push_back(i-lb);
            x.push_back(v[i]);
          }
        }
    }

    //! Constructor for creating a sparse vector our of a dense vector v. Only the non-zero elements of v are stored explicitly.
    sivector(const rvector& v) : lb(Lb(v)), ub(Ub(v)), n(VecLen(v)) {
        for(int i=lb ; i<=ub ; i++) {
          if(v[i] != 0.0) {
            p.push_back(i-lb);
            x.push_back(interval(v[i]));
          }
        }
    }

    //! Creates a sparse vector of dimension n with nnz non zeros, who are defined by the elements of index and values
    sivector(const int n, const int nnz, const intvector& index, const ivector& values) : lb(1), ub(n) {
      this->n = n;
      for(int i=0 ; i<nnz ; i++) {
        if(values[Lb(values)+i] != 0.0) {
          p.push_back(index[Lb(index)+i]);
          x.push_back(values[Lb(values)+i]);
        }
      }
    }

    //! Creates a sparse vector of dimension n with nnz non zeros, who are defined by the elements of index and values
    sivector(const int n, const int nnz, const int* index, const interval* values) : lb(1), ub(n) {
      this->n = n;
      for(int i=0 ; i<nnz ; i++) {
        if(values[i] != 0.0) {
          p.push_back(index[i]);
          x.push_back(values[i]);
        }
      }
    }

    //! Creates a sparse interval vector out of a sparse real vector
    sivector(const srvector& v) : p(v.p), lb(v.lb), ub(v.ub), n(v.n) {
      x.reserve(v.get_nnz());
      for(int i=0 ; i<v.get_nnz() ; i++) 
        x.push_back(interval(v.x[i]));
    }

    //! Creates a sparse interval vector out of a sparse real vector slice
    sivector(const srvector_slice&);
    //! Creates a sparse interval vector out of a sparse interval vector slice    
    sivector(const sivector_slice&);
    //! Creates a sparse interval vector out of a row or column of a sparse real matrix
    sivector(const srmatrix_subv& A);
    //! Creates a sparse interval vector out of a row or column of a sparse interval matrix
    sivector(const simatrix_subv& A);

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
    std::vector<interval>& values() {
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
    const std::vector<interval>& values() const {
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


    /* sivector& operator=(const sivector& v) {
      p = v.p;
      x = v.x;
      return *this;
    } */

    //! Assign a sparse real vector to a sparse interval vector
    sivector& operator=(const srvector& v) {
      n = v.n;
      p = v.p;
      x.clear();
      x.reserve(v.get_nnz());
      for(unsigned int i=0 ; i<v.x.size() ; i++)
        x[i] = interval(v.x[i]);
      return *this;
    } 

    //! Assigns v to all elements of the vector (resulting in a dense vector!)
    sivector& operator=(const real& v) {
      return sp_vs_assign<sivector,real,interval>(*this,v);
    }

    //! Assigns v to all elements of the vector (resulting in a dense vector!)
    sivector& operator=(const interval& v) {
      return sp_vs_assign<sivector,interval,interval>(*this,v);
    }

    //! Assign the dense vector v to a sparse vector. Only the non zero elements of v are used.
    sivector& operator=(const rvector& v) {
      return spf_vv_assign<sivector,rvector,interval>(*this,v);
    }

    //! Assign the dense vector v to a sparse vector. Only the non zero elements of v are used.
    sivector& operator=(const ivector& v) {
      return spf_vv_assign<sivector,ivector,interval>(*this,v);
    }

    //! Assign the dense vector slice v to a sparse vector. Only the non zero elements of v are used.
    sivector& operator=(const rvector_slice& v) {
      return spf_vv_assign<sivector,rvector_slice,interval>(*this,v);
    }

    //! Assign the dense vector slice v to a sparse vector. Only the non zero elements of v are used.
    sivector& operator=(const ivector_slice& v) {
      return spf_vv_assign<sivector,ivector_slice,interval>(*this,v);
    }

    //! Assign the sparse vector slice v to a sparse vector.
    sivector& operator=(const sivector_slice&);
    //! Assign the sparse vector slice v to a sparse vector.
    sivector& operator=(const srvector_slice&);

    //! Returns a reference to the i-th (according to the currently used indexing) element of the vector.
    /*! If the i-th element is explicitly stored, a reference to the value is returned. If is not explicitly stored,
     *  it will be added to the data structure as a zero element. The returned reference then points to this new 
     *  element. Hence ths []-operator should only be used for write access to the elements of a sparse vector. 
     *  Use the ()-operator for read access. */
    interval& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("sivector::operator[](const int)"));
#endif
      int k;

      for(k=0 ; k<get_nnz() && p[k]<=i-lb ; k++) {
        if(p[k] == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, interval(0.0));

      return x[k];
    }

    //! Returns a copy of the i-th (according to the currently used indexing) element of the vector.
    /*! If the i-th element is explicitly stored, a copy to the value is returned. If is not explicitly stored,
     *  zero will be returned. This is the const-Version of this operator, added for convenience. It is suggested to use thei
     * ()-operator for read access. 
     */
    interval operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("sivector::operator[](const int)"));
#endif
      return (*this)(i);
    }

    //! Returns a copy of the i-th (according to the currently used indexing) element of the vector.
    /*! If the i-th element is explicitly stored, a copy of this value is returned. Otherwise, 0.0
     *  will be returned. Unlike with the []-operator, the data structure remains unchanged either way.
     *  Thus this operator should always be used for read-only access to the elements of a sparse vector.
     */
    const interval operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("sivector::operator()(const int)"));
#endif
      interval r(0.0);

      for(int k=0 ; k<get_nnz() && p[k]<=i-lb ; k++) {
        if(p[k] == i-lb) 
          r = x[k];
      }

      return r; 
    }

    //! Returns a vector whose elemnts are rearranged according to a given permutation vector
    /*! For a permutation vector p with p[Lb(p)+i]=j, the j-th element of the given sparse vector
     *  will be the i-th element of the returned permuted vector.
     */
    sivector operator()(const intvector& per) {
      sivector v(n,get_nnz());
      intvector pinv = perminv(per);

      std::map<int,interval> work;
      for(int i=0 ; i<get_nnz() ; i++) {
         work.insert(std::make_pair(pinv[Lb(pinv)+p[i]], x[i]));
      }
 
      for(std::map<int,interval>::iterator it=work.begin() ; it!=work.end() ; it++) {
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
    sivector operator()(const intmatrix& P) {
      intvector p = permvec(P);
      return (*this)(p);
    }

    //! Returns a slice of the vector from the i-th to the j-th (according to the currently used indexing) element.
    /*! This operator can be used for read and write access to a slice of the sparse vector.
     */
    sivector_slice operator()(const int, const int);

    //! Operator for multiplication with a scalar, result is assigned to the vector
    sivector& operator*=(const real& s) {
      return sp_vs_multassign(*this,s);
    }

    //! Operator for multiplication with an interval, result is assigned to the vector
    sivector& operator*=(const interval& s) {
      return sp_vs_multassign(*this,s);
    }

    //! Operator for division of each element of the vector with a scalar, result is assigned to the vector
    sivector& operator/=(const real& s) {
      return sp_vs_divassign(*this,s);
    }

    //! Operator for division of each element of the vector with an interval, result is assigned to the vector
    sivector& operator/=(const interval& s) {
      return sp_vs_divassign(*this,s);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    sivector& operator+=(const rvector& v) {
      return spf_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    sivector& operator+=(const ivector& v) {
      return spf_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    sivector& operator+=(const rvector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    sivector& operator+=(const ivector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    sivector& operator+=(const srvector& v) {
      return spsp_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    sivector& operator+=(const sivector& v) {
      return spsp_vv_addassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    sivector& operator-=(const rvector& v) {
      return spf_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    sivector& operator-=(const ivector& v) {
      return spf_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    sivector& operator-=(const rvector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    sivector& operator-=(const ivector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    sivector& operator-=(const srvector& v) {
      return spsp_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    sivector& operator-=(const sivector& v) {
      return spsp_vv_subassign(*this,v);
    }

    //! Operator for element-wise convex hull with another vector, result is assigned to the vector
    sivector& operator|=(const rvector& v) {
      return spf_vv_hullassign(*this,v);
    }

    //! Operator for element-wise convex hull with another vector, result is assigned to the vector
    sivector& operator|=(const ivector& v) {
      return spf_vv_hullassign(*this,v);
    }

    //! Operator for element-wise convex hull with another vector, result is assigned to the vector
    sivector& operator|=(const rvector_slice& v) {
      return spf_vv_hullassign(*this,v);
    }

    //! Operator for element-wise convex hull with another vector, result is assigned to the vector
    sivector& operator|=(const ivector_slice& v) {
      return spf_vv_hullassign(*this,v);
    }

    //! Operator for element-wise convex hull with another vector, result is assigned to the vector
    sivector& operator|=(const srvector& v) {
      return spsp_vv_hullassign(*this,v);
    }

    //! Operator for element-wise convex hull with another vector, result is assigned to the vector
    sivector& operator|=(const sivector& v) {
      return spsp_vv_hullassign(*this,v);
    }

    //! Operator for element-wise intersection with another vector, result is assigned to the vector
    sivector& operator&=(const ivector_slice& v) {
      return spf_vv_intersectassign(*this,v);
    }

    //! Operator for element-wise intersection with another vector, result is assigned to the vector
    sivector& operator&=(const sivector& v) {
      return spsp_vv_intersectassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    sivector& operator+=(const srvector_slice&);
    //! Operator for element-wise addition with a vector, result is assigned to the vector
    sivector& operator+=(const sivector_slice&);
    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    sivector& operator-=(const srvector_slice&);
    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    sivector& operator-=(const sivector_slice&);

    friend void SetLb(sivector&, const int);
    friend void SetUb(sivector&, const int);
    friend int Lb(const sivector&);
    friend int Ub(const sivector&);
    friend srvector Inf(const sivector&);
    friend srvector Sup(const sivector&);
    friend sivector Re(const scivector&);
    friend sivector Im(const scivector&);
    friend sivector abs(const sivector&);
    friend sivector abs(const sivector_slice&);
    friend srvector mid(const sivector&);
    friend srvector diam(const sivector&);
    friend sivector abs(const scivector&);
    friend sivector abs(const scivector_slice&);
    friend srvector absmin(const sivector&);
    friend srvector absmax(const sivector&);
    friend int VecLen(const sivector&);
    friend sivector Blow(const sivector&, const real&);

    friend class srvector_slice;
    friend class sivector_slice;
    friend class scivector_slice;
    friend class scivector;
    friend class ivector;
    friend class ivector_slice;
    friend class civector;
    friend class civector_slice;

#include "vector_friend_declarations.inl"
};

inline ivector::ivector(const sivector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new interval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline ivector::ivector(const srvector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new interval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline ivector& ivector::operator=(const sivector& v) {
  return fsp_vv_assign<ivector,sivector,interval>(*this,v);
}

inline ivector& ivector::operator=(const sivector_slice& v) {
  return fsl_vv_assign<ivector,sivector_slice,interval>(*this,v);
}

inline ivector& ivector::operator=(const srvector& v) {
  return fsp_vv_assign<ivector,srvector,interval>(*this,v);
}

inline ivector& ivector::operator=(const srvector_slice& v) {
  return fsl_vv_assign<ivector,srvector_slice,interval>(*this,v);
}

//! Sets the lower index bound of the vector v to i
/** 
 * After setting the lower index bound to i, the indexing of the vector is i-based.
 */
inline void SetLb(sivector& v, const int i) {
  v.lb = i;
  v.ub = v.lb + v.n - 1;
}

//! Sets the upper index bound of the vector v to i
/** 
 * After setting the upper index bound to i, the indexing of the vector of dimension n is (i-n+1)-based.
 */
inline void SetUb(sivector& v, const int j) {
  v.ub = j;
  v.lb = v.ub - v.n + 1;
}

//! Returns the lower index bound of the vector v
inline int Lb(const sivector& v) {
  return v.lb;
}

//! Returns the upper index bound of the vector v
inline int Ub(const sivector& v) {
  return v.ub;
}

//! Resizes the vector to length 0 (all elements are deleted)
inline void Resize(sivector& v) {
  sp_v_resize(v);
}

//! Resizes the vector to length n.
/**
 * All elements of the vector that can still be stored after the resizing are copied into the resized vector.
 */
inline void Resize(sivector& v, const int n) {
  sp_v_resize(v,n);
}

//! Resizes the vector to length u-l+1.
/**
 * The new vector has lower index bound l and upper index bound u. 
 * All elements of the vector that can still be stored after the resizing are copied into the resized vector.
 */
inline void Resize(sivector& v, const int l, const int u) {
  sp_v_resize(v,l,u);
}

//! Returns the infimum of the interval vector as a new sparse point vector
inline srvector Inf(const sivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = Inf(v.x[i]);
  return res;
}

//! Returns the supremum of the interval vector as a new sparse point vector
inline srvector Sup(const sivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = Sup(v.x[i]);
  return res;
}

//! Computes the component-wise absolute values as the interval hull of \f$  \{ |v| \mid v \in [v] \} \f$ for a vector v
inline sivector abs(const sivector& v) {
  sivector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x.push_back(abs(v.x[i]));
  return res;
}

//! Computes the component-wise minimum absolute values \f$  \min\limits_{v \in [v]} (|v|) \f$ for a vector v
inline srvector absmin(const sivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x.push_back(AbsMin(v.x[i]));
  res.dropzeros();
  return res;
}

//! Computes the component-wise maximum absolute values \f$  \max\limits_{v \in [v]} (|v|) \f$ for a vector v
inline srvector absmax(const sivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x.push_back(AbsMax(v.x[i]));
  res.dropzeros();
  return res;
}

//! Compute the midpoint vector of v
inline srvector mid(const sivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++) {
    res.x.push_back(mid(v.x[i]));
  }
  return res;
}

//! Computes the diameter of v
inline srvector diam(const sivector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x.push_back(diam(v.x[i]));
  return res;
}

//! Returns the length of the vector (the dimension)
inline int VecLen(const sivector& v) {
  return v.n;
}

//! Performs an epsilon inflation of the vector v
inline sivector Blow(const sivector& v, const real& eps) {
  sivector res(v);
  for(unsigned int i=0 ; i<v.x.size() ; i++)
    res.x[i] = Blow(v.x[i],eps);
  return res;
}

//! Checks if all elements of v1 lie in the interior of v2
inline bool in (const sivector& v1, const sivector& v2) {
  for(int i=0 ; i<VecLen(v1) ; i++)
    if(!in(v1(i+Lb(v1)), v2(i+Lb(v2)))) return false;
  return true;
}

//! Checks if all elements of v1 are euqal to 0
inline bool Zero(const sivector& v1) {
  for(int i=0 ; i<VecLen(v1) ; i++)
    if(v1(i+Lb(v1)) != 0.0) return false;
  return true;
}

//! Unary operator, returns -v
inline sivector operator-(const sivector& v) {
  return sp_v_negative(v);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector& v1, const ivector& v2) {
  return spf_vv_mult<sivector,ivector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector& v1, const rvector& v2) {
  return spf_vv_mult<sivector,rvector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srvector& v1, const ivector& v2) {
  return spf_vv_mult<srvector,ivector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const rvector& v1, const sivector& v2) {
  return fsp_vv_mult<rvector,sivector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const ivector& v1, const srvector& v2) {
  return fsp_vv_mult<ivector,srvector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const ivector& v1, const sivector& v2) {
  return fsp_vv_mult<ivector,sivector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_mult<sivector,rvector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_mult<sivector,ivector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_mult<srvector,ivector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_mult<ivector_slice,srvector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_mult<ivector_slice,sivector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_mult<rvector_slice,sivector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector& v1, const srvector& v2) {
  return spsp_vv_mult<sivector,srvector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srvector& v1, const sivector& v2) {
  return spsp_vv_mult<srvector,sivector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector& v1, const sivector& v2) {
  return spsp_vv_mult<sivector,sivector,interval,sparse_idot>(v1,v2);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline sivector operator*(const sivector& v, const real& s) {
  return sp_vs_mult<sivector,real,sivector>(v,s);
}

//! Multiplies all elements of v with the interval s and returns the result as a new vector
inline sivector operator*(const sivector& v, const interval& s) {
  return sp_vs_mult<sivector,interval,sivector>(v,s);
}

//! Multiplies all elements of v with the interval s and returns the result as a new vector
inline sivector operator*(const srvector& v, const interval& s) {
  return sp_vs_mult<srvector,interval,sivector>(v,s);
}

//! Divides all elements of v by the scalar s and returns the result as a new vector
inline sivector operator/(const sivector& v, const real& s) {
  return sp_vs_div<sivector,real,sivector>(v,s);
}

//! Divides all elements of v by the interval s and returns the result as a new vector
inline sivector operator/(const sivector& v, const interval& s) {
  return sp_vs_div<sivector,interval,sivector>(v,s);
}

//! Divides all elements of v by the interval s and returns the result as a new vector
inline sivector operator/(const srvector& v, const interval& s) {
  return sp_vs_div<srvector,interval,sivector>(v,s);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline sivector operator*(const real& s, const sivector& v) {
  return sp_sv_mult<real,sivector,sivector>(s,v);
}

//! Multiplies all elements of v with the interval s and returns the result as a new vector
inline sivector operator*(const interval& s, const sivector& v) {
  return sp_sv_mult<interval,sivector,sivector>(s,v);
}

//! Multiplies all elements of v with the interval s and returns the result as a new vector
inline sivector operator*(const interval& s, const srvector& v) {
  return sp_sv_mult<interval,srvector,sivector>(s,v);
}

//! Element-wise addition of the vectors v1 and v2
inline ivector operator+(const ivector& v1, const srvector& v2) {
  return fsp_vv_add<ivector,srvector,ivector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline ivector operator+(const rvector& v1, const sivector& v2) {
  return fsp_vv_add<rvector,sivector,ivector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline ivector operator+(const ivector& v1, const sivector& v2) {
  return fsp_vv_add<ivector,sivector,ivector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline ivector operator+(const sivector& v1, const rvector& v2) {
  return spf_vv_add<sivector,rvector,ivector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline ivector operator+(const srvector& v1, const ivector& v2) {
  return spf_vv_add<srvector,ivector,ivector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline ivector operator+(const sivector& v1, const ivector& v2) {
  return spf_vv_add<sivector,ivector,ivector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline ivector operator+(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_add<ivector_slice,srvector,ivector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline ivector operator+(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_add<rvector_slice,sivector,ivector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline ivector operator+(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_add<ivector_slice,sivector,ivector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline ivector operator+(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_add<sivector,rvector_slice,ivector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline ivector operator+(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_add<srvector,ivector_slice,ivector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline ivector operator+(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_add<sivector,ivector_slice,ivector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline sivector operator+(const sivector& v1, const srvector& v2) {
  return spsp_vv_add<sivector,srvector,sivector,interval>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline sivector operator+(const srvector& v1, const sivector& v2) {
  return spsp_vv_add<srvector,sivector,sivector,interval>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline sivector operator+(const sivector& v1, const sivector& v2) {
  return spsp_vv_add<sivector,sivector,sivector,interval>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline ivector operator-(const ivector& v1, const srvector& v2) {
  return fsp_vv_sub<ivector,srvector,ivector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline ivector operator-(const rvector& v1, const sivector& v2) {
  return fsp_vv_sub<rvector,sivector,ivector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline ivector operator-(const ivector& v1, const sivector& v2) {
  return fsp_vv_sub<ivector,sivector,ivector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline ivector operator-(const sivector& v1, const rvector& v2) {
  return spf_vv_sub<sivector,rvector,ivector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline ivector operator-(const srvector& v1, const ivector& v2) {
  return spf_vv_sub<srvector,ivector,ivector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline ivector operator-(const sivector& v1, const ivector& v2) {
  return spf_vv_sub<sivector,ivector,ivector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline ivector operator-(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_sub<ivector_slice,srvector,ivector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline ivector operator-(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_sub<rvector_slice,sivector,ivector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline ivector operator-(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_sub<ivector_slice,sivector,ivector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline ivector operator-(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_sub<sivector,rvector_slice,ivector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline ivector operator-(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_sub<srvector,ivector_slice,ivector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline ivector operator-(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_sub<sivector,ivector_slice,ivector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline sivector operator-(const sivector& v1, const srvector& v2) {
  return spsp_vv_sub<sivector,srvector,sivector,interval>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline sivector operator-(const srvector& v1, const sivector& v2) {
  return spsp_vv_sub<srvector,sivector,sivector,interval>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline sivector operator-(const sivector& v1, const sivector& v2) {
  return spsp_vv_sub<sivector,sivector,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const rvector& v1, const srvector& v2) {
  return fsp_vv_hull<rvector,srvector,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const srvector& v1, const rvector& v2) {
  return spf_vv_hull<srvector,rvector,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const rvector_slice& v1, const srvector& v2) {
  return fsp_vv_hull<rvector_slice,srvector,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const srvector& v1, const rvector_slice& v2) {
  return spf_vv_hull<srvector,rvector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline sivector operator|(const srvector& v1, const srvector& v2) {
  return spsp_vv_hull<srvector,srvector,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const ivector& v1, const srvector& v2) {
  return fsp_vv_hull<ivector,srvector,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const rvector& v1, const sivector& v2) {
  return fsp_vv_hull<rvector,sivector,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const ivector& v1, const sivector& v2) {
  return fsp_vv_hull<ivector,sivector,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const sivector& v1, const rvector& v2) {
  return spf_vv_hull<sivector,rvector,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const srvector& v1, const ivector& v2) {
  return spf_vv_hull<srvector,ivector,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const sivector& v1, const ivector& v2) {
  return spf_vv_hull<sivector,ivector,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_hull<ivector_slice,srvector,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_hull<rvector_slice,sivector,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_hull<ivector_slice,sivector,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_hull<sivector,rvector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_hull<srvector,ivector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline ivector operator|(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_hull<sivector,ivector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline sivector operator|(const sivector& v1, const srvector& v2) {
  return spsp_vv_hull<sivector,srvector,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline sivector operator|(const srvector& v1, const sivector& v2) {
  return spsp_vv_hull<srvector,sivector,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of the vectors v1 and v2
inline sivector operator|(const sivector& v1, const sivector& v2) {
  return spsp_vv_hull<sivector,sivector,sivector,interval>(v1,v2);
}

//! Element-wise intersection of the vectors v1 and v2
inline sivector operator&(const ivector& v1, const sivector& v2) {
  return fsp_vv_intersect<ivector,sivector,ivector>(v1,v2);
}

//! Element-wise intersection of the vectors v1 and v2
inline sivector operator&(const sivector& v1, const ivector& v2) {
  return spf_vv_intersect<sivector,ivector,ivector>(v1,v2);
}

//! Element-wise intersection of the vectors v1 and v2
inline sivector operator&(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_intersect<ivector_slice,sivector,ivector>(v1,v2);
}

//! Element-wise intersection of the vectors v1 and v2
inline sivector operator&(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_intersect<sivector,ivector_slice,ivector>(v1,v2);
}

//! Element-wise intersection of the vectors v1 and v2
inline sivector operator&(const sivector& v1, const sivector& v2) {
  return spsp_vv_intersect<sivector,sivector,sivector,interval>(v1,v2);
}

inline ivector& ivector::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline ivector& ivector::operator+=(const sivector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator+=(const sivector& v2) {
  return fsp_vv_addassign(*this,v2);
}
 
inline ivector& ivector::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline ivector& ivector::operator-=(const sivector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator-=(const sivector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline ivector& ivector::operator|=(const srvector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline ivector& ivector::operator|=(const sivector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator|=(const srvector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator|=(const sivector& v2) {
  return fsp_vv_hullassign(*this,v2);
}

inline ivector& ivector::operator&=(const sivector& v2) {
  return fsp_vv_intersectassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator&=(const sivector& v2) {
  return fsp_vv_intersectassign(*this,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const sivector& v1, const sivector& v2) {
  return spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const sivector& v1, const srvector& v2) {
  return spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const srvector& v1, const sivector& v2) {
  return spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const sivector& v1, const rvector& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const srvector& v1, const ivector& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const sivector& v1, const ivector& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const ivector& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const rvector& v1, const sivector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const ivector& v1, const sivector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector& v1, const srvector& v2) {
  return !spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector& v1, const sivector& v2) {
  return !spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector& v1, const sivector& v2) {
  return !spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector& v1, const rvector& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector& v1, const ivector& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector& v1, const ivector& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const ivector& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const rvector& v1, const sivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const ivector& v1, const sivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector& v1, const rvector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector& v1, const ivector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector& v1, const ivector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const ivector_slice& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const rvector_slice& v1, const sivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const ivector_slice& v1, const sivector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const sivector& v1, const sivector& v2) {
  return spsp_vv_less<sivector,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector& v1, const sivector& v2) {
  return spsp_vv_less<srvector,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector& v1, const ivector& v2) {
  return spf_vv_less<srvector,ivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const sivector& v1, const ivector& v2) {
  return spf_vv_less<sivector,ivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const rvector& v1, const sivector& v2) {
  return fsp_vv_less<rvector,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const ivector& v1, const sivector& v2) {
  return fsp_vv_less<ivector,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_less<srvector,ivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_less<sivector,ivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_less<rvector_slice,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_less<ivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const sivector& v1, const sivector& v2) {
  return spsp_vv_leq<sivector,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector& v1, const sivector& v2) {
  return spsp_vv_leq<srvector,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector& v1, const ivector& v2) {
  return spf_vv_leq<srvector,ivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const sivector& v1, const ivector& v2) {
  return spf_vv_leq<sivector,ivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const rvector& v1, const sivector& v2) {
  return fsp_vv_leq<rvector,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const ivector& v1, const sivector& v2) {
  return fsp_vv_leq<ivector,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector& v1, const ivector_slice& v2) {
  return spf_vv_leq<srvector,ivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_leq<sivector,ivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const rvector_slice& v1, const sivector& v2) {
  return fsp_vv_leq<rvector_slice,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_leq<ivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector& v1, const sivector& v2) {
  return spsp_vv_greater<sivector,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector& v1, const srvector& v2) {
  return spsp_vv_greater<sivector,srvector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector& v1, const rvector& v2) {
  return spf_vv_greater<sivector,rvector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector& v1, const ivector& v2) {
  return spf_vv_greater<sivector,ivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const ivector& v1, const srvector& v2) {
  return fsp_vv_greater<ivector,srvector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const ivector& v1, const sivector& v2) {
  return fsp_vv_greater<ivector,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_greater<sivector,rvector_slice,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_greater<sivector,ivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_greater<ivector_slice,srvector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_greater<ivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector& v1, const sivector& v2) {
  return spsp_vv_geq<sivector,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector& v1, const srvector& v2) {
  return spsp_vv_geq<sivector,srvector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector& v1, const rvector& v2) {
  return spf_vv_geq<sivector,rvector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector& v1, const ivector& v2) {
  return spf_vv_geq<sivector,ivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const ivector& v1, const srvector& v2) {
  return fsp_vv_geq<ivector,srvector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const ivector& v1, const sivector& v2) {
  return fsp_vv_geq<ivector,sivector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector& v1, const rvector_slice& v2) {
  return spf_vv_geq<sivector,rvector_slice,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector& v1, const ivector_slice& v2) {
  return spf_vv_geq<sivector,ivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const ivector_slice& v1, const srvector& v2) {
  return fsp_vv_geq<ivector_slice,srvector,interval>(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const ivector_slice& v1, const sivector& v2) {
  return fsp_vv_geq<ivector_slice,sivector,interval>(v1,v2);
}

//! Output operator for sparse vector v.
/**
 *  The output format is set by global flags, default is dense output.
 *  Use cout << SparseInOut; for sparse output or cout << MatrixMarketInOut; for
 *  output in matrix market format.
 */
inline std::ostream& operator<<(std::ostream& os, const sivector& v) {
  return sp_v_output<sivector,interval>(os,v);
}

//! Input operator for sparse vector v.
/**
 *  The input format is set by global flags, default is dense input.
 *  Use cout << SparseInOut; for sparse input or cout << MatrixMarketInOut; for
 *  input in matrix market format.
 */
inline std::istream& operator>>(std::istream& is, sivector& v) {
  return sp_v_input<sivector,interval>(is,v);
}

//! Helper class for slices of sparse vectors.
/**
 * This class stores a reference to a sparse vector and operates on a slice of it. This class
 * is used internally by C-XSC, it should normally not be necessary for the user to use it explicitly.
 * 
 */
class sivector_slice {
  private:
    std::vector<int>& p;
    std::vector<interval>& x;
    sivector& orig;
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
    sivector_slice(sivector& v, int l, int u) : p(v.p), x(v.x), orig(v), lb(l), ub(u), n(u-l+1) {
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
    interval& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("sivector_slice::operator[](const int)"));
#endif
      int k;

      for(k=start ; k<end+1 && p[k]-start<=i-lb ; k++) {
        if(p[k]-offset == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, interval(0.0));
      end++;

      return x[k];
    }

    //! Returns a copy of the i-th (according to the currently used indexing) element of the vector slice.
    /*! If the i-th element is explicitly stored, a copy to the value is returned. If is not explicitly stored,
     *  zero will be returned. This is the const-version of this operator, added for convenience. It is suggested to use the
     * ()-operator for read access. 
     */
    interval operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("sivector_slice::operator[](const int)"));
#endif
      return (*this)(i);
    }

    //! Returns a copy of the i-th (according to the currently used indexing) element of the vector slice.
    /*! If the i-th element is explicitly stored, a copy of this value is returned. Otherwise, 0.0
     *  will be returned. Unlike with the []-operator, the data structure remains unchanged either way.
     *  Thus this operator should always be used for read-only access to the elements of a sparse vector slice.
     */
    const interval operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("srvector_slice::operator()(const int)"));
#endif
      interval r(0.0);

      for(int k=start ; k<end && p[k]-start<=i-lb ; k++) {
        if(p[k]-start == i-lb) 
          r = x[k];
      }

      return r; 
    }

    //! Assigns v to all elements of the vector slice
    sivector_slice& operator=(const real& v) {
      return sl_vs_assign<sivector_slice,real,interval,std::vector<interval>::iterator>(*this,v);
    }

    //! Assigns v to all elements of the vector slice
    sivector_slice& operator=(const interval& v) {
      return sl_vs_assign<sivector_slice,interval,interval,std::vector<interval>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    sivector_slice& operator=(const srvector_slice& v) {
      return slsl_vv_assign<sivector_slice,srvector_slice,interval,std::vector<interval>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    sivector_slice& operator=(const sivector_slice& v) {
      return slsl_vv_assign<sivector_slice,sivector_slice,interval,std::vector<interval>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    sivector_slice& operator=(const srvector& v) {
      return slsp_vv_assign<sivector_slice,srvector,interval,std::vector<interval>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    sivector_slice& operator=(const sivector& v) {
      return slsp_vv_assign<sivector_slice,sivector,interval,std::vector<interval>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    sivector_slice& operator=(const rvector& v) {
      return slf_vv_assign<sivector_slice,rvector,interval,std::vector<interval>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    sivector_slice& operator=(const ivector& v) {
      return slf_vv_assign<sivector_slice,ivector,interval,std::vector<interval>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    sivector_slice& operator=(const rvector_slice& v) {
      return slf_vv_assign<sivector_slice,rvector_slice,interval,std::vector<interval>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    sivector_slice& operator=(const ivector_slice& v) {
      return slf_vv_assign<sivector_slice,ivector_slice,interval,std::vector<interval>::iterator>(*this,v);
    }

    //! Operator for multiplication with a scalar, result is assigned to the vector slice
    sivector_slice& operator*=(const real& s) {
      return sl_vs_multassign(*this,s);
    }

    //! Operator for multiplication with an interval, result is assigned to the vector slice
    sivector_slice& operator*=(const interval& s) {
      return sl_vs_multassign(*this,s);
    }

    //! Operator for division of each element of the vector slice with a scalar, result is assigned to the vector slice
    sivector_slice& operator/=(const real& s) {
      return sl_vs_divassign(*this,s);
    }

    //! Operator for division of each element of the vector slice with an interval, result is assigned to the vector slice
    sivector_slice& operator/=(const interval& s) {
      return sl_vs_divassign(*this,s);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    sivector_slice& operator+=(const rvector& v) {
      return slf_vv_addassign<sivector_slice,rvector,interval>(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    sivector_slice& operator+=(const ivector& v) {
      return slf_vv_addassign<sivector_slice,ivector,interval>(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    sivector_slice& operator+=(const rvector_slice& v) {
      return slf_vv_addassign<sivector_slice,rvector_slice,interval>(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    sivector_slice& operator+=(const ivector_slice& v) {
      return slf_vv_addassign<sivector_slice,ivector_slice,interval>(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    sivector_slice& operator+=(const srvector& v) {
      return slsp_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    sivector_slice& operator+=(const sivector& v) {
      return slsp_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    sivector_slice& operator+=(const srvector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    sivector_slice& operator+=(const sivector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    sivector_slice& operator-=(const rvector& v) {
      return slf_vv_subassign<sivector_slice,rvector,interval>(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    sivector_slice& operator-=(const ivector& v) {
      return slf_vv_subassign<sivector_slice,ivector,interval>(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    sivector_slice& operator-=(const rvector_slice& v) {
      return slf_vv_subassign<sivector_slice,rvector_slice,interval>(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    sivector_slice& operator-=(const ivector_slice& v) {
      return slf_vv_subassign<sivector_slice,ivector_slice,interval>(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    sivector_slice& operator-=(const srvector& v) {
      return slsp_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    sivector_slice& operator-=(const sivector& v) {
      return slsp_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    sivector_slice& operator-=(const srvector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    sivector_slice& operator-=(const sivector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }

    //! Operator for element-wise convex hull with a vector, result is assigned to the vector slice
    sivector_slice& operator|=(const rvector& v) {
      return slf_vv_hullassign<sivector_slice,rvector,interval>(*this,v);
    }

    //! Operator for element-wise convex hull with a vector, result is assigned to the vector slice
    sivector_slice& operator|=(const ivector& v) {
      return slf_vv_hullassign<sivector_slice,ivector,interval>(*this,v);
    }

    //! Operator for element-wise convex hull with a vector, result is assigned to the vector slice
    sivector_slice& operator|=(const rvector_slice& v) {
      return slf_vv_hullassign<sivector_slice,rvector_slice,interval>(*this,v);
    }

    //! Operator for element-wise convex hull with a vector, result is assigned to the vector slice
    sivector_slice& operator|=(const ivector_slice& v) {
      return slf_vv_hullassign<sivector_slice,ivector_slice,interval>(*this,v);
    }

    //! Operator for element-wise convex hull with a vector, result is assigned to the vector slice
    sivector_slice& operator|=(const srvector& v) {
      return slsp_vv_hullassign(*this,v);
    }

    //! Operator for element-wise convex hull with a vector, result is assigned to the vector slice
    sivector_slice& operator|=(const sivector& v) {
      return slsp_vv_hullassign(*this,v);
    }

    //! Operator for element-wise convex hull with a vector, result is assigned to the vector slice
    sivector_slice& operator|=(const srvector_slice& v) {
      return slsl_vv_hullassign(*this,v);
    }

    //! Operator for element-wise convex hull with a vector, result is assigned to the vector slice
    sivector_slice& operator|=(const sivector_slice& v) {
      return slsl_vv_hullassign(*this,v);
    }

    //! Operator for element-wise intersection with a vector, result is assigned to the vector slice
    sivector_slice& operator&=(const ivector& v) {
      return slf_vv_intersectassign<sivector_slice,ivector,interval>(*this,v);
    }

    //! Operator for element-wise intersection with a vector, result is assigned to the vector slice
    sivector_slice& operator&=(const ivector_slice& v) {
      return slf_vv_intersectassign<sivector_slice,ivector_slice,interval>(*this,v);
    }

    //! Operator for element-wise intersection with a vector, result is assigned to the vector slice
    sivector_slice& operator&=(const sivector& v) {
      return slsp_vv_intersectassign(*this,v);
    }

    //! Operator for element-wise intersection with a vector, result is assigned to the vector slice
    sivector_slice& operator&=(const sivector_slice& v) {
      return slsl_vv_intersectassign(*this,v);
    }

    friend int Lb(const sivector_slice&);
    friend int Ub(const sivector_slice&);
    friend srvector Inf(const sivector_slice&);
    friend srvector Sup(const sivector_slice&);
    friend sivector abs(const sivector_slice&);
    friend srvector mid(const sivector_slice&);
    friend srvector diam(const sivector_slice&);
    friend int VecLen(const sivector_slice&);

//     friend srvector operator*(const srmatrix&, const srvector_slice&); //ok
//     friend srvector operator*(const srmatrix_slice&, const srvector_slice&); //ok

    friend class srvector;
    friend class sivector;
    friend class scivector;
    friend class ivector;
    friend class ivector_slice;
    friend class civector;
    friend class civector_slice;

#include "vector_friend_declarations.inl"
};

inline ivector::ivector(const srvector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new interval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline ivector::ivector(const sivector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new interval[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline ivector_slice& ivector_slice::operator=(const srvector& v) {
  *this = rvector(v);
  return *this;
}

inline ivector_slice& ivector_slice::operator=(const srvector_slice& v) {
  *this = rvector(v);
  return *this;
}

inline ivector_slice& ivector_slice::operator=(const sivector& v) {
  *this = ivector(v);
  return *this;
}

inline ivector_slice& ivector_slice::operator=(const sivector_slice& v) {
  *this = ivector(v);
  return *this;
}

inline sivector::sivector(const srvector_slice& s) : lb(s.lb), ub(s.ub), n(s.n)  {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(interval(s.x[i]));
  }

}

inline sivector::sivector(const sivector_slice& s) : lb(s.lb), ub(s.ub), n(s.n) {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(s.x[i]);
  }

}

inline sivector& sivector::operator=(const srvector_slice& v) {
  return spsl_vv_assign<sivector,srvector_slice,interval>(*this,v);
}

inline sivector& sivector::operator=(const sivector_slice& v) {
  return spsl_vv_assign<sivector,sivector_slice,interval>(*this,v);
}

inline sivector_slice sivector::operator()(const int i, const int j) {
#if(CXSC_INDEX_CHECK)
  if(i<lb || j>ub) cxscthrow(ELEMENT_NOT_IN_VEC("sivector::operator()(const int,const int)"));
#endif
  return sivector_slice(*this,i,j);
}

//! Returns the vector -v
inline sivector operator-(const sivector_slice& v) {
  return sl_v_negative<sivector_slice,sivector>(v);
}

//! Returns the lower index bound of the vector slice v
inline int Lb(const sivector_slice& v) {
  return v.lb;
}

//! Returns the upper index bound of the vector slice v
inline int Ub(const sivector_slice& v) {
  return v.ub;
}

//! Returns the infimum vector slice v
inline srvector Inf(const sivector_slice& v) {
  return Inf(sivector(v));
}

//! Returns the supremum of the vector slice v
inline srvector Sup(const sivector_slice& v) {
  return Sup(sivector(v));
}

//! Computes the component-wise absolute values as the interval hull of \f$  \{ |v| \mid v \in [v] \} \f$ for a vector v
inline sivector abs(const sivector_slice& v) {
  sivector res(v.n, v.nnz);
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=v.start ; i<=v.end ; i++)
    res.x.push_back(abs(v.x[i]));
  return res;
}

//! Computes the midpoint vector of v
inline srvector mid(const sivector_slice& v) {
  srvector res(v.n, v.nnz);
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=v.start ; i<=v.end ; i++)
    res.x.push_back(mid(v.x[i]));
  return res;
}

//! Computes the diameter of v
inline srvector diam(const sivector_slice& v) {
  srvector res(v.n, v.nnz);
  res.lb = v.lb;
  res.ub = v.ub;
  res.p = v.p;
  for(int i=v.start ; i<v.end ; i++)
    res.x.push_back(diam(v.x[i]));
  return res;
}

//! Returns the length of the vector slice
inline int VecLen(const sivector_slice& v) {
  return v.n;
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_mult<sivector_slice,rvector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_mult<srvector_slice,ivector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_mult<sivector_slice,ivector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_mult<ivector,srvector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_mult<rvector,sivector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_mult<ivector,sivector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_mult<sivector_slice,rvector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_mult<srvector_slice,ivector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_mult<sivector_slice,ivector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_mult<ivector_slice,srvector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_mult<rvector_slice,sivector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_mult<ivector_slice,sivector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_mult<sivector,srvector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_mult<srvector,sivector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_mult<sivector,sivector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_mult<sivector_slice,srvector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_mult<srvector_slice,sivector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_mult<sivector_slice,sivector,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_mult<sivector_slice,srvector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const srvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_mult<srvector_slice,sivector_slice,interval,sparse_idot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline interval operator*(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_mult<sivector_slice,sivector_slice,interval,sparse_idot>(v1,v2);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline sivector operator*(const sivector_slice& v, const real& s) {
  return sp_vs_mult<sivector_slice,real,sivector>(v,s);
}

//! Multiplies all elements of v with the interval s and returns the result as a new vector
inline sivector operator*(const sivector_slice& v, const interval& s) {
  return sp_vs_mult<sivector_slice,interval,sivector>(v,s);
}

//! Multiplies all elements of v with the interval s and returns the result as a new vector
inline sivector operator*(const srvector_slice& v, const interval& s) {
  return sp_vs_mult<srvector_slice,interval,sivector>(v,s);
}

//! Divides all elements of v with by the scalar s and returns the result as a new vector
inline sivector operator/(const sivector_slice& v, const real& s) {
  return sp_vs_div<sivector_slice,real,sivector>(v,s);
}

//! Divides all elements of v with by the interval s and returns the result as a new vector
inline sivector operator/(const sivector_slice& v, const interval& s) {
  return sp_vs_div<sivector_slice,interval,sivector>(v,s);
}

//! Divides all elements of v with by the interval s and returns the result as a new vector
inline sivector operator/(const srvector_slice& v, const interval& s) {
  return sp_vs_div<srvector_slice,interval,sivector>(v,s);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline sivector operator*(const real& s, const sivector_slice& v) {
  return sp_sv_mult<real,sivector_slice,sivector>(s,v);
}

//! Multiplies all elements of v with the interval s and returns the result as a new vector
inline sivector operator*(const interval& s, const sivector_slice& v) {
  return sp_sv_mult<interval,sivector_slice,sivector>(s,v);
}

//! Multiplies all elements of v with the interval s and returns the result as a new vector
inline sivector operator*(const interval& s, const srvector_slice& v) {
  return sp_sv_mult<interval,srvector_slice,sivector>(s,v);
}

//! Element-wise addition of v1 and v2
inline ivector operator+(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_add<ivector,srvector_slice,ivector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline ivector operator+(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_add<rvector,sivector_slice,ivector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline ivector operator+(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_add<ivector,sivector_slice,ivector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline ivector operator+(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_add<sivector_slice,rvector,ivector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline ivector operator+(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_add<srvector_slice,ivector,ivector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline ivector operator+(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_add<sivector_slice,ivector,ivector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline ivector operator+(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_add<ivector_slice,srvector_slice,ivector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline ivector operator+(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_add<rvector_slice,sivector_slice,ivector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline ivector operator+(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_add<ivector_slice,sivector_slice,ivector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline ivector operator+(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_add<sivector_slice,rvector_slice,ivector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline ivector operator+(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_add<srvector_slice,ivector_slice,ivector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline ivector operator+(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_add<sivector_slice,ivector_slice,ivector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline sivector operator+(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_add<sivector_slice,srvector_slice,sivector,interval>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline sivector operator+(const srvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_add<srvector_slice,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline sivector operator+(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_add<sivector_slice,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline sivector operator+(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_add<sivector,srvector_slice,sivector,interval>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline sivector operator+(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_add<srvector,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline sivector operator+(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_add<sivector,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline sivector operator+(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_add<sivector_slice,srvector,sivector,interval>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline sivector operator+(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_add<srvector_slice,sivector,sivector,interval>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline sivector operator+(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_add<sivector_slice,sivector,sivector,interval>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline ivector operator-(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_sub<ivector,srvector_slice,ivector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline ivector operator-(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_sub<rvector,sivector_slice,ivector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline ivector operator-(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_sub<ivector,sivector_slice,ivector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline ivector operator-(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_sub<sivector_slice,rvector,ivector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline ivector operator-(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_sub<srvector_slice,ivector,ivector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline ivector operator-(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_sub<sivector_slice,ivector,ivector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline ivector operator-(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_sub<ivector_slice,srvector_slice,ivector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline ivector operator-(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_sub<rvector_slice,sivector_slice,ivector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline ivector operator-(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_sub<ivector_slice,sivector_slice,ivector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline ivector operator-(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_sub<sivector_slice,rvector_slice,ivector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline ivector operator-(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_sub<srvector_slice,ivector_slice,ivector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline ivector operator-(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_sub<sivector_slice,ivector_slice,ivector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline sivector operator-(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_sub<sivector_slice,srvector_slice,sivector,interval>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline sivector operator-(const srvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_sub<srvector_slice,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline sivector operator-(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_sub<sivector_slice,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline sivector operator-(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_sub<sivector,srvector_slice,sivector,interval>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline sivector operator-(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_sub<srvector,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline sivector operator-(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_sub<sivector,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline sivector operator-(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_sub<sivector_slice,srvector,sivector,interval>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline sivector operator-(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_sub<srvector_slice,sivector,sivector,interval>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline sivector operator-(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_sub<sivector_slice,sivector,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const rvector& v1, const srvector_slice& v2) {
  return fsl_vv_hull<rvector,srvector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const srvector_slice& v1, const rvector& v2) {
  return slf_vv_hull<srvector_slice,rvector,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const rvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_hull<rvector_slice,srvector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const srvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_hull<srvector_slice,rvector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline sivector operator|(const srvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_hull<srvector_slice,srvector_slice,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline sivector operator|(const srvector& v1, const srvector_slice& v2) {
  return spsl_vv_hull<srvector,srvector_slice,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline sivector operator|(const srvector_slice& v1, const srvector& v2) {
  return slsp_vv_hull<srvector_slice,srvector,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_hull<ivector,srvector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_hull<rvector,sivector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_hull<ivector,sivector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_hull<sivector_slice,rvector,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_hull<srvector_slice,ivector,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_hull<sivector_slice,ivector,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_hull<ivector_slice,srvector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_hull<rvector_slice,sivector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_hull<ivector_slice,sivector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_hull<sivector_slice,rvector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_hull<srvector_slice,ivector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline ivector operator|(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_hull<sivector_slice,ivector_slice,ivector>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline sivector operator|(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_hull<sivector_slice,srvector_slice,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline sivector operator|(const srvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_hull<srvector_slice,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline sivector operator|(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_hull<sivector_slice,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline sivector operator|(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_hull<sivector,srvector_slice,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline sivector operator|(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_hull<srvector,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline sivector operator|(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_hull<sivector,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline sivector operator|(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_hull<sivector_slice,srvector,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline sivector operator|(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_hull<srvector_slice,sivector,sivector,interval>(v1,v2);
}

//! Element-wise convex hull of v1 and v2
inline sivector operator|(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_hull<sivector_slice,sivector,sivector,interval>(v1,v2);
}

//! Element-wise intersection of v1 and v2
inline ivector operator&(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_intersect<ivector,sivector_slice,ivector>(v1,v2);
}

//! Element-wise intersection of v1 and v2
inline ivector operator&(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_intersect<sivector_slice,ivector,ivector>(v1,v2);
}

//! Element-wise intersection of v1 and v2
inline ivector operator&(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_intersect<ivector_slice,sivector_slice,ivector>(v1,v2);
}

//! Element-wise intersection of v1 and v2
inline ivector operator&(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_intersect<sivector_slice,ivector_slice,ivector>(v1,v2);
}

//! Element-wise intersection of v1 and v2
inline sivector operator&(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_intersect<sivector_slice,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise intersection of v1 and v2
inline sivector operator&(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_intersect<sivector,sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise intersection of v1 and v2
inline sivector operator&(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_intersect<sivector_slice,sivector,sivector,interval>(v1,v2);
}

inline ivector& ivector::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline ivector& ivector::operator+=(const sivector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator+=(const sivector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline sivector& sivector::operator+=(const srvector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline sivector& sivector::operator+=(const sivector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline ivector& ivector::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline ivector& ivector::operator-=(const sivector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator-=(const sivector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline sivector& sivector::operator-=(const srvector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

inline sivector& sivector::operator-=(const sivector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

inline ivector& ivector::operator|=(const srvector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline ivector& ivector::operator|=(const sivector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator|=(const srvector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator|=(const sivector_slice& v2) {
  return fsl_vv_hullassign(*this,v2);
}

inline ivector& ivector::operator&=(const sivector_slice& v2) {
  return fsl_vv_intersectassign(*this,v2);
}

inline ivector_slice& ivector_slice::operator&=(const sivector_slice& v2) {
  return fsl_vv_intersectassign(*this,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector_slice& v1, const srvector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector_slice& v1, const sivector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector_slice& v1, const sivector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector_slice& v1, const rvector& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector_slice& v1, const ivector& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector_slice& v1, const ivector& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const ivector& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const rvector& v1, const sivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const ivector& v1, const sivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector_slice& v1, const srvector& v2) {
  return !slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector_slice& v1, const sivector& v2) {
  return !slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector_slice& v1, const sivector& v2) {
  return !slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector& v1, const srvector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector& v1, const sivector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector& v1, const sivector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector_slice& v1, const rvector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector_slice& v1, const ivector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const sivector_slice& v1, const ivector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const ivector_slice& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const rvector_slice& v1, const sivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const ivector_slice& v1, const sivector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_less<srvector_slice,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_less<sivector_slice,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_less<srvector_slice,sivector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_less<sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_less<srvector,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_less<sivector,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_less<srvector_slice,ivector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_less<sivector_slice,ivector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_less<rvector,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_less<ivector,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_less<srvector_slice,ivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_less<sivector_slice,ivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_less<rvector_slice,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than the respective elements of v2.
 */
inline bool operator<(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_less<ivector_slice,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_leq<sivector_slice,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector_slice& v1, const sivector& v2) {
  return slsp_vv_leq<srvector_slice,sivector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_leq<sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector& v1, const sivector_slice& v2) {
  return spsl_vv_leq<srvector,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_leq<sivector,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector_slice& v1, const ivector& v2) {
  return slf_vv_leq<srvector_slice,ivector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_leq<sivector_slice,ivector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const rvector& v1, const sivector_slice& v2) {
  return fsl_vv_leq<rvector,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_leq<ivector,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const srvector_slice& v1, const ivector_slice& v2) {
  return slf_vv_leq<srvector_slice,ivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_leq<sivector_slice,ivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const rvector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_leq<rvector_slice,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are smaller than or equal to the respective elements of v2.
 */
inline bool operator<=(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_leq<ivector_slice,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_greater<sivector_slice,srvector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_greater<sivector_slice,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_greater<sivector_slice,srvector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_greater<sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_greater<sivector,srvector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_greater<sivector,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_greater<sivector_slice,rvector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_greater<sivector_slice,ivector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_greater<ivector,srvector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_greater<ivector,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_greater<sivector_slice,rvector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_greater<sivector_slice,ivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_greater<ivector_slice,srvector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than the respective elements of v2.
 */
inline bool operator>(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_greater<ivector_slice,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_geq<sivector_slice,srvector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector_slice& v1, const sivector_slice& v2) {
  return slsl_vv_geq<sivector_slice,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector_slice& v1, const srvector& v2) {
  return slsp_vv_geq<sivector_slice,srvector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector_slice& v1, const sivector& v2) {
  return slsp_vv_geq<sivector_slice,sivector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector& v1, const srvector_slice& v2) {
  return spsl_vv_geq<sivector,srvector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector& v1, const sivector_slice& v2) {
  return spsl_vv_geq<sivector,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector_slice& v1, const rvector& v2) {
  return slf_vv_geq<sivector_slice,rvector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector_slice& v1, const ivector& v2) {
  return slf_vv_geq<sivector_slice,ivector,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const ivector& v1, const srvector_slice& v2) {
  return fsl_vv_geq<ivector,srvector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const ivector& v1, const sivector_slice& v2) {
  return fsl_vv_geq<ivector,sivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector_slice& v1, const rvector_slice& v2) {
  return slf_vv_geq<sivector_slice,rvector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const sivector_slice& v1, const ivector_slice& v2) {
  return slf_vv_geq<sivector_slice,ivector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const ivector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_geq<ivector_slice,srvector_slice,interval>(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are larger than or equal to the respective elements of v2.
 */
inline bool operator>=(const ivector_slice& v1, const sivector_slice& v2) {
  return fsl_vv_geq<ivector_slice,sivector_slice,interval>(v1,v2);
}

//! Output operator for sparse vector slice v.
/**
 *  The output format is set by global flags, default is dense output.
 *  Use cout << SparseInOut; for sparse output or cout << MatrixMarketInOut; for
 *  output in matrix market format.
 */
inline std::ostream& operator<<(std::ostream& os, const sivector_slice& v) {
  return sl_v_output<sivector_slice,interval>(os,v);
}

//! Input operator for sparse vector slice v.
/**
 *  The input format is set by global flags, default is dense input.
 *  Use cout << SparseInOut; for sparse input or cout << MatrixMarketInOut; for
 *  input in matrix market format.
 */
inline std::istream& operator>>(std::istream& is, sivector_slice& v) {
  return sl_v_input<sivector_slice,interval>(is,v);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector& x, const sivector& y) {
  spsp_vv_accu<idotprecision,sivector,sivector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector& x, const srvector& y) {
  spsp_vv_accu<idotprecision,sivector,srvector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector& x, const sivector& y) {
  spsp_vv_accu<idotprecision,srvector,sivector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector& x, const ivector& y) {
  spf_vv_accu<idotprecision,sivector,ivector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector& x, const rvector& y) {
  spf_vv_accu<idotprecision,sivector,rvector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector& x, const ivector& y) {
  spf_vv_accu<idotprecision,srvector,ivector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector& x, const ivector_slice& y) {
  spf_vv_accu<idotprecision,sivector,ivector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector& x, const rvector_slice& y) {
  spf_vv_accu<idotprecision,sivector,rvector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector& x, const ivector_slice& y) {
  spf_vv_accu<idotprecision,srvector,ivector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const ivector& x, const sivector& y) {
  fsp_vv_accu<idotprecision,ivector,sivector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const ivector& x, const srvector& y) {
  fsp_vv_accu<idotprecision,ivector,srvector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const rvector& x, const sivector& y) {
  fsp_vv_accu<idotprecision,rvector,sivector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const ivector_slice& x, const sivector& y) {
  fsp_vv_accu<idotprecision,ivector_slice,sivector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const ivector_slice& x, const srvector& y) {
  fsp_vv_accu<idotprecision,ivector_slice,srvector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const rvector_slice& x, const sivector& y) {
  fsp_vv_accu<idotprecision,rvector_slice,sivector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector_slice& x, const ivector& y) {
  slf_vv_accu<idotprecision,sivector_slice,ivector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector_slice& x, const rvector& y) {
  slf_vv_accu<idotprecision,sivector_slice,rvector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector_slice& x, const ivector& y) {
  slf_vv_accu<idotprecision,srvector_slice,ivector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector_slice& x, const ivector_slice& y) {
  slf_vv_accu<idotprecision,sivector_slice,ivector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector_slice& x, const rvector_slice& y) {
  slf_vv_accu<idotprecision,sivector_slice,rvector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector_slice& x, const ivector_slice& y) {
  slf_vv_accu<idotprecision,srvector_slice,ivector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const ivector& x, const sivector_slice& y) {
  fsl_vv_accu<idotprecision,ivector,sivector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const ivector& x, const srvector_slice& y) {
  fsl_vv_accu<idotprecision,ivector,srvector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const rvector& x, const sivector_slice& y) {
  fsl_vv_accu<idotprecision,rvector,sivector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const ivector_slice& x, const sivector_slice& y) {
  fsl_vv_accu<idotprecision,ivector_slice,sivector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const ivector_slice& x, const srvector_slice& y) {
  fsl_vv_accu<idotprecision,ivector_slice,srvector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const rvector_slice& x, const sivector_slice& y) {
  fsl_vv_accu<idotprecision,rvector_slice,sivector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector_slice& x, const sivector_slice& y) {
  slsl_vv_accu<idotprecision,sivector_slice,sivector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector_slice& x, const srvector_slice& y) {
  slsl_vv_accu<idotprecision,sivector_slice,srvector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector_slice& x, const sivector_slice& y) {
  slsl_vv_accu<idotprecision,srvector_slice,sivector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector& x, const sivector_slice& y) {
  spsl_vv_accu<idotprecision,sivector,sivector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector& x, const srvector_slice& y) {
  spsl_vv_accu<idotprecision,sivector,srvector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector& x, const sivector_slice& y) {
  spsl_vv_accu<idotprecision,srvector,sivector_slice,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector_slice& x, const sivector& y) {
  slsp_vv_accu<idotprecision,sivector_slice,sivector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const sivector_slice& x, const srvector& y) {
  slsp_vv_accu<idotprecision,sivector_slice,srvector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(idotprecision& dot, const srvector_slice& x, const sivector& y) {
  slsp_vv_accu<idotprecision,srvector_slice,sivector,sparse_idot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector& x, const srvector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector& x, const ivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector& x, const rvector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector& x, const ivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector& x, const ivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector& x, const rvector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector& x, const ivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const ivector& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const ivector& x, const srvector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const ivector_slice& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const ivector_slice& x, const srvector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector_slice& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector_slice& x, const ivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector_slice& x, const rvector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector_slice& x, const ivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector_slice& x, const ivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector_slice& x, const rvector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector_slice& x, const ivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const ivector& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const ivector& x, const srvector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const ivector_slice& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const ivector_slice& x, const srvector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector_slice& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector_slice& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector_slice& x, const srvector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector_slice& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector& x, const srvector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector& x, const sivector_slice& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector_slice& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const sivector_slice& x, const srvector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector_slice& x, const sivector& y) {
  idotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  SetRe(dot, Re(dot) + tmp);
}

} //namespace cxsc

#include "sparsevector.inl"

#endif 
