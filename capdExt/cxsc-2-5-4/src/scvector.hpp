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

/* CVS $Id: scvector.hpp,v 1.17 2014/01/30 17:23:48 cxsc Exp $ */

#ifndef _CXSC_SCVECTOR_HPP_INCLUDED
#define _CXSC_SCVECTOR_HPP_INCLUDED

#include <complex.hpp>
#include <cvector.hpp>
#include <vector>
#include <iostream>
#include <srvector.hpp>
#include <sparsecdot.hpp>
#include <sparsevector.hpp>

namespace cxsc {

class srvector_slice;
class srmatrix;
class srmatrix_slice;
class srmatrix_subv;
class scvector_slice;
class scmatrix;
class scmatrix_slice;
class scmatrix_subv;
class scivector;
class scivector_slice;
class scimatrix;
class scimatrix_slice;
class scimatrix_subv;

//! A sparse complex vector
/*!
This data type represents a sparse complex vector. Only the non zero elements are stored explicitly with their value and
the respective index. All operators are overloaded to take advantage of the sparsity.
*/
class scvector {
  private:
    std::vector<int> p;
    std::vector<complex> x;
    int lb;
    int ub;
    int n; 

  public:
    //! Default constructor, creates an empty vector of size 0    
    scvector() : lb(0), ub(-1) , n(0) {
    }

    //! Constructor for creating an empty vector of size s
    explicit scvector(const int s) : lb(1), ub(s), n(s) {
	p.reserve((int)(s*0.1));
	x.reserve((int)(s*0.1));
    }

    //! Constructor for creating an empty vector of size s and reserving memory for b elements
    scvector(const int s, const int b) : lb(1), ub(s), n(s) {
	p.reserve(b);
	x.reserve(b);
    }

    //! Constructor for creating a sparse vector our of a dense vector v. Only the non-zero elements of v are stored explicitly.
    scvector(const cvector& v) : lb(Lb(v)), ub(Ub(v)), n(VecLen(v)) {
        for(int i=lb ; i<=ub ; i++) {
          if(v[i] != 0.0) {
            p.push_back(i-lb);
            x.push_back(v[i]);
          }
        }
    }

    //! Constructor for creating a sparse vector our of a dense vector v. Only the non-zero elements of v are stored explicitly.
    scvector(const rvector& v) : lb(Lb(v)), ub(Ub(v)), n(VecLen(v)) {
        for(int i=lb ; i<=ub ; i++) {
          if(v[i] != 0.0) {
            p.push_back(i-lb);
            x.push_back(complex(v[i]));
          }
        }
    }

    //! Creates a sparse vector of dimension n with nnz non zeros, who are defined by the elements of index and values
    scvector(const int n, const int nnz, const intvector& index, const cvector& values) : lb(1), ub(n) {
      this->n = n;
      for(int i=0 ; i<nnz ; i++) {
        if(values[i+Lb(values)] != 0.0) {
          p.push_back(index[i+Lb(index)]);
          x.push_back(values[i+Lb(values)]);
        }
      }
    }

    //! Creates a sparse vector of dimension n with nnz non zeros, who are defined by the elements of index and values
    scvector(const int n, const int nnz, const int* index, const complex* values) : lb(1), ub(n) {
      this->n = n;
      for(int i=0 ; i<nnz ; i++) {
        if(values[i] != 0.0) {
          p.push_back(index[i]);
          x.push_back(values[i]);
        }
      } 
    }

    //! Constructor for creating a sparse complex vector our of a sparse real vector
    scvector(const srvector& v) : p(v.p), lb(v.lb), ub(v.ub), n(v.n) {
      x.reserve(v.get_nnz());
      for(int i=0 ; i<v.get_nnz() ; i++) 
        x.push_back(complex(v.x[i]));
    }

    //! Creates a sparse vector out of a sparse vector slice
    scvector(const srvector_slice&);
    //! Creates a sparse vector out of a sparse vector slice
    scvector(const scvector_slice&);
    //! Creates a sparse vector out of a row or column of a sparse matrix
    scvector(const srmatrix_subv& A);
    //! Creates a sparse vector out of a row or column of a sparse matrix
    scvector(const scmatrix_subv& A);

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
     *  refers to the i-th element of the STL-vector storing the values of the elements). */    std::vector<complex>& values() {
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
     *  refers to the i-th element of the STL-vector storing the values of the elements). */    const std::vector<complex>& values() const {
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

    //! Assigns a sparse real vector to a sparse complex vector
    scvector& operator=(const srvector& v) {
      n = v.n;
      p = v.p;
      x.clear();
      x.reserve(v.get_nnz());
      for(unsigned int i=0 ; i<v.x.size() ; i++)
        x[i] = complex(v.x[i]);
      return *this;
    } 

    //! Assigns v to all elements of the vector (resulting in a dense vector!)
    scvector& operator=(const real& v) {
      return sp_vs_assign<scvector,real,complex>(*this,v);
    }

    //! Assigns v to all elements of the vector (resulting in a dense vector!)
    scvector& operator=(const complex& v) {
      return sp_vs_assign<scvector,complex,complex>(*this,v);
    }

    //! Assign the dense vector v to a sparse vector. Only the non zero elements of v are used.
    scvector& operator=(const rvector& v) {
      return spf_vv_assign<scvector,rvector,complex>(*this,v);
    }

    //! Assign the dense vector v to a sparse vector. Only the non zero elements of v are used.
    scvector& operator=(const cvector& v) {
      return spf_vv_assign<scvector,cvector,complex>(*this,v);
    }

    //! Assign the dense vector slice v to a sparse vector. Only the non zero elements of v are used.
    scvector& operator=(const rvector_slice& v) {
      return spf_vv_assign<scvector,rvector_slice,complex>(*this,v);
    }

    //! Assign the dense vector slice v to a sparse vector. Only the non zero elements of v are used.
    scvector& operator=(const cvector_slice& v) {
      return spf_vv_assign<scvector,cvector_slice,complex>(*this,v);
    }

    //! Assign the sparse vector slice v to a sparse vector.
    scvector& operator=(const scvector_slice&);
    //! Assign the sparse vector slice v to a sparse vector.   
    scvector& operator=(const srvector_slice&);

    //! Returns a reference to the i-th (according to the currently used indexing) element of the vector.
    /*! If the i-th element is explicitly stored, a reference to the value is returned. If is not explicitly stored,
     *  it will be added to the data structure as a zero element. The returned reference then points to this new 
     *  element. Hence ths []-operator should only be used for write access to the elements of a sparse vector. 
     *  Use the ()-operator for read access. */    
    complex& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector::operator[](const int)"));
#endif
      int k;

      for(k=0 ; k<get_nnz() && p[k]<=i-lb ; k++) {
        if(p[k] == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, complex(0.0));

      return x[k];
    }

    //! Returns a copy of the i-th (according to the currently used indexing) element of the vector.
    /*! If the i-th element is explicitly stored, a copy to the value is returned. If is not explicitly stored,
     *  zero will be returned. This is the const-Version of this operator, added for convenience. It is suggested to use thei
     * ()-operator for read access. 
     */
    complex operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector::operator[](const int)"));
#endif
      return (*this)(i);
    }

    //! Returns a copy of the i-th (according to the currently used indexing) element of the vector.
    /*! If the i-th element is explicitly stored, a copy of this value is returned. Otherwise, 0.0
     *  will be returned. Unlike with the []-operator, the data structure remains unchanged either way.
     *  Thus this operator should always be used for read-only access to the elements of a sparse vector.
     */
    const complex operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector::operator()(const int)"));
#endif
      complex r(0.0);

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
    scvector operator()(const intvector& per) {
      scvector v(n,get_nnz());
      intvector pinv = perminv(per);

      std::map<int,complex> work;
      for(int i=0 ; i<get_nnz() ; i++)
         work.insert(std::make_pair(pinv[Lb(pinv)+p[i]], x[i]));
 
      for(std::map<int,complex>::iterator it=work.begin() ; it!=work.end() ; it++) {
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
    scvector operator()(const intmatrix& P) {
      intvector p = permvec(P);
      return (*this)(p);
    }

    //! Returns a slice of the vector from the i-th to the j-th (according to the currently used indexing) element.
    /*! This operator can be used for read and write access to a slice of the sparse vector.
     */
    scvector_slice operator()(const int, const int);

    //! Operator for multiplication with a scalar, result is assigned to the vector    
    scvector& operator*=(const real& s) {
      return sp_vs_multassign(*this,s);
    }

    //! Operator for multiplication with a scalar, result is assigned to the vector
    scvector& operator*=(const complex& s) {
      return sp_vs_multassign(*this,s);
    }

    //! Operator for division of each element of the vector with a scalar, result is assigned to the vector
    scvector& operator/=(const real& s) {
      return sp_vs_divassign(*this,s);
    }

    //! Operator for division of each element of the vector with a scalar, result is assigned to the vector
    scvector& operator/=(const complex& s) {
      return sp_vs_divassign(*this,s);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    scvector& operator+=(const rvector& v) {
      return spf_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    scvector& operator+=(const cvector& v) {
      return spf_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    scvector& operator+=(const rvector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    scvector& operator+=(const cvector_slice& v) {
      return spf_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    scvector& operator+=(const srvector& v) {
      return spsp_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    scvector& operator+=(const scvector& v) {
      return spsp_vv_addassign(*this,v);
    }


    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    scvector& operator-=(const rvector& v) {
      return spf_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    scvector& operator-=(const cvector& v) {
      return spf_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    scvector& operator-=(const rvector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    scvector& operator-=(const cvector_slice& v) {
      return spf_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    scvector& operator-=(const srvector& v) {
      return spsp_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    scvector& operator-=(const scvector& v) {
      return spsp_vv_subassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector
    scvector& operator+=(const srvector_slice&);
    //! Operator for element-wise addition with a vector, result is assigned to the vector
    scvector& operator+=(const scvector_slice&);
    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    scvector& operator-=(const srvector_slice&);
    //! Operator for element-wise subtraction with a vector, result is assigned to the vector
    scvector& operator-=(const scvector_slice&);

    friend void SetLb(scvector&, const int);
    friend void SetUb(scvector&, const int);
    friend int Lb(const scvector&);
    friend int Ub(const scvector&);
    friend srvector Re(const scvector&);
    friend srvector Im (const scvector&);
    friend scvector Inf(const scivector&);
    friend scvector Sup (const scivector&);
    friend scvector mid(const scivector&);
    friend scvector diam(const scivector&);
    friend scvector mid(const scivector_slice&);
    friend scvector diam(const scivector_slice&);
    friend int VecLen(const scvector&);
    friend srvector abs(const scvector&);

    friend class srvector_slice;
    friend class scvector_slice;
    friend class scivector_slice;
    friend class scivector;
    friend class cvector;
    friend class cvector_slice;
    friend class civector;
    friend class civector_slice;

#include "vector_friend_declarations.inl"
};

inline cvector::cvector(const scvector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new complex[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline cvector::cvector(const srvector& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new complex[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=0 ; i<v.get_nnz() ; i++)
    dat[v.p[i]] = v.x[i];
}

inline cvector& cvector::operator=(const scvector& v) {
  return fsp_vv_assign<cvector,scvector,complex>(*this,v);
}

inline cvector& cvector::operator=(const scvector_slice& v) {
  return fsl_vv_assign<cvector,scvector_slice,complex>(*this,v);
}

inline cvector& cvector::operator=(const srvector& v) {
  return fsp_vv_assign<cvector,srvector,complex>(*this,v);
}

inline cvector& cvector::operator=(const srvector_slice& v) {
  return fsl_vv_assign<cvector,srvector_slice,complex>(*this,v);
}

//! Sets the lower index bound of the vector v to i
/** 
 * After setting the lower index bound to i, the indexing of the vector is i-based.
 */
inline void SetLb(scvector& v, const int i) {
  v.lb = i;
  v.ub = v.lb + v.n - 1;
}

//! Sets the upper index bound of the vector v to i
/** 
 * After setting the upper index bound to i, the indexing of the vector of dimension n is (i-n+1)-based.
 */
inline void SetUb(scvector& v, const int j) {
  v.ub = j;
  v.lb = v.ub - v.n + 1;
}

//! Returns the lower index bound of the vector v
inline int Lb(const scvector& v) {
  return v.lb;
}

//! Returns the upper index bound of the vector v
inline int Ub(const scvector& v) {
  return v.ub;
}

//! Returns the real part of the complex vector v
inline srvector Re(const scvector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = Re(v.x[i]);
  return res;
}

//! Returns the imaginary part of the complex vector v
inline srvector Im(const scvector& v) {
  srvector res(v.n, v.get_nnz());
  res.lb = v.lb;
  res.ub = v.ub;
  res.p  = v.p;
  for(int i=0 ; i<v.get_nnz() ; i++)
    res.x[i] = Im(v.x[i]);
  return res;
}

//! Returns the vector of component-wise absolute values of v
inline srvector abs(const scvector& v) {
  srvector ret(VecLen(v));
  const std::vector<int>& rv = v.row_indices();
  const std::vector<complex>& xv = v.values();
  std::vector<int>& r = ret.row_indices();  
  std::vector<real>& x = ret.values();
  
  for(unsigned int i=0 ; i<xv.size() ; i++) {
    x.push_back(abs(xv[i]));
    r.push_back(rv[i]);
  }
  
  return ret;
}

//! Returns the length of the vector (the dimension)
inline int VecLen(const scvector& v) {
  return v.n;
}

//! Resizes the vector to length 0 (all elements are deleted)
inline void Resize(scvector& v) {
  sp_v_resize(v);
}

//! Resizes the vector to length n.
/**
 * All elements of the vector that can still be stored after the resizing are copied into the resized vector.
 */
inline void Resize(scvector& v, const int n) {
  sp_v_resize(v,n);
}

//! Resizes the vector to length u-l+1.
/**
 * The new vector has lower index bound l and upper index bound u. 
 * All elements of the vector that can still be stored after the resizing are copied into the resized vector.
 */
inline void Resize(scvector& v, const int l, const int u) {
  sp_v_resize(v,l,u);
}

//! Unary operator, returns -v
inline scvector operator-(const scvector& v) {
  return sp_v_negative(v);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector& v1, const cvector& v2) {
  return spf_vv_mult<scvector,cvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector& v1, const rvector& v2) {
  return spf_vv_mult<scvector,rvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const srvector& v1, const cvector& v2) {
  return spf_vv_mult<srvector,cvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const rvector& v1, const scvector& v2) {
  return fsp_vv_mult<rvector,scvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const cvector& v1, const srvector& v2) {
  return fsp_vv_mult<cvector,srvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const cvector& v1, const scvector& v2) {
  return fsp_vv_mult<cvector,scvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector& v1, const rvector_slice& v2) {
  return spf_vv_mult<scvector,rvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector& v1, const cvector_slice& v2) {
  return spf_vv_mult<scvector,cvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const srvector& v1, const cvector_slice& v2) {
  return spf_vv_mult<srvector,cvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const cvector_slice& v1, const srvector& v2) {
  return fsp_vv_mult<cvector_slice,srvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const cvector_slice& v1, const scvector& v2) {
  return fsp_vv_mult<cvector_slice,scvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const rvector_slice& v1, const scvector& v2) {
  return fsp_vv_mult<rvector_slice,scvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector& v1, const srvector& v2) {
  return spsp_vv_mult<scvector,srvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const srvector& v1, const scvector& v2) {
  return spsp_vv_mult<srvector,scvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector& v1, const scvector& v2) {
  return spsp_vv_mult<scvector,scvector,complex,sparse_cdot>(v1,v2);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline scvector operator*(const scvector& v, const real& s) {
  return sp_vs_mult<scvector,real,scvector>(v,s);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline scvector operator*(const scvector& v, const complex& s) {
  return sp_vs_mult<scvector,complex,scvector>(v,s);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline scvector operator*(const srvector& v, const complex& s) {
  return sp_vs_mult<srvector,complex,scvector>(v,s);
}

//! Divides all elements of v by the scalar s and returns the result as a new vector
inline scvector operator/(const scvector& v, const real& s) {
  return sp_vs_div<scvector,real,scvector>(v,s);
}

//! Divides all elements of v by the scalar s and returns the result as a new vector
inline scvector operator/(const scvector& v, const complex& s) {
  return sp_vs_div<scvector,complex,scvector>(v,s);
}

//! Divides all elements of v by the scalar s and returns the result as a new vector
inline scvector operator/(const srvector& v, const complex& s) {
  return sp_vs_div<srvector,complex,scvector>(v,s);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline scvector operator*(const real& s, const scvector& v) {
  return sp_sv_mult<real,scvector,scvector>(s,v);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline scvector operator*(const complex& s, const scvector& v) {
  return sp_sv_mult<complex,scvector,scvector>(s,v);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline scvector operator*(const complex& s, const srvector& v) {
  return sp_sv_mult<complex,srvector,scvector>(s,v);
}

//! Element-wise addition of the vectors v1 and v2
inline cvector operator+(const cvector& v1, const srvector& v2) {
  return fsp_vv_add<cvector,srvector,cvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline cvector operator+(const rvector& v1, const scvector& v2) {
  return fsp_vv_add<rvector,scvector,cvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline cvector operator+(const cvector& v1, const scvector& v2) {
  return fsp_vv_add<cvector,scvector,cvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline cvector operator+(const scvector& v1, const rvector& v2) {
  return spf_vv_add<scvector,rvector,cvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline cvector operator+(const srvector& v1, const cvector& v2) {
  return spf_vv_add<srvector,cvector,cvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline cvector operator+(const scvector& v1, const cvector& v2) {
  return spf_vv_add<scvector,cvector,cvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline cvector operator+(const cvector_slice& v1, const srvector& v2) {
  return fsp_vv_add<cvector_slice,srvector,cvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline cvector operator+(const rvector_slice& v1, const scvector& v2) {
  return fsp_vv_add<rvector_slice,scvector,cvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline cvector operator+(const cvector_slice& v1, const scvector& v2) {
  return fsp_vv_add<cvector_slice,scvector,cvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline cvector operator+(const scvector& v1, const rvector_slice& v2) {
  return spf_vv_add<scvector,rvector_slice,cvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline cvector operator+(const srvector& v1, const cvector_slice& v2) {
  return spf_vv_add<srvector,cvector_slice,cvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline cvector operator+(const scvector& v1, const cvector_slice& v2) {
  return spf_vv_add<scvector,cvector_slice,cvector>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline scvector operator+(const scvector& v1, const srvector& v2) {
  return spsp_vv_add<scvector,srvector,scvector,complex>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline scvector operator+(const srvector& v1, const scvector& v2) {
  return spsp_vv_add<srvector,scvector,scvector,complex>(v1,v2);
}

//! Element-wise addition of the vectors v1 and v2
inline scvector operator+(const scvector& v1, const scvector& v2) {
  return spsp_vv_add<scvector,scvector,scvector,complex>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline cvector operator-(const cvector& v1, const srvector& v2) {
  return fsp_vv_sub<cvector,srvector,cvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline cvector operator-(const rvector& v1, const scvector& v2) {
  return fsp_vv_sub<rvector,scvector,cvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline cvector operator-(const cvector& v1, const scvector& v2) {
  return fsp_vv_sub<cvector,scvector,cvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline cvector operator-(const scvector& v1, const rvector& v2) {
  return spf_vv_sub<scvector,rvector,cvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline cvector operator-(const srvector& v1, const cvector& v2) {
  return spf_vv_sub<srvector,cvector,cvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline cvector operator-(const scvector& v1, const cvector& v2) {
  return spf_vv_sub<scvector,cvector,cvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline cvector operator-(const cvector_slice& v1, const srvector& v2) {
  return fsp_vv_sub<cvector_slice,srvector,cvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline cvector operator-(const rvector_slice& v1, const scvector& v2) {
  return fsp_vv_sub<rvector_slice,scvector,cvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline cvector operator-(const cvector_slice& v1, const scvector& v2) {
  return fsp_vv_sub<cvector_slice,scvector,cvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline cvector operator-(const scvector& v1, const rvector_slice& v2) {
  return spf_vv_sub<scvector,rvector_slice,cvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline cvector operator-(const srvector& v1, const cvector_slice& v2) {
  return spf_vv_sub<srvector,cvector_slice,cvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline cvector operator-(const scvector& v1, const cvector_slice& v2) {
  return spf_vv_sub<scvector,cvector_slice,cvector>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline scvector operator-(const scvector& v1, const srvector& v2) {
  return spsp_vv_sub<scvector,srvector,scvector,complex>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline scvector operator-(const srvector& v1, const scvector& v2) {
  return spsp_vv_sub<srvector,scvector,scvector,complex>(v1,v2);
}

//! Element-wise subtraction of the vectors v1 and v2
inline scvector operator-(const scvector& v1, const scvector& v2) {
  return spsp_vv_sub<scvector,scvector,scvector,complex>(v1,v2);
}

inline cvector& cvector::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline cvector& cvector::operator+=(const scvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator+=(const srvector& v2) {
  return fsp_vv_addassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator+=(const scvector& v2) {
  return fsp_vv_addassign(*this,v2);
}
 
inline cvector& cvector::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline cvector& cvector::operator-=(const scvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator-=(const srvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator-=(const scvector& v2) {
  return fsp_vv_subassign(*this,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const scvector& v1, const scvector& v2) {
  return spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const scvector& v1, const srvector& v2) {
  return spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const srvector& v1, const scvector& v2) {
  return spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const scvector& v1, const rvector& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const srvector& v1, const cvector& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const scvector& v1, const cvector& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const cvector& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const rvector& v1, const scvector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const cvector& v1, const scvector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const scvector& v1, const rvector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const srvector& v1, const cvector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const scvector& v1, const cvector_slice& v2) {
  return spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const cvector_slice& v1, const srvector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const rvector_slice& v1, const scvector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are euqal to the respective elements of v2.
 */
inline bool operator==(const cvector_slice& v1, const scvector& v2) {
  return fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector& v1, const srvector& v2) {
  return !spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector& v1, const scvector& v2) {
  return !spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector& v1, const scvector& v2) {
  return !spsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector& v1, const rvector& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector& v1, const cvector& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector& v1, const cvector& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const cvector& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const rvector& v1, const scvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const cvector& v1, const scvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector& v1, const rvector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector& v1, const cvector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector& v1, const cvector_slice& v2) {
  return !spf_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const cvector_slice& v1, const srvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const rvector_slice& v1, const scvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Element-wise comparison of the vectors v1 and v2. 
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const cvector_slice& v1, const scvector& v2) {
  return !fsp_vv_comp(v1,v2);
}

//! Output operator for sparse vector v.
/**
 *  The output format is set by global flags, default is dense output.
 *  Use cout << SparseInOut; for sparse output or cout << MatrixMarketInOut; for
 *  output in matrix market format.
 */
inline std::ostream& operator<<(std::ostream& os, const scvector& v) {
  return sp_v_output<scvector,complex>(os,v);
}

//! Input operator for sparse vector v.
/**
 *  The input format is set by global flags, default is dense input.
 *  Use cout << SparseInOut; for sparse input or cout << MatrixMarketInOut; for
 *  input in matrix market format.
 */
inline std::istream& operator>>(std::istream& is, scvector& v) {
  return sp_v_input<scvector,complex>(is,v);
}


//! Helper class for slices of sparse vectors.
/**
 * This class stores a reference to a sparse vector and operates on a slice of it. This class
 * is used internally by C-XSC, it should normally not be necessary for the user to use it explicitly.
 * 
 */
class scvector_slice {
  private:
    std::vector<int>& p;
    std::vector<complex>& x;
    scvector& orig;
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
    scvector_slice(scvector& v, int l, int u) : p(v.p), x(v.x), orig(v), lb(l), ub(u), n(u-l+1)  {
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
    complex& operator[](const int i) {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector_slice::operator[](const int)"));
#endif
      int k;

      for(k=start ; k<end+1 && p[k]-start<=i-lb ; k++) {
        if(p[k]-offset == i-lb) 
          return x[k];
      }

      p.insert(p.begin() + k, i-lb);
      x.insert(x.begin() + k, complex(0.0));
      end++;

      return x[k];
    }

    //! Returns a copy of the i-th (according to the currently used indexing) element of the vector slice.
    /*! If the i-th element is explicitly stored, a copy to the value is returned. If is not explicitly stored,
     *  zero will be returned. This is the const-version of this operator, added for convenience. It is suggested to use the
     * ()-operator for read access. 
     */
    complex operator[](const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector_slice::operator[](const int)"));
#endif
      return (*this)(i);
    }

    //! Returns a copy of the i-th (according to the currently used indexing) element of the vector slice.
    /*! If the i-th element is explicitly stored, a copy of this value is returned. Otherwise, 0.0
     *  will be returned. Unlike with the []-operator, the data structure remains unchanged either way.
     *  Thus this operator should always be used for read-only access to the elements of a sparse vector slice.
     */
    const complex operator()(const int i) const {
#if(CXSC_INDEX_CHECK)
      if(i<lb || i>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector_slice::operator()(const int)"));
#endif
      complex r(0.0);

      for(int k=start ; k<end && p[k]-start<=i-lb ; k++) {
        if(p[k]-start == i-lb) 
          r = x[k];
      }

      return r; 
    }

    //! Assigns v to all elements of the vector slice 
    scvector_slice& operator=(const real& v) {
      return sl_vs_assign<scvector_slice,real,complex,std::vector<complex>::iterator>(*this,v);
    }

    //! Assigns v to all elements of the vector slice
    scvector_slice& operator=(const complex& v) {
      return sl_vs_assign<scvector_slice,complex,complex,std::vector<complex>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    scvector_slice& operator=(const srvector_slice& v) {
      return slsl_vv_assign<scvector_slice,srvector_slice,complex,std::vector<complex>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    scvector_slice& operator=(const scvector_slice& v) {
      return slsl_vv_assign<scvector_slice,scvector_slice,complex,std::vector<complex>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    scvector_slice& operator=(const srvector& v) {
      return slsp_vv_assign<scvector_slice,srvector,complex,std::vector<complex>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    scvector_slice& operator=(const scvector& v) {
      return slsp_vv_assign<scvector_slice,scvector,complex,std::vector<complex>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    scvector_slice& operator=(const rvector& v) {
      return slf_vv_assign<scvector_slice,rvector,complex,std::vector<complex>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    scvector_slice& operator=(const cvector& v) {
      return slf_vv_assign<scvector_slice,cvector,complex,std::vector<complex>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    scvector_slice& operator=(const rvector_slice& v) {
      return slf_vv_assign<scvector_slice,rvector_slice,complex,std::vector<complex>::iterator>(*this,v);
    }

    //! Overwrite the vector slice with the elements of v
    scvector_slice& operator=(const cvector_slice& v) {
      return slf_vv_assign<scvector_slice,cvector_slice,complex,std::vector<complex>::iterator>(*this,v);
    }

    //! Operator for multiplication with a scalar, result is assigned to the vector slice
    scvector_slice& operator*=(const real& s) {
      return sl_vs_multassign(*this,s);
    }

    //! Operator for multiplication with a scalar, result is assigned to the vector slice
    scvector_slice& operator*=(const complex& s) {
      return sl_vs_multassign(*this,s);
    }

    //! Operator for division of each element of the vector slice with a scalar, result is assigned to the vector slice
    scvector_slice& operator/=(const real& s) {
      return sl_vs_divassign(*this,s);
    }

    //! Operator for division of each element of the vector slice with a scalar, result is assigned to the vector slice
    scvector_slice& operator/=(const complex& s) {
      return sl_vs_divassign(*this,s);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    scvector_slice& operator+=(const rvector& v) {
      return slf_vv_addassign<scvector_slice,rvector,complex>(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    scvector_slice& operator+=(const cvector& v) {
      return slf_vv_addassign<scvector_slice,cvector,complex>(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    scvector_slice& operator+=(const rvector_slice& v) {
      return slf_vv_addassign<scvector_slice,rvector_slice,complex>(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    scvector_slice& operator+=(const cvector_slice& v) {
      return slf_vv_addassign<scvector_slice,cvector_slice,complex>(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    scvector_slice& operator+=(const srvector& v) {
      return slsp_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    scvector_slice& operator+=(const scvector& v) {
      return slsp_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    scvector_slice& operator+=(const srvector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    //! Operator for element-wise addition with a vector, result is assigned to the vector slice
    scvector_slice& operator+=(const scvector_slice& v) {
      return slsl_vv_addassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    scvector_slice& operator-=(const rvector& v) {
      return slf_vv_subassign<scvector_slice,rvector,complex>(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    scvector_slice& operator-=(const cvector& v) {
      return slf_vv_subassign<scvector_slice,cvector,complex>(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    scvector_slice& operator-=(const rvector_slice& v) {
      return slf_vv_subassign<scvector_slice,rvector_slice,complex>(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    scvector_slice& operator-=(const cvector_slice& v) {
      return slf_vv_subassign<scvector_slice,cvector_slice,complex>(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    scvector_slice& operator-=(const srvector& v) {
      return slsp_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    scvector_slice& operator-=(const scvector& v) {
      return slsp_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    scvector_slice& operator-=(const srvector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }

    //! Operator for element-wise subtraction with a vector, result is assigned to the vector slice
    scvector_slice& operator-=(const scvector_slice& v) {
      return slsl_vv_subassign(*this,v);
    }

    friend int Lb(const scvector_slice&);
    friend int Ub(const scvector_slice&);
    friend srvector Re(const scvector_slice&);
    friend srvector Im(const scvector_slice&);
    friend int VecLen(const scvector_slice&);

//     friend srvector operator*(const srmatrix&, const srvector_slice&); //ok
//     friend srvector operator*(const srmatrix_slice&, const srvector_slice&); //ok

    friend class srvector;
    friend class scvector;
    friend class scivector;
    friend class cvector;
    friend class cvector_slice;
    friend class civector;
    friend class civector_slice;

#include "vector_friend_declarations.inl"
};

inline cvector::cvector(const scvector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new complex[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline cvector::cvector(const srvector_slice& v) {
  l = v.lb;
  u = v.ub;
  size = v.n;
  dat = new complex[v.n];
  for(int i=0 ; i<v.n ; i++)
    dat[i] = 0.0;
  for(int i=v.start ; i<=v.end ; i++)
    dat[v.p[i]] = v.x[i];
}

inline cvector_slice& cvector_slice::operator=(const srvector& v) {
  *this = rvector(v);
  return *this;
}

inline cvector_slice& cvector_slice::operator=(const srvector_slice& v) {
  *this = rvector(v);
  return *this;
}

inline cvector_slice& cvector_slice::operator=(const scvector& v) {
  *this = cvector(v);
  return *this;
}

inline cvector_slice& cvector_slice::operator=(const scvector_slice& v) {
  *this = cvector(v);
  return *this;
}

inline scvector::scvector(const srvector_slice& s) : lb(s.lb), ub(s.ub), n(s.n)  {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(complex(s.x[i]));
  }

}

inline scvector::scvector(const scvector_slice& s) : lb(s.lb), ub(s.ub), n(s.n) {
  p.reserve(s.nnz);
  x.reserve(s.nnz);

  for(int i=s.start ; i<=s.end ; i++) {
    p.push_back(s.p[i]-s.offset);
    x.push_back(s.x[i]);
  }

}

inline scvector& scvector::operator=(const srvector_slice& v) {
  return spsl_vv_assign<scvector,srvector_slice,complex>(*this,v);
}

inline scvector& scvector::operator=(const scvector_slice& v) {
  return spsl_vv_assign<scvector,scvector_slice,complex>(*this,v);
}

inline scvector_slice scvector::operator()(const int i, const int j) {
#if(CXSC_INDEX_CHECK)
  if(i<lb || j>ub) cxscthrow(ELEMENT_NOT_IN_VEC("scvector_slice::operator()(const int, const int)"));
#endif
  return scvector_slice(*this,i,j);
}

//! Returns the vector -v
inline scvector operator-(const scvector_slice& v) {
  return sl_v_negative<scvector_slice,scvector>(v);
}

//! Returns the lower index bound of the vector slice v
inline int Lb(const scvector_slice& v) {
  return v.lb;
}

//! Returns the upper index bound of the vector slice v
inline int Ub(const scvector_slice& v) {
  return v.ub;
}

//! Returns the real part of the complex vector slice
inline srvector Re(const scvector_slice& v) {
  return Re(scvector(v));
}

//! Returns the imaginary part of the complex vector slice
inline srvector Im(const scvector_slice& v) {
  return Im(scvector(v));
}

//! Returns the length of the vector slice
inline int VecLen(const scvector_slice& v) {
  return v.n;
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector_slice& v1, const rvector& v2) {
  return slf_vv_mult<scvector_slice,rvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const srvector_slice& v1, const cvector& v2) {
  return slf_vv_mult<srvector_slice,cvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector_slice& v1, const cvector& v2) {
  return slf_vv_mult<scvector_slice,cvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const cvector& v1, const srvector_slice& v2) {
  return fsl_vv_mult<cvector,srvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const rvector& v1, const scvector_slice& v2) {
  return fsl_vv_mult<rvector,scvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const cvector& v1, const scvector_slice& v2) {
  return fsl_vv_mult<cvector,scvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_mult<scvector_slice,rvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const srvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_mult<srvector_slice,cvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_mult<scvector_slice,cvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const cvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_mult<cvector_slice,srvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const rvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_mult<rvector_slice,scvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const cvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_mult<cvector_slice,scvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector& v1, const srvector_slice& v2) {
  return spsl_vv_mult<scvector,srvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const srvector& v1, const scvector_slice& v2) {
  return spsl_vv_mult<srvector,scvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector& v1, const scvector_slice& v2) {
  return spsl_vv_mult<scvector,scvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector_slice& v1, const srvector& v2) {
  return slsp_vv_mult<scvector_slice,srvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const srvector_slice& v1, const scvector& v2) {
  return slsp_vv_mult<srvector_slice,scvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector_slice& v1, const scvector& v2) {
  return slsp_vv_mult<scvector_slice,scvector,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_mult<scvector_slice,srvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const srvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_mult<srvector_slice,scvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Computes the dot product v1*v2
/**
 * Note that the precision used for the computation is set by the global variable opdotprec. 
 * By default it is set to 0, meaning maximum accuracy (this is also the slowest option).
 * To use standard floating point operations, set opdotprec=1. Setting opdotprec to K>=2
 * uses (simulated) K-fold double precision.
 */
inline complex operator*(const scvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_mult<scvector_slice,scvector_slice,complex,sparse_cdot>(v1,v2);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline scvector operator*(const scvector_slice& v, const real& s) {
  return sp_vs_mult<scvector_slice,real,scvector>(v,s);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline scvector operator*(const scvector_slice& v, const complex& s) {
  return sp_vs_mult<scvector_slice,complex,scvector>(v,s);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline scvector operator*(const srvector_slice& v, const complex& s) {
  return sp_vs_mult<srvector_slice,complex,scvector>(v,s);
}

//! Divides all elements of v by the scalar s and returns the result as a new vector
inline scvector operator/(const scvector_slice& v, const real& s) {
  return sp_vs_div<scvector_slice,real,scvector>(v,s);
}

//! Divides all elements of v by the scalar s and returns the result as a new vector
inline scvector operator/(const scvector_slice& v, const complex& s) {
  return sp_vs_div<scvector_slice,complex,scvector>(v,s);
}

//! Divides all elements of v by the scalar s and returns the result as a new vector
inline scvector operator/(const srvector_slice& v, const complex& s) {
  return sp_vs_div<srvector_slice,complex,scvector>(v,s);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline scvector operator*(const real& s, const scvector_slice& v) {
  return sp_sv_mult<real,scvector_slice,scvector>(s,v);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline scvector operator*(const complex& s, const scvector_slice& v) {
  return sp_sv_mult<complex,scvector_slice,scvector>(s,v);
}

//! Multiplies all elements of v with the scalar s and returns the result as a new vector
inline scvector operator*(const complex& s, const srvector_slice& v) {
  return sp_sv_mult<complex,srvector_slice,scvector>(s,v);
}

//! Element-wise addition of v1 and v2
inline cvector operator+(const cvector& v1, const srvector_slice& v2) {
  return fsl_vv_add<cvector,srvector_slice,cvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline cvector operator+(const rvector& v1, const scvector_slice& v2) {
  return fsl_vv_add<rvector,scvector_slice,cvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline cvector operator+(const cvector& v1, const scvector_slice& v2) {
  return fsl_vv_add<cvector,scvector_slice,cvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline cvector operator+(const scvector_slice& v1, const rvector& v2) {
  return slf_vv_add<scvector_slice,rvector,cvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline cvector operator+(const srvector_slice& v1, const cvector& v2) {
  return slf_vv_add<srvector_slice,cvector,cvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline cvector operator+(const scvector_slice& v1, const cvector& v2) {
  return slf_vv_add<scvector_slice,cvector,cvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline cvector operator+(const cvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_add<cvector_slice,srvector_slice,cvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline cvector operator+(const rvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_add<rvector_slice,scvector_slice,cvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline cvector operator+(const cvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_add<cvector_slice,scvector_slice,cvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline cvector operator+(const scvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_add<scvector_slice,rvector_slice,cvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline cvector operator+(const srvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_add<srvector_slice,cvector_slice,cvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline cvector operator+(const scvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_add<scvector_slice,cvector_slice,cvector>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline scvector operator+(const scvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_add<scvector_slice,srvector_slice,scvector,complex>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline scvector operator+(const srvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_add<srvector_slice,scvector_slice,scvector,complex>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline scvector operator+(const scvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_add<scvector_slice,scvector_slice,scvector,complex>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline scvector operator+(const scvector& v1, const srvector_slice& v2) {
  return spsl_vv_add<scvector,srvector_slice,scvector,complex>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline scvector operator+(const srvector& v1, const scvector_slice& v2) {
  return spsl_vv_add<srvector,scvector_slice,scvector,complex>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline scvector operator+(const scvector& v1, const scvector_slice& v2) {
  return spsl_vv_add<scvector,scvector_slice,scvector,complex>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline scvector operator+(const scvector_slice& v1, const srvector& v2) {
  return slsp_vv_add<scvector_slice,srvector,scvector,complex>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline scvector operator+(const srvector_slice& v1, const scvector& v2) {
  return slsp_vv_add<srvector_slice,scvector,scvector,complex>(v1,v2);
}

//! Element-wise addition of v1 and v2
inline scvector operator+(const scvector_slice& v1, const scvector& v2) {
  return slsp_vv_add<scvector_slice,scvector,scvector,complex>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline cvector operator-(const cvector& v1, const srvector_slice& v2) {
  return fsl_vv_sub<cvector,srvector_slice,cvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline cvector operator-(const rvector& v1, const scvector_slice& v2) {
  return fsl_vv_sub<rvector,scvector_slice,cvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline cvector operator-(const cvector& v1, const scvector_slice& v2) {
  return fsl_vv_sub<cvector,scvector_slice,cvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline cvector operator-(const scvector_slice& v1, const rvector& v2) {
  return slf_vv_sub<scvector_slice,rvector,cvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline cvector operator-(const srvector_slice& v1, const cvector& v2) {
  return slf_vv_sub<srvector_slice,cvector,cvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline cvector operator-(const scvector_slice& v1, const cvector& v2) {
  return slf_vv_sub<scvector_slice,cvector,cvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline cvector operator-(const cvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_sub<cvector_slice,srvector_slice,cvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline cvector operator-(const rvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_sub<rvector_slice,scvector_slice,cvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline cvector operator-(const cvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_sub<cvector_slice,scvector_slice,cvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline cvector operator-(const scvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_sub<scvector_slice,rvector_slice,cvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline cvector operator-(const srvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_sub<srvector_slice,cvector_slice,cvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline cvector operator-(const scvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_sub<scvector_slice,cvector_slice,cvector>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline scvector operator-(const scvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_sub<scvector_slice,srvector_slice,scvector,complex>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline scvector operator-(const srvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_sub<srvector_slice,scvector_slice,scvector,complex>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline scvector operator-(const scvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_sub<scvector_slice,scvector_slice,scvector,complex>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline scvector operator-(const scvector& v1, const srvector_slice& v2) {
  return spsl_vv_sub<scvector,srvector_slice,scvector,complex>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline scvector operator-(const srvector& v1, const scvector_slice& v2) {
  return spsl_vv_sub<srvector,scvector_slice,scvector,complex>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline scvector operator-(const scvector& v1, const scvector_slice& v2) {
  return spsl_vv_sub<scvector,scvector_slice,scvector,complex>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline scvector operator-(const scvector_slice& v1, const srvector& v2) {
  return slsp_vv_sub<scvector_slice,srvector,scvector,complex>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline scvector operator-(const srvector_slice& v1, const scvector& v2) {
  return slsp_vv_sub<srvector_slice,scvector,scvector,complex>(v1,v2);
}

//! Element-wise subtraction of v1 and v2
inline scvector operator-(const scvector_slice& v1, const scvector& v2) {
  return slsp_vv_sub<scvector_slice,scvector,scvector,complex>(v1,v2);
}

inline cvector& cvector::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline cvector& cvector::operator+=(const scvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator+=(const srvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator+=(const scvector_slice& v2) {
  return fsl_vv_addassign(*this,v2);
}

inline scvector& scvector::operator+=(const srvector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline scvector& scvector::operator+=(const scvector_slice& v2) {
  return spsl_vv_addassign(*this,v2);
}

inline cvector& cvector::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline cvector& cvector::operator-=(const scvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator-=(const srvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline cvector_slice& cvector_slice::operator-=(const scvector_slice& v2) {
  return fsl_vv_subassign(*this,v2);
}

inline scvector& scvector::operator-=(const srvector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

inline scvector& scvector::operator-=(const scvector_slice& v2) {
  return spsl_vv_subassign(*this,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const scvector_slice& v1, const srvector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const scvector_slice& v1, const scvector_slice& v2) {
  return slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const scvector_slice& v1, const srvector& v2) {
  return slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector_slice& v1, const scvector& v2) {
  return slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const scvector_slice& v1, const scvector& v2) {
  return slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const scvector& v1, const srvector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector& v1, const scvector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const scvector& v1, const scvector_slice& v2) {
  return spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const scvector_slice& v1, const rvector& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector_slice& v1, const cvector& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const scvector_slice& v1, const cvector& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const cvector& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const rvector& v1, const scvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const cvector& v1, const scvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const scvector_slice& v1, const rvector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const srvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const scvector_slice& v1, const cvector_slice& v2) {
  return slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const cvector_slice& v1, const srvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const rvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are equal to the respective elements of v2.
 */
inline bool operator==(const cvector_slice& v1, const scvector_slice& v2) {
  return fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector_slice& v1, const srvector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector_slice& v1, const scvector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector_slice& v1, const scvector_slice& v2) {
  return !slsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector_slice& v1, const rvector& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector_slice& v1, const cvector& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector_slice& v1, const cvector& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const cvector& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const rvector& v1, const scvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const cvector& v1, const scvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector_slice& v1, const srvector& v2) {
  return !slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector_slice& v1, const scvector& v2) {
  return !slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector_slice& v1, const scvector& v2) {
  return !slsp_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector& v1, const srvector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector& v1, const scvector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector& v1, const scvector_slice& v2) {
  return !spsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector_slice& v1, const rvector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const srvector_slice& v1, const cvector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const scvector_slice& v1, const cvector_slice& v2) {
  return !slf_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const cvector_slice& v1, const srvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const rvector_slice& v1, const scvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Element-wise comparison of v1 and v2
/** 
 * Returns true only if all elements of v1 are not equal to the respective elements of v2.
 */
inline bool operator!=(const cvector_slice& v1, const scvector_slice& v2) {
  return !fsl_vv_comp(v1,v2);
}

//! Output operator for sparse vector slice v.
/**
 *  The output format is set by global flags, default is dense output.
 *  Use cout << SparseInOut; for sparse output or cout << MatrixMarketInOut; for
 *  output in matrix market format.
 */
inline std::ostream& operator<<(std::ostream& os, const scvector_slice& v) {
  return sl_v_output<scvector_slice,complex>(os,v);
}

//! Input operator for sparse vector slice v.
/**
 *  The input format is set by global flags, default is dense input.
 *  Use cout << SparseInOut; for sparse input or cout << MatrixMarketInOut; for
 *  input in matrix market format.
 */
inline std::istream& operator>>(std::istream& is, scvector_slice& v) {
  return sl_v_input<scvector_slice,complex>(is,v);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector& x, const scvector& y) {
  spsp_vv_accu<cdotprecision,scvector,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector& x, const srvector& y) {
  spsp_vv_accu<cdotprecision,scvector,srvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector& x, const scvector& y) {
  spsp_vv_accu<cdotprecision,srvector,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector& x, const cvector& y) {
  spf_vv_accu<cdotprecision,scvector,cvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector& x, const rvector& y) {
  spf_vv_accu<cdotprecision,scvector,rvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector& x, const cvector& y) {
  spf_vv_accu<cdotprecision,srvector,cvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector& x, const cvector_slice& y) {
  spf_vv_accu<cdotprecision,scvector,cvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector& x, const rvector_slice& y) {
  spf_vv_accu<cdotprecision,scvector,rvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector& x, const cvector_slice& y) {
  spf_vv_accu<cdotprecision,srvector,cvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const cvector& x, const scvector& y) {
  fsp_vv_accu<cdotprecision,cvector,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const cvector& x, const srvector& y) {
  fsp_vv_accu<cdotprecision,cvector,srvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const rvector& x, const scvector& y) {
  fsp_vv_accu<cdotprecision,rvector,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const cvector_slice& x, const scvector& y) {
  fsp_vv_accu<cdotprecision,cvector_slice,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const cvector_slice& x, const srvector& y) {
  fsp_vv_accu<cdotprecision,cvector_slice,srvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const rvector_slice& x, const scvector& y) {
  fsp_vv_accu<cdotprecision,rvector_slice,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector_slice& x, const cvector& y) {
  slf_vv_accu<cdotprecision,scvector_slice,cvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector_slice& x, const rvector& y) {
  slf_vv_accu<cdotprecision,scvector_slice,rvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector_slice& x, const cvector& y) {
  slf_vv_accu<cdotprecision,srvector_slice,cvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector_slice& x, const cvector_slice& y) {
  slf_vv_accu<cdotprecision,scvector_slice,cvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector_slice& x, const rvector_slice& y) {
  slf_vv_accu<cdotprecision,scvector_slice,rvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector_slice& x, const cvector_slice& y) {
  slf_vv_accu<cdotprecision,srvector_slice,cvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const cvector& x, const scvector_slice& y) {
  fsl_vv_accu<cdotprecision,cvector,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const cvector& x, const srvector_slice& y) {
  fsl_vv_accu<cdotprecision,cvector,srvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const rvector& x, const scvector_slice& y) {
  fsl_vv_accu<cdotprecision,rvector,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const cvector_slice& x, const scvector_slice& y) {
  fsl_vv_accu<cdotprecision,cvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const cvector_slice& x, const srvector_slice& y) {
  fsl_vv_accu<cdotprecision,cvector_slice,srvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const rvector_slice& x, const scvector_slice& y) {
  fsl_vv_accu<cdotprecision,rvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector_slice& x, const scvector_slice& y) {
  slsl_vv_accu<cdotprecision,scvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector_slice& x, const srvector_slice& y) {
  slsl_vv_accu<cdotprecision,scvector_slice,srvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector_slice& x, const scvector_slice& y) {
  slsl_vv_accu<cdotprecision,srvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector& x, const scvector_slice& y) {
  spsl_vv_accu<cdotprecision,scvector,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector& x, const srvector_slice& y) {
  spsl_vv_accu<cdotprecision,scvector,srvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector& x, const scvector_slice& y) {
  spsl_vv_accu<cdotprecision,srvector,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector_slice& x, const scvector& y) {
  slsp_vv_accu<cdotprecision,scvector_slice,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const scvector_slice& x, const srvector& y) {
  slsp_vv_accu<cdotprecision,scvector_slice,srvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cdotprecision& dot, const srvector_slice& x, const scvector& y) {
  slsp_vv_accu<cdotprecision,srvector_slice,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector& x, const scvector& y) {
  spsp_vv_accuapprox<cdotprecision,scvector,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector& x, const srvector& y) {
  spsp_vv_accuapprox<cdotprecision,scvector,srvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector& x, const scvector& y) {
  spsp_vv_accuapprox<cdotprecision,srvector,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector& x, const cvector& y) {
  spf_vv_accuapprox<cdotprecision,scvector,cvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector& x, const rvector& y) {
  spf_vv_accuapprox<cdotprecision,scvector,rvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector& x, const cvector& y) {
  spf_vv_accuapprox<cdotprecision,srvector,cvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector& x, const cvector_slice& y) {
  spf_vv_accuapprox<cdotprecision,scvector,cvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector& x, const rvector_slice& y) {
  spf_vv_accuapprox<cdotprecision,scvector,rvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector& x, const cvector_slice& y) {
  spf_vv_accuapprox<cdotprecision,srvector,cvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const cvector& x, const scvector& y) {
  fsp_vv_accuapprox<cdotprecision,cvector,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const cvector& x, const srvector& y) {
  fsp_vv_accuapprox<cdotprecision,cvector,srvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const rvector& x, const scvector& y) {
  fsp_vv_accuapprox<cdotprecision,rvector,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const cvector_slice& x, const scvector& y) {
  fsp_vv_accuapprox<cdotprecision,cvector_slice,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const cvector_slice& x, const srvector& y) {
  fsp_vv_accuapprox<cdotprecision,cvector_slice,srvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const rvector_slice& x, const scvector& y) {
  fsp_vv_accuapprox<cdotprecision,rvector_slice,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const cvector& y) {
  slf_vv_accuapprox<cdotprecision,scvector_slice,cvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const rvector& y) {
  slf_vv_accuapprox<cdotprecision,scvector_slice,rvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const cvector& y) {
  slf_vv_accuapprox<cdotprecision,srvector_slice,cvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const cvector_slice& y) {
  slf_vv_accuapprox<cdotprecision,scvector_slice,cvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const rvector_slice& y) {
  slf_vv_accuapprox<cdotprecision,scvector_slice,rvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const cvector_slice& y) {
  slf_vv_accuapprox<cdotprecision,srvector_slice,cvector_slice,sparse_cdot>(dot,x,y);
}

inline void accumulate_approx(cdotprecision& dot, const cvector& x, const scvector_slice& y) {
  fsl_vv_accuapprox<cdotprecision,cvector,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const cvector& x, const srvector_slice& y) {
  fsl_vv_accuapprox<cdotprecision,cvector,srvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const rvector& x, const scvector_slice& y) {
  fsl_vv_accuapprox<cdotprecision,rvector,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const cvector_slice& x, const scvector_slice& y) {
  fsl_vv_accuapprox<cdotprecision,cvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const cvector_slice& x, const srvector_slice& y) {
  fsl_vv_accuapprox<cdotprecision,cvector_slice,srvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const rvector_slice& x, const scvector_slice& y) {
  fsl_vv_accuapprox<cdotprecision,rvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const scvector_slice& y) {
  slsl_vv_accuapprox<cdotprecision,scvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const srvector_slice& y) {
  slsl_vv_accuapprox<cdotprecision,scvector_slice,srvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const scvector_slice& y) {
  slsl_vv_accuapprox<cdotprecision,srvector_slice,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector& x, const scvector_slice& y) {
  spsl_vv_accuapprox<cdotprecision,scvector,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector& x, const srvector_slice& y) {
  spsl_vv_accuapprox<cdotprecision,scvector,srvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector& x, const scvector_slice& y) {
  spsl_vv_accuapprox<cdotprecision,srvector,scvector_slice,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const scvector& y) {
  slsp_vv_accuapprox<cdotprecision,scvector_slice,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const scvector_slice& x, const srvector& y) {
  slsp_vv_accuapprox<cdotprecision,scvector_slice,srvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object. This version does not compute an error bound if the precision is not equal
 * to 0. This is faster, but no reliable enclosure of the computed result can be given.
 */
inline void accumulate_approx(cdotprecision& dot, const srvector_slice& x, const scvector& y) {
  slsp_vv_accuapprox<cdotprecision,srvector_slice,scvector,sparse_cdot>(dot,x,y);
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector& x, const srvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector& x, const cvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector& x, const rvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector& x, const cvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector& x, const cvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector& x, const rvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector& x, const cvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const cvector& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const cvector& x, const srvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const cvector_slice& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const cvector_slice& x, const srvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector_slice& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector_slice& x, const cvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector_slice& x, const rvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector_slice& x, const cvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector_slice& x, const cvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector_slice& x, const rvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector_slice& x, const cvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const cvector& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const cvector& x, const srvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const cvector_slice& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const cvector_slice& x, const srvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const rvector_slice& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector_slice& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector_slice& x, const srvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector_slice& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector& x, const srvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector& x, const scvector_slice& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector_slice& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const scvector_slice& x, const srvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

//! The accurate scalar product of the last two arguments added to the value of the first argument.
/**
 * The precision for the dotproduct can be set by calling the set_dotprec member function of the dotprecision object.
 */
inline void accumulate(cidotprecision& dot, const srvector_slice& x, const scvector& y) {
  cdotprecision tmp(0.0);
  tmp.set_k(dot.get_k());
  accumulate(tmp,x,y);
  dot += tmp;
}

} //namespace cxsc

#include "sparsevector.inl"

#endif
