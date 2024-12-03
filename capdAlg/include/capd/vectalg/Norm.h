/// @addtogroup vectalg
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Norm.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include <string>

#ifndef _CAPD_VECTALG_NORM_H_ 
#define _CAPD_VECTALG_NORM_H_ 

namespace capd{
namespace vectalg{

/**
 * A general abstract norm
 */
template<typename VectorType,typename MatrixType>
class Norm{
public:
  typedef typename VectorType::ScalarType ScalarType;

  /// Computes the norm of the vector
  virtual ScalarType operator()(const VectorType &iv) const = 0;
  /// Computes the norm of the matrix
  virtual ScalarType operator()(const MatrixType &iv) const = 0;
  /// Returns the name of the norm
  virtual std::string show(void) const = 0;
  /// Makes virtual copy of the object
  virtual Norm *clone(void) const = 0;
  /// Destructor
  virtual ~Norm() {}
};
/**
 *  Euclidean norm
 */
template<typename VectorType,typename MatrixType>
class EuclNorm : public Norm<VectorType,MatrixType>{
public:
  typedef typename VectorType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  virtual ScalarType operator()(const VectorType &iv) const;
  virtual ScalarType operator()(const MatrixType &iv) const;
  virtual std::string show(void) const;
  Norm<VectorType,MatrixType> * clone(void) const;
};
/**
 *  \f$ L_\infty \f$ norm (max norm)
 */
template<typename VectorType,typename MatrixType>
class MaxNorm : public Norm<VectorType,MatrixType>{
public:
  typedef typename VectorType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  virtual ScalarType operator()(const VectorType &iv) const;
  virtual ScalarType operator()(const MatrixType &iv) const;
  virtual std::string show(void) const;
  Norm<VectorType,MatrixType> * clone(void) const;
};
/**
 * \f$ L_1 \f$ norm 
 */
template<typename VectorType,typename MatrixType>
class SumNorm : public Norm<VectorType,MatrixType>{
public:
  typedef typename VectorType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  virtual ScalarType operator()(const VectorType &iv) const;
  virtual ScalarType operator()(const MatrixType &iv) const;
  virtual std::string show(void) const;
  Norm<VectorType,MatrixType> * clone(void) const;
};
/**
 * Euclidean Logarithmic Norm
 */
template<typename VectorType,typename MatrixType>
class EuclLNorm : public Norm<VectorType,MatrixType>{
public:
  typedef typename VectorType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  virtual ScalarType operator()(const VectorType &iv) const;
  virtual ScalarType operator()(const MatrixType &iv) const;
  virtual std::string show(void) const;
  Norm<VectorType,MatrixType> * clone(void) const;
};

/**
 *  \f$ L_\infty \f$  Logarithmic Norm
 */
template<typename VectorType,typename MatrixType>
class MaxLNorm : public Norm<VectorType,MatrixType>{
public:
  typedef typename VectorType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  virtual ScalarType operator()(const VectorType &iv) const;
  virtual ScalarType operator()(const MatrixType &iv) const;
  virtual std::string show(void) const;
  Norm<VectorType,MatrixType> * clone(void) const;
};

/**
 * \f$ L_1 \f$ Logarytmic norm
 */
template<typename VectorType,typename MatrixType>
class SumLNorm : public Norm<VectorType,MatrixType>{
public:
  typedef typename VectorType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  virtual ScalarType operator()(const VectorType &iv) const;
  virtual ScalarType operator()(const MatrixType &iv) const;
  virtual std::string show(void) const;
  Norm<VectorType,MatrixType> * clone(void) const;
};

}} // namespace capd::vectalg

#endif // _CAPD_VECTALG_NORM_H_ 

/// @}
