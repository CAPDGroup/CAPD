/////////////////////////////////////////////////////////////////////////////
/// @file CoordWiseReorganization.h
///
/// @author kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_COORDWISEREORGANIZATION_H_
#define _CAPD_DYNSET_COORDWISEREORGANIZATION_H_

#include "capd/dynset/reorganization/FactorPolicy.h"
#include "capd/dynset/DefaultPolicy.h"
#include "capd/vectalg/Vector_Interval.hpp"

namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/**
 *  In this reorganization policy column vector B_i corresponding to the biggest coordinate in r (r_i)
 *  replaces vector in C_j the 'closest' to B_i  and the r_i is moved to r0
 */
template <typename BasePolicy = DefaultPolicy>
class CoordWiseReorganization : public FactorPolicy<BasePolicy>
{
public:

  CoordWiseReorganization() : m_minReorganizationSize(0.0){ }

  /// it finds index of the biggest coordinate of vector v
  template<class VectorType>
  static int findBiggestCoord(const VectorType & v) {
    typedef typename VectorType::ScalarType::BoundType BoundType;
    int best = 0;
    BoundType max = v[0].rightBound();
    for(typename VectorType::size_type i=1; i<v.dimension(); i++){
      BoundType value = v[i].rightBound();
      if(max < value){
        best = i;
        max = value;
      }
    }
    return best;
  }
  /// returns index of the column vector in base which best fits to given vector v
  template<class MatrixType, class VectorType>
  static int findBestFittingVector(const MatrixType & base, const VectorType & v) {
    VectorType vInBase = capd::matrixAlgorithms::gauss(base, v);
    return findBiggestCoord(vInBase);
  }

  /// reorganize given set
  template<class SetType>
  void reorganize(SetType & result) const{
    typedef typename SetType::ScalarType ScalarType;
    typedef typename SetType::VectorType VectorType;
       // index of coordinate with biggest diameter in r
     int rIndex = findBiggestCoord(diam(result.get_r()));
       // index of coordinate in r0 and column in C which will be replaced
     int r0Index = findBestFittingVector(result.get_C(), result.getColumn_B(rIndex));

     ScalarType removedValue = result.getElement_r0(r0Index);
     VectorType removedVector = result.getColumn_C(r0Index);

       // replacing coordinate in r0 and column in C with data from r and B
     result.setColumn_C(r0Index, result.getColumn_B(rIndex));
     result.setElement_r0(r0Index, result.getElement_r(rIndex));

       // we express removed vector in the new base and add it to r0
     removedVector = capd::matrixAlgorithms::gauss(result.get_C(), removedVector) * removedValue;
     result.set_r0(result.get_r0()+removedVector);

       // we set r[rIndex] to 0
     result.setElement_r(rIndex, ScalarType(0.0));
  }

  /// reorganize given set
  template<class Matrix, class Vector>
  void reorganize(Matrix& B, Matrix& /*invB*/, Vector& r, Matrix& C, Vector& r0) const {
    typedef typename Vector::ScalarType ScalarType;
    typedef Vector VectorType;

    // index of coordinate with biggest diameter in r
     int rIndex = findBiggestCoord(diam(r));
       // index of coordinate in r0 and column in C which will be replaced
     int r0Index = findBestFittingVector(C, B.getColumn(rIndex));

     ScalarType removedValue = r0[r0Index];
     VectorType removedVector = C.getColumn(r0Index);

       // replacing coordinate in r0 and column in C with data from r and B
     C.getColumnn(r0Index) = B.getColumn(rIndex);
     r0[r0Index]= r[rIndex];

       // we express removed vector in the new base and add it to r0
     removedVector = capd::matrixAlgorithms::gauss(C, removedVector) * removedValue;
     r0 = r0 +removedVector;

       // we set r[rIndex] to 0
     r[rIndex]= ScalarType(0.0);
  }

//
  template<class SetType>
  bool isReorganizationNeeded(const SetType & result) const{
    typedef typename SetType::ScalarType ScalarType;
    ScalarType s = capd::vectalg::maxDiam(result.get_r()).rightBound();
    return (this->isReorganizationEnabled() &&
        (s > m_minReorganizationSize) &&
        (s > this->getC0Factor() * capd::vectalg::maxDiam(result.get_invB() *( result.get_C()*result.get_r0()))));
  }

  /// makes reorganization if needed.
  /// return true if reorganization was performed
  template<class SetType>
  bool reorganizeIfNeeded(SetType & result) const{
    if(isReorganizationNeeded(result)){
      reorganize(result);
      return true;
    }
    return false;
  }

  template<class Matrix, class Vector>
  bool reorganizeIfNeeded(Matrix& B, Matrix& invB, Vector& r, Matrix& C, Vector& r0) const
  {
    if(this->isReorganizationNeeded(r,r0)){
      this->reorganize(B,invB,r,C,r0);
      return true;
    }
    return false;
  }

  std::string name() const {
    return "coordwise reorganization";
  }
protected :
  double m_minReorganizationSize;
};

/// @}
}} // end of capd::dynset

#endif //_CAPD_DYNSET_COORDWISEREORGANIZATION_H_
