

/////////////////////////////////////////////////////////////////////////////
/// @file DiscreteDynSys.h
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

#ifndef  _CAPD_DYNSYS_DISCRETEDYNSYS_H_
#define  _CAPD_DYNSYS_DISCRETEDYNSYS_H_
#include "capd/dynsys/DynSys.hpp"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

/**
 *  DiscreteDynSys is a proxy to convert any Map into discrete Dynamical System
 *
 *  This conversion is needed to use Map as dicrete dynamical systems
 *  to move various set representations.  This class implements both
 *  DynSys and C1Map interface.
 *
 *  MapT should be type implementing C1Map interface i.e.
 *  if f is of type MapT then the following expresions should be valid:
 *  f(x) - returns value of a map in x,
 *  f[x] - returns derivative of a map in x
 *  f(x, deriv) returns value of a map and in deriv its derivative.
 */
template <typename MapT>
class DiscreteDynSys : public DynSys<typename MapT::MatrixType>{
public:
	  typedef MapT MapType;
	  typedef typename MapType::MatrixType MatrixType;
	  typedef typename DynSys<MatrixType>::VectorType VectorType;
	  typedef typename DynSys<MatrixType>::ScalarType ScalarType;
	  typedef typename DynSys<MatrixType>::NormType NormType;

	  ///  map should implement C1Map concept
    DiscreteDynSys(const MapType & map) : pMap(&map){
    }

    /// value of a map for a point iv
	  VectorType Phi(const ScalarType& t, const VectorType &iv){
		return (*pMap)(t, iv);
	  }
	  /// value of a map for a point iv
	  VectorType operator()(const ScalarType& t, const VectorType &iv){
	  	return (*pMap)(t, iv);
	  }

	  /// derivative of a map in a point iv
	  MatrixType JacPhi(const ScalarType& t, const VectorType &iv){
		  return pMap->derivative(t, iv);
	  }
	  /// derivative of a map in a point iv
	  MatrixType operator[](const VectorType &iv){
	  		  return (*pMap)[iv];
	  }
    /// computes derivative of the map for a given vector and time
    MatrixType derivative(ScalarType t, const VectorType& u) const{
      return pMap->derivative(t, u);
    }                   
	  /// value and derivative of a map in a point iv
	  VectorType operator()(const ScalarType& t, const VectorType &iv, MatrixType & derivative){
	  	  	return (*pMap)(t, iv, derivative);
	  }
	  /// it computes image of the set (in give representation) using set.move function
	  template <typename SetType>
	  SetType & operator()(SetType & set) {
		    set.move(*this);
	        return set;
	  }

	  VectorType & operator()(const ScalarType& t, VectorType &iv) {
	        return (*pMap)(iv);
	  }
	  /// because we have explicit formula for a map not only numerical method
	  /// the remainder term is equal to zero
	  VectorType Remainder(const ScalarType& t, const VectorType &iv, VectorType &enc){
		  VectorType zero(iv.dimension());
		  return zero;
	  }
	  virtual void encloseC0Map(
	      const ScalarType& t,
	      const VectorType& x,
	      const VectorType& xx,
	      VectorType& o_phi,
	      VectorType& o_rem,
	      VectorType& o_enc,
	      MatrixType& o_jacPhi
	  ){
	    o_phi = (*pMap)(x);
	    o_rem.clear();
	    o_enc.clear();
	    o_jacPhi = (*pMap)[xx];
	  }
	  /// empty implemented - throws exception on use
	  ScalarType Lipschitz(const ScalarType& t, const VectorType &iv, NormType &n){
		  throw "DiscreteDynSys::Lipschitz should not be used!";
		  return ScalarType(0.0);
	  }

	  VectorType enclosure(const ScalarType& t, const VectorType& x){
      throw "DiscreteDynSys::enclosure should not be used!";
      return x;
	  }

	  /// returns name of the class (for backward compatibility)
	  virtual std::string type(void){
	    return "DiscreteDynSys";
	  }
	  virtual ~DiscreteDynSys(){}
protected:
	const MapType * pMap;
};

/// Makes DiscreteDynSys object from given map.
/// Template parameters are recognized automatically.
template <typename MapType>
DiscreteDynSys<MapType> makeDiscreteDynSys(const MapType & map){
	return DiscreteDynSys<MapType>(map);
}
/// @}
}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_DISCRETEDYNSYS_H_



