#include <boost/test/unit_test.hpp>

#include "capd/capdlib.h"
#include "capd/poincare/TimeMap.hpp"
#include "capd/basicalg/TypeTraits.h"
#include <iomanip>
using namespace capd;

// Test for rigorous computations (sets)
template <class SetT, bool is_interval  = capd::TypeTraits<typename SetT::ScalarType>::isInterval  >
class SetProperties{
  public:
    static void check(SetT & set){
      typedef typename SetT::VectorType Vector;
      Vector solution = (Vector)set;
      Vector enclosure = set.getLastEnclosure();
      BOOST_WARN_MESSAGE(subset(solution,enclosure), set.name() << std::setprecision(20) << " : solution is not subset of enclosure \n solution  : " << solution 
                                                                                                          << "\n enclosure : " << enclosure);
      BOOST_REQUIRE_MESSAGE(!intersectionIsEmpty(solution,enclosure), set.name() << std::setprecision(20) << " : solution and enclosure have empty intersection\n solution  : " << solution 
                                                                                                          << "\n enclosure : " << enclosure);
    //  solution = set.currentSet();
    //  BOOST_REQUIRE_MESSAGE(subset(solution,enclosure), "solution is not subset of solution \n solution" << solution << " enclosure : " << enclosure); 
    }
};

// Test for nonrigorous sets (vectors)
template <class SetT>
class SetProperties<SetT, false>{
  public:
  static void check(SetT&){
  }
};
template<class TM, class V>
void checkVariableTime(TM& tm, V& s){
  tm.stopAfterStep(true);
  typename TM::ScalarType T = 0.;
  do{
    tm(1.,s);
    BOOST_REQUIRE(T<tm.getCurrentTime());
    BOOST_REQUIRE(tm.getSolver().getStep()>0);
    T += tm.getSolver().getStep();
    BOOST_REQUIRE(T==tm.getCurrentTime());
    SetProperties<V>::check(s);
  }while(!tm.completed());
  BOOST_REQUIRE(tm.getSolver().getStep()>0);
  BOOST_REQUIRE(tm.getSolver().getStep()!=1);
}

template<class TM, class V>
void checkFixedTime(TM& tm, V& x, double fixedTimeStep){
  tm.stopAfterStep(true);
  double s = 1./fixedTimeStep;
  BOOST_REQUIRE_MESSAGE(s!=(int)s,"the test has incorrect settings. Fixed time step cannot divide the time of integration.");
  typename TM::ScalarType T = 0.;
  do{
    tm(1.,x);
    BOOST_REQUIRE(T<tm.getCurrentTime());
    if(!tm.completed()){
      BOOST_REQUIRE(tm.getSolver().getStep()==fixedTimeStep);
    } else {
      BOOST_REQUIRE(tm.getSolver().getStep()<fixedTimeStep);
    }
    T += tm.getSolver().getStep();
    BOOST_REQUIRE_MESSAGE(T==tm.getCurrentTime(), "T=" << T << ", tm.getCurrentTime()=" << tm.getCurrentTime());
    SetProperties<V>::check(x);
  }while(!tm.completed());
  BOOST_REQUIRE(tm.getSolver().getStep()>0);
  BOOST_REQUIRE(tm.getSolver().getStep()<fixedTimeStep);
}
