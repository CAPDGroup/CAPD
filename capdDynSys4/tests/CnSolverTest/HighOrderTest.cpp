#include <boost/test/unit_test.hpp>
#include "capd/capdlib.h"
#include "capd/dynset/C1GraphicalSet.h"
using namespace capd;
BOOST_AUTO_TEST_SUITE(HighOrderTest)
      
template <typename T>
class TimeCount{
  public:
    TimeCount() : numberOfMoves(0), totalTime(0.){
    }
    void reset(){
      numberOfMoves = 0;
      totalTime = capd::TypeTraits<T>::zero();
      
    }
//    template <typename SetType>
    void show(C1Set & set){
     totalTime = set.getCurrentTime();
     ++ numberOfMoves;
    }
    friend std::ostream & operator << (std::ostream & str, const TimeCount & t){
      str << "total time = " << t.totalTime << "  number of steps = " << t.numberOfMoves ;
      return str;
    }
  protected:
    int numberOfMoves; 
    T totalTime;
    
//    static int totalNumberOfMoves;
//    static T totalTime;
};
   
bool testCnTimeMap(int order, interval time, IVector initial_conditions){
  
  IMap vf("var:x,y,px,py;par:c;fun:px,py,-2*x^2*x-c*x*y^2,-2*y^2*y-c*x^2*y;",3);
  vf.setParameter("c",interval(-0.9));
  ICnOdeSolver  solver(vf, order);
  ICnTimeMap tmap(solver);
  typedef capd::dynset::C1GraphicalSet<IMatrix, TimeCount<interval> > SetType;
  //CnRect2Set set(initial_conditions, 3);
  C1Rect2Set set(initial_conditions);
  TimeCount<interval> output;
  SetType gset(set, output);
  IVector value = tmap(time, gset);
  std::cout << "\n result " << value << "\n time spend : " << output << std::endl ;
  return true;
}    
bool testCnPoincareMap(int order, IVector initial_conditions){
  
  IMap vf("var:x,y,px,py;par:c;fun:px,py,-2*x^2*x-c*x*y^2,-2*y^2*y-c*x^2*y;",3);
  vf.setParameter("c",interval(-0.9));
  ICnOdeSolver  solver(vf, order);
  
  ICoordinateSection section(4,1);
  ICnPoincareMap pm(solver, section, capd::poincare::MinusPlus);
  CnRect2Set set(initial_conditions, 3);
  IVector value = pm(set);
  std::cout << "\n result " << value ;
  return true;
}    

BOOST_AUTO_TEST_CASE(Order80Test){
  IVector ic("{[3.527185573233460e-1 ,3.527185573461631e-1 ],[0 ,0 ],[5.478828384836637e-1 ,5.478828385134561e-1 ],[8.272523511223313e-1 ,8.272523511444839e-1 ]}");
  interval time = 10;
  int order = 80;
  testCnTimeMap(order, time, ic);
  testCnPoincareMap(order, ic);
}

BOOST_AUTO_TEST_CASE(Order40Test){
  
  IVector ic("{[3.527185573233460e-1 ,3.527185573461631e-1 ],[0 ,0 ],[5.478828384836637e-1 ,5.478828385134561e-1 ],[8.272523511223313e-1 ,8.272523511444839e-1 ]}");
  interval time = 20;
  int order = 40;
  testCnTimeMap(order, time, ic);
  testCnPoincareMap(order, ic);
}

BOOST_AUTO_TEST_SUITE_END()
