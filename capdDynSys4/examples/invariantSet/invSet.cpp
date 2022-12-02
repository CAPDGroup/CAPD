#include "capd/capdlib.h"
#include "capd/invset/MapGraph.h"
#include "capd/invset/ForbiddenSet.h"
#include "capd/invset/Scope.h"
#include "capd/invset/invariantSets.h"
#include "capd/krak/krak.h"
#ifdef  _OPENMP
#include <omp.h>
#endif
using namespace capd;
using namespace capd::invset;
using namespace capd::krak;

/*
 * We define map for which we will compute invariant set
 */
class RotatedHenonMap
{
  interval c;
public:
  // public types needed for invariant set computation 
  typedef IVector VectorType;
  typedef IMatrix MatrixType;
  
  RotatedHenonMap(interval c = 1.5) : c(c){
  }
  
  //Function is given by  f(x,y) = (c*(1-x*x)+2*x+y, -x) where c is a parameter
  IVector operator()(const IVector& u) const {
    IVector result(2);
    result[0] = c*(1-power(u[0],2))+2*u[0]+u[1] ;
    result[1] = -u[0];
    return result;
  }
  
  // The most important function: 
  // For given vector u it computes its image and adds it to result
  // Image can be a single interval vector or collection of interval vectors. 
  void eval(const IVector & u, std::list<IVector>& result) const {
    result.push_back((*this)(u));
  }
};
 
 inline void plot(Frame& f, const IVector& vect,  int color){
    f.boxFill(vect[0].leftBound(),vect[1].leftBound(),vect[0].rightBound(),vect[1].rightBound(),color);
    f.box(vect[0].leftBound(),vect[1].leftBound(),vect[0].rightBound(),vect[1].rightBound(),color);
 }

// Routine that plots graph on a screen (one dot for each vertex)
template<class GraphT>
inline void showGraph(const GraphT& G){
    //rootFrame.Clear();
    rootFrame.Clear();
    typename GraphT::VertexType denominator = GraphT::resolutionToDenominator(G.getResolution());
    for(auto vertex: G){
      // we transform vertex to coordinates of its centre 
      IVector x = GraphT::vertexToVector(vertex.first,denominator);
      plot(rootFrame, x, RED);
      // simplified plotting
      //DVector x = capd::vectalg::midObject<DVector>(GraphT::vertexToVector(vertex.first,denominator));
      //  rootFrame.dot( x[0], x[1], RED);
      
    }
}

int main(int, char**){
#ifdef  _OPENMP
  std::cout << "Number of threads : " <<  omp_get_max_threads() << " " << omp_get_num_threads() << std::endl;
#endif
  // One can set up the amount of information send to console
  capd::setLogLevel(capd::LogLevel1);
  //capd::setLogLevel(capd::Debug);
   int width = 800, height = 800, bgcolor = WHITE, fgcolor = BLACK;
  openGW(width, height, bgcolor, fgcolor);
  
  typedef invset::MapGraph<RotatedHenonMap> HenonMapGraph;
  RotatedHenonMap map;
  HenonMapGraph::ResolutionType resolution(2);
  resolution[0] = 4;  resolution[1] = 4;
  
  HenonMapGraph G(map,resolution);

 // From analytic estimates we know that invariant set is contained in [-1, 13/6] x [-13/6, 1]
  double minx = -1., maxx = 13/6;
  double miny = -13/6, maxy = 1.;
  double border = 0.1;
  
 
  // We set domain for invariant set computation
  IVector domain(2);
  domain[0] = interval(minx,maxx);
  domain[1] = interval(miny,maxy);

 
 
 // Since plotting is the most time consuming operation 
 // we can  compute invariant set in given resolution and plot it afterwards
  computeInvariantSet(G, domain, 8);
  rootFrame.setWorldCoord(minx-border,miny-border,maxx+border,maxy+border);
  showGraph(G);
 
  // We can compute invariant set, plotting it in each resolution 
 // computeInvariantSet(G, domain, 8, showGraph<HenonMapGraph>);
  
  // We can save graph to a file
  // G.save("data.txt");
 std::cout << "Done! Hit any key!\n"; 
 waitBt();
 closeGW();
  return 0;
}