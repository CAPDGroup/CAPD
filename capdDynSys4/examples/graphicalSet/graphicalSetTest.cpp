
#include <cmath>
#include <stdexcept>
#include <fstream>

#include "capd/krak/krak.h"
#include "capd/capdlib.h"

#include "capd/dynset/C1GraphicalSet.h"
#include "capd/dynset/C0GraphicalSet.h"

using namespace capd;
using capd::dynset::C0GraphicalSet;
using capd::dynset::C1GraphicalSet;

template <typename SetType>
class GraphicalOutput{
public:
  GraphicalOutput(Frame &fr, int color = GREEN, int freqency=1)
	  : m_fr(fr), m_color(color), m_freqency(freqency), m_count(0){
  }
  void show(const SetType & set){
	  if((m_count++ % m_freqency) ==0){
		  IVector w = IVector(set);
		  m_fr.boxFill(w[0].leftBound(),w[1].leftBound(),w[0].rightBound(),w[1].rightBound(),m_color);
		  m_fr.box(w[0].leftBound(),w[1].leftBound(),w[0].rightBound(),w[1].rightBound());
	  }
  }
private:
  Frame &m_fr;
  int m_color;
  int m_freqency;
  int m_count;
};

typedef GraphicalOutput<C0Set> C0Output;

double minx=-2.6;
double maxx=2.6;
double miny=-2.6;
double maxy=2.6;

Frame fr[4], txt;

void axes(Frame &fr, double minx, double miny, double maxx, double maxy)
{
   fr.setWorldCoord(minx,miny,maxx,maxy);
   fr.Clear();
   fr.line(0.0,miny,0.0,maxy,BLACK);
   fr.line(minx,0.0,maxx,0.0,BLACK);
}

// uncomment the following line if you wish to wait for a key after each step
#define _WAIT_FOR_KEY_

int readKey() {
  int s=NO_KEY;  UserMove um;
#ifdef _WAIT_FOR_KEY_
  while(NO_KEY==s)
#endif
  {   GetUserMove(um);
      s=um.key;
  }
  return s;
}

// -------------------------------------------------------------------------

void initGraph()
{
   openGW(860,610,WHITE,BLACK);
   fr[0] = Frame(5,5,400,300,WHITE,BLACK);
   fr[1] = Frame(400,5,800,300,WHITE,BLACK);
   fr[2] = Frame(5,300,400,600,WHITE,BLACK);
   fr[3] = Frame(400,300,800,600,WHITE,BLACK);
   txt = Frame(5,400,800,860,WHITE,BLACK);
   rootFrame.Clear();
   for(int i =0; i<4; i++)
    axes(fr[i],minx,miny,maxx,maxy);
}


// -------------------------------------------------------------------------
void drawRectangle(std::vector<IVector> & cor, Frame & fr, int color){
  for(int i = 0; i<4; i++){
    for(int j=0; j<4; j++){
      fr.line( cor[i][0].mid().leftBound(), cor[i][1].mid().leftBound(),cor[j][0].mid().leftBound(),cor[j][1].mid().leftBound(), color);
    }
    fr.boxFill(cor[i][0].leftBound(), cor[i][1].leftBound(),cor[i][0].rightBound(), cor[i][1].rightBound(),color,SOLID_P);
  }
}

template<typename T>
void corners(T& head, const T & tail, int i, int dim, std::vector<T> & cor)
{
  if(i < dim)
  {
    head[i] = left(tail[i]);
    corners(head, tail, i+1, dim, cor);
    head[i] = right(tail[i]);
    corners(head, tail, i+1, dim, cor);
    head[i] = tail[i];
  }
  else
  {
    cor.push_back(head);
  }
}


std::vector<IVector> getCorners(const IVector &center, const IMatrix &B, const IVector &r)
{
  typedef IVector VectorType;
  std::vector<VectorType> cor;
  VectorType v = r;
  corners(v, r, 0, v.dimension(), cor);
  for(std::vector<VectorType>::iterator it = cor.begin(); it != cor.end(); ++it){
    *it = center + B * *it;
  }
  return cor;
}



// -------------------------------------------------------------------------

void makeTests(IOdeSolver &T, const IVector & x0)
{

  GraphicalOutput<C0Set> out0(fr[0], RED);
  GraphicalOutput<C0Set> out1(fr[1], BLUE);
  GraphicalOutput<C0Set> out2(fr[2], GREEN);
  GraphicalOutput<C0Set> out3(fr[3], BLACK);
  
  
  fr[0] << "C0GraphicalSet(Intv2Set)";
  C0GraphicalSet<C0Intv2Set,C0Output> set0(C0Intv2Set(x0), out0);

  fr[1] << "C0GraphicalSet(C0RectSet)";
  C0GraphicalSet<C0RectSet,C0Output> set1(C0RectSet(x0), out1);

  fr[2] << "C0GraphicalSet(C0Rect2Set)";
  C0GraphicalSet<C0Rect2Set,C0Output> set2(C0Rect2Set(x0), out2);

  fr[3] << "C0GraphicalSet(C0TripletonSet)";
  C0GraphicalSet<C0TripletonSet,C0Output> set3(C0TripletonSet(x0), out3);
 

  int s;
  do{
    set0.move(T);
    set1.move(T);
    set2.move(T);
    set3.move(T);
    s=readKey();
  }while(s!=EscKey);
  rootFrame.Clear();
  rootFrame << "C0GraphicalSet(Intv2Set)     size = " << maxWidth(IVector(set0)) << "\n\n"
            << "C0GraphicalSet(C0RectSet)    size = " << maxWidth(IVector(set1)) << "\n\n"
            << "C0GraphicalSet(C0Rect2Set)   size = " << maxWidth(IVector(set2)) << "\n\n"
            << "C1GraphicalSet(C1Rect2Set)   size = " << maxWidth(IVector(set3)) << "\n\n";
}

void harmonicOscilator(){
  int order = 5;
  double step = 6.28 / 30.;
  double sizeOfSet = 0.15;
  DVector x(2);  x[0] = 1.0;  x[1] = 0.0;


  IVector v(2);
  v[0] = interval(x[0]-sizeOfSet,x[0]+sizeOfSet);
  v[1] = interval(x[1]-sizeOfSet,x[1]+sizeOfSet);

  IMap vfield("var:x,y;fun:-y,x;");
  IOdeSolver solver(vfield,order);
  solver.setStep(step);

  makeTests(solver, v);
}
// -------------------------------------------------------------------------

int main(int, char**)
{
  initGraph();

  try
  {
    harmonicOscilator();
    waitBt();

  }catch(std::exception& e)
  {
    std::ofstream plik;
    plik.open("report");
    plik << e.what();
    plik.close();

   rootFrame << "\n\nException caught! See 'report' file for details.";
   waitBt();
  }
  closeGW();
  return 0;
}
