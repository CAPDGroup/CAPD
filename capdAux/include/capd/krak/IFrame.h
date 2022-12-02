#ifndef CAPD_KRAK_IFRAME_H
#define CAPD_KRAK_IFRAME_H


#include <capd/krak/defines.h>


namespace capd{
namespace krak{


  /*
    Interface class for a frame.
    With that interface we can develop a code with dependencies to GUI and do not require GUI libraries.
    TODO add more methods.
   */
  class IFrame {

  public:
    virtual void draw(int i,int j,int color=FRAME_FG) = 0;
    virtual void line(int i1,int j1,int i2,int j2,int color=FRAME_FG) = 0;
    virtual void Xcrss(double x,double y,int size=1, int color=FRAME_FG) = 0;
    virtual void box(int lti,int ltj,int rbi,int rbj,int color=FRAME_FG) = 0;
    virtual void boxFill(int lti,int ltj,int rbi,int rbj,int col,int pattern=SOLID_P) = 0;
    virtual void jump(int i,int j) = 0;

    virtual void setWorldCoord(double swx,double swy,double nex,double ney) = 0;
    virtual void jump(double x,double y) = 0;
    virtual void draw(double x,double y,int color=FRAME_FG) = 0;
    virtual void drawText(const char *c,double x,double y,int color=FRAME_FG) = 0;
    virtual void dot(double x,double y,int color=FRAME_FG) = 0;
    virtual void circle(double x,double y,int r, int color=FRAME_FG) = 0;
    virtual void line(double x1,double y1,double x2,double y2, int color=FRAME_FG) = 0;
    virtual void box(double swx,double swy, double nex, double ney,int color=FRAME_FG) = 0;
    virtual void boxFill(double swx,double swy, double nex, double ney,int col,int pattern=SOLID_P) = 0;
};
}
}


#endif
