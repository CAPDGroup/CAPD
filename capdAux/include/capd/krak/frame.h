
/////////////////////////////////////////////////////////////////////////////
/// @file frame.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

// Documentation added by Mikolaj Zalewski, June 2000.
#ifndef _CAPD_KRAK_FRAME_H_
#define _CAPD_KRAK_FRAME_H_

#include "capd/krak/IFrame.h"

#include "capd/krak/primitiv.h"
#include "capd/krak/usermove.h"
#include "capd/krak/manip.h"


namespace capd{
namespace krak{

/// @addtogroup krak
/// @{

  class Frame: public IFrame {
public:
   int ltj,lti,rbj,rbi,cj,ci;
   int imarg,jmarg;
   double swx,swy,nex,ney;
   int cRow,cCol,lRow,lCol;
   int bgColor,fgColor;
   int prec;

   void initFrm(
         int arglti,int argltj,int argrbi, int argrbj,
         int bgc=WHITE, int fgc=BLACK,
         int im=fontHght/2, int jm=fontWdth/2
   );

   Frame(void);
   Frame(
         int lti,int ltj,int rbi,int rbj,
         int bgc=WHITE,int fgc=BLACK,
         int im=fontWdth/2, int jm=fontHght/2
   );

   Frame(
         const Frame &prntFrm,const capd::krak::At &lt,const capd::krak::At &rb,
         int bgc=WHITE,int fgc=BLACK,
         int im=fontWdth/2, int jm=fontHght/2
   );

   Frame(
         Frame &prntFrm,int lti,int ltj,int rbi,int rbj,
         int bgc=WHITE,int fgc=BLACK,
         int im=fontHght/2, int jm=fontWdth/2
   );

   Frame(
         Frame &prntFrm,int lti,int ltj,int rbi,int rbj,
         double swx,double swy,double nex,double ney,
         int bgc=WHITE,int fgc=BLACK,
         int im=fontHght/2, int jm=fontWdth/2
   );

   void Locate(Frame &prntFrm, capd::krak::At &lt, capd::krak::At &rb);
   void Locate(Frame &prntFrm,int lti,int ltj,int rbi,int rbj);
   void setWorldCoord(double swx,double swy,double nex,double ney);

   void adjust(void);

   int NoCol(void);
   int NoRow(void);

   void Clear(void);
   void Clear(int color);
   void Bound(int color=BLACK);


   void SetBgColor(int c);
   void SetFgColor(int c);

   int getRow(capd::krak::Pxl &pxl);
   int getCol(capd::krak::Pxl &pxl);

   int x2i(double x);
   int y2j(double y);

   double x2p(double x);
   double y2q(double y);

   double i2x(int i);
   double j2y(int j);

   void jump(int i,int j);
   void draw(int i,int j,int color=FRAME_FG);
   void drawText(const char *c,int i,int j,int color=FRAME_FG);
   void dot(int i,int j,int color=FRAME_FG);
   void circle(int i,int j,int r, int color=FRAME_FG);
   void line(int i1,int j1,int i2,int j2,int color=FRAME_FG);
   void box(int lti,int ltj,int rbi,int rbj,int color=FRAME_FG);
   void boxFill(int lti,int ltj,int rbi,int rbj,int col,int pattern=SOLID_P);
#if defined (WIN95) || defined (WXWIN)
   void polygon(int coords[],int nPoints,int color=FRAME_FG);
   void polygonFill(int coords[],int nPoints,int col,int pattern=SOLID_P);
   void arc(int lti,int ltj,int rbi,int rbj,int bi,int bj,int ei,int ej,int color=FRAME_FG);
   void arcFill(int lti,int ltj,int rbi,int rbj,
      int bi,int bj,int ei,int ej,int col,int pattern=SOLID_P
   );
   void ellipse(int lti,int ltj,int rbi,int rbj,int color=FRAME_FG);
   void ellipseFill(int lti,int ltj,int rbi,int rbj,int col,int pattern=SOLID_P);
#endif

   void jump(double x,double y);
   void draw(double x,double y,int color=FRAME_FG);
   void drawText(const char *c,double x,double y,int color=FRAME_FG);
   void dot(double x,double y,int color=FRAME_FG);
   void circle(double x,double y,int r, int color=FRAME_FG);
   void line(double x1,double y1,double x2,double y2, int color=FRAME_FG);
   void box(double swx,double swy, double nex, double ney,int color=FRAME_FG);
   void boxFill(double swx,double swy, double nex, double ney,int col,int pattern=SOLID_P);
#if defined (WIN95) || defined (WXWIN)
   void polygon(double coords[],int nPoints,int color=FRAME_FG);
   void polygonFill(double coords[],int nPoints,int col,int pattern=SOLID_P);
   void arc(double swx,double swy, double nex, double ney,
      double bx,double by,double ex,double ey,int color=FRAME_FG
   );
   void arcFill(double swx,double swy, double nex, double ney,
      double bx,double by,double ex,double ey,int col,int pattern=SOLID_P
   );
   void ellipse(double swx,double swy, double nex, double ney,int color=FRAME_FG);
   void ellipseFill(double swx,double swy, double nex, double ney,int col,int pattern=SOLID_P);
#endif

   void Xcrss(double x,double y,int size=1, int color=FRAME_FG);

   int precision(int p){prec = p; return p;}

   Frame &operator<<(char c);
   Frame &operator<<(int n);
   Frame &operator<<(long n);
   Frame &operator<<(double r);

   Frame &operator<<(const capd::krak::frstring &a_string);
   Frame &operator<<(const char *text);
   Frame &operator<<(const capd::krak::colstring &a_colstring);

   Frame &operator<<(const capd::krak::FgColor &c);
   Frame &operator<<(const capd::krak::BgColor &c);
   friend Frame &operator<<(Frame &f, const capd::krak::At &at);
   Frame &operator<<(capd::krak::Tab tab);


   Frame &operator>>(const capd::krak::At &at);
   Frame &operator>>(const capd::krak::FgColor &c);
   Frame &operator>>(const capd::krak::BgColor &c);

   int isInside(capd::krak::Pxl &p);

   Frame &operator>>(capd::krak::frstring &a_string);
   Frame &operator>>(int &n);
   Frame &operator>>(long &n);
   Frame &operator>>(double &r);
   int Edit(At at, int no_col, capd::krak::frstring &s);
};

/// @}
}} // the end of the namespace capd::krak

//###################### TEMPLATES DEFINITIONS ####################################


// A universal template for outputting anything to a frame via a string.
namespace capd{
namespace krak{
template <class type>
capd::krak::Frame &operator << (capd::krak::Frame &f, const type &x)
{
   std::ostringstream s;
   s.precision(f.prec);
   s << x << std::ends;
   f << s.str().c_str();
   return f;
} /* operator << */
}}

// ###############################  inline definitions ########################


inline void capd::krak::Frame::adjust(void){
   lti = lti/capd::krak::fontWdth * capd::krak::fontWdth;
   ltj = ltj/capd::krak::fontHght * capd::krak::fontHght;
   rbi = rbi/capd::krak::fontWdth * capd::krak::fontWdth;
   rbj = rbj/capd::krak::fontHght * capd::krak::fontHght;
   lRow= (rbj - ltj-jmarg-jmarg)/capd::krak::fontHght-1;
   lCol= (rbi - lti-imarg-imarg)/capd::krak::fontWdth-1;
}

/**
  Return the number of columns in the frame
*/
inline int capd::krak::Frame::NoCol(void)
{
   return (rbi - lti - imarg - imarg)/capd::krak::fontWdth;
}

/**
  Return the number of rows in the frame
*/
inline int capd::krak::Frame::NoRow(void)
{
  return (rbj - ltj - jmarg - jmarg)/capd::krak::fontHght;
}

/**
  Sets to world coordinates.
  @param swx,swy the coordinates of the bottom-left (southwest) corner
  @param nex,ney the coordinates of the upper-right (northeast) corner
*/
inline void capd::krak::Frame::setWorldCoord(double swx,double swy,
                                 double nex,double ney)
{
   dscrFrm(this,swx,swy,nex,ney);
}

/**
  Clears the frame with the background color and moves
  the current position to (0,0)
*/
inline void capd::krak::Frame::Clear(void)
{
   capd::krak::Rct r;

   SetRct(&r,lti,ltj,rbi,rbj);
   FillRct(&r,SOLID_P,bgColor);
   cRow=cCol=0;
}

/**
  Clears the frame with the color \c color and moves
  the current position to (0,0)
*/
inline void capd::krak::Frame::Clear(int color)
{
   capd::krak::Rct r;

   SetRct(&r,lti,ltj,rbi,rbj);
   FillRct(&r,SOLID_P,color);
   cRow=cCol=0;
}

/**
  Draw a boundary around the frame in color \c color
*/
inline void capd::krak::Frame::Bound(int color)
{
   SetFgCol(color);
   Rctngl(lti,ltj,rbi,rbj);
}

/**
  Sets the background color to \c c.
*/

inline void capd::krak::Frame::SetBgColor(int c)
{
   bgColor=c;
}

/**
  Sets the foreground color to \c c.
*/
inline void capd::krak::Frame::SetFgColor(int c)
{
   fgColor=c;
}

/**
  returns the character row that correspond to the pixel \c pxl
*/
inline int capd::krak::Frame::getRow(capd::krak::Pxl &pxl)
{
   return (pxl.j-ltj-jmarg)/capd::krak::fontHght;
}

/**
  returns the character column that correspond to the pixel \c pxl
*/
inline int capd::krak::Frame::getCol(capd::krak::Pxl &pxl)
{
   return (pxl.i-lti-imarg)/capd::krak::fontWdth;
}

/**
Changes the foreground color of the frame, like:
<pre>
frm<< FgColor(RED)<<"red "<< FgColor(BLUE)<<"blue"
</pre>
*/
inline capd::krak::Frame &capd::krak::Frame::operator<<(const capd::krak::FgColor &c)
{
   SetFgColor(c.color);
   return *this;
}

/**
Changes the background color of the frame, like:
<pre>
frm<< BgColor(YELLOW)<<"yellow"<< BgColor(GREEN)<<"green"
</pre>
*/
inline capd::krak::Frame &capd::krak::Frame::operator<<(const capd::krak::BgColor &c)
{
   SetBgColor(c.color);
   return *this;
}

namespace capd{
namespace krak{

/**
  Moves the current position to the cell refered by \c at, like:
<pre>
frm<<At(30,30)<<"AAAAA";
</pre>
*/
inline capd::krak::Frame &operator<<(capd::krak::Frame &f, const capd::krak::At &at)
{
   f.cCol=at.col;
   f.cRow=at.row;
   return f;
}
}}

/**
  Moves the current position to the column refered by \c tab
*/
inline capd::krak::Frame &capd::krak::Frame::operator<<(capd::krak::Tab tab)
{
   cCol=tab.col;
   return *this;
}

/**
  Translates world coordinate \c x to a device coordinate
*/
inline int capd::krak::Frame::x2i(double x)
{
   return (int)(((x-swx)/(nex-swx))*(rbi-lti)+lti);
}

/**
  Translates world coordinate \c y to a device coordinate
*/
inline int capd::krak::Frame::y2j(double y)
{
   return (int)((1-(y-swy)/(ney-swy))*(rbj-ltj)+ltj);
}

/**
  Same as \ref x2i() but returns a double
*/
inline double capd::krak::Frame::x2p(double x)
{
   return (((x-swx)/(nex-swx))*(rbi-lti)+lti);
}

/**
  Same as \ref y2j() but returns a double
*/
inline double capd::krak::Frame::y2q(double y)
{
   return ((1-(y-swy)/(ney-swy))*(rbj-ltj)+ltj);
}

/**
  Translates device coordinate \c i to a world coordinate
*/
inline double capd::krak::Frame::i2x(int i)
{
   return (double)(((double)i-lti)/(rbi-lti))*(nex-swx)+swx;
}

/**
  Translates device coordinate \c j to a world coordinate
*/
inline double capd::krak::Frame::j2y(int j)
{
   return (double)(1-((double)j-ltj)/(rbj-ltj))*(ney-swy)+swy;
}

/**
  Moves the current position to (\c i, \c j), in device coordinates
*/
inline void capd::krak::Frame::jump(int i,int j)
{
   ci=i;cj=j;
}

/**
  Moves the current position to (\c x, \c y), in world coordinates
*/
inline void capd::krak::Frame::jump(double x,double y)
{
   jump(x2i(x),y2j(y));
}

/**
  Draws a line from the current position to (\c x, \c y) (in world
  coordinates) using the color \c color
*/
inline void capd::krak::Frame::draw(double x,double y,int color)
{
   draw(x2i(x),y2j(y),color);
}

/**
  Draws a text from the position (\c x, \c y) (in world
  coordinates) using the color \c color
*/
inline void capd::krak::Frame::drawText(const char *c,double x,double y,int color)
{
   drawText(c,x2i(x),y2j(y),color);
}

/**
  draws a box (an empty rectangle) with the left top corner in
  (\c lti, \c ltj) and the right botton corner in (\c rbi, \c rbj), in
  device coordinates.
@param color The color of the box. If it is FRAME_FG then the frame
  foreground color is used
*/
inline void capd::krak::Frame::box(int lti,int ltj,int rbi,int rbj,int color)
{
   jump(lti,ltj);
   draw(lti,rbj,color);
   draw(rbi,rbj,color);
   draw(rbi,ltj,color);
   draw(lti,ltj,color);
}

/**
  draws a box (an empty rectangle) with the left botton corner in
  (\c swx, \c swy) and the right top corner in (\c nex, \c ney), in
  world coordinates.
@param color The color of the box. If it is FRAME_FG then the frame
  foreground color is used
*/
inline void capd::krak::Frame::box(double swx,double swy, double nex, double ney, int color)
{
   box(x2i(swx),y2j(swy),x2i(nex),y2j(ney),color);
}

/**
  Draws a filled box with the botton left corner at (\c swx, \c swy) and
  the right top corner at (\c nex, \c ney). The box is filled with color
  \c color
*/
inline void capd::krak::Frame::boxFill(int lti,int ltj,int rbi,int rbj,int col,int pattern)
{
   capd::krak::Rct r;
   SetRct(&r,lti,ltj,rbi,rbj);
   FillRct(&r,pattern,col);
}

/**
  Draws a filled box with the botton left corner at (\c swx, \c swy) and
  the right top corner at (\c nex, \c ney). The box is filled with color
  \c color
*/
inline void capd::krak::Frame::boxFill(
   double swx,double swy, double nex, double ney,int col,int pattern)
{
   capd::krak::Rct r;
#ifndef LINUX
   SetRct(&r,x2i(swx),y2j(swy),x2i(nex),y2j(ney));
#else
   SetRct(&r,x2i(swx),y2j(ney),x2i(nex),y2j(swy));
#endif
   FillRct(&r,pattern,col);
}

#if defined (WIN95) || defined (WXWIN)
/**
  Draws an arc of an ellipse bound by the rectangle of world coordinates coordinates
  (\c swx,\c swy) and (\c nex,\c ney). The arc is indicated by the points
  (\c bx,\c by) and (\c ex,\c ey). The arc is drawn in the color i@(color)
*/
inline void capd::krak::Frame::arc(double swx,double swy, double nex, double ney,
      double bx,double by,double ex,double ey,int color)
{
   arc(x2i(swx),y2j(swy),x2i(nex),y2j(ney),x2i(bx),y2j(by),x2i(ex),y2j(ey),color);
}

/**
  Fills and area between an arc and its chord. The arch is an arc
  of an ellipse  bound by the rectangle of world coordinates coordinates
  (\c swx,\c swy) and (\c nex,\c ney). The arc is indicated by the points
  (\c bx,\c by) and (\c ex,\c ey). The arc is drawn in the color i@(color)
*/
inline void capd::krak::Frame::arcFill(double swx,double swy, double nex, double ney,
      double bx,double by,double ex,double ey,int color,int pattern)
{
   arcFill(x2i(swx),y2j(swy),x2i(nex),y2j(ney),x2i(bx),y2j(by),x2i(ex),y2j(ey),color,pattern);
}


/**
  Draws an ellipse  bound by the rectangle of world coordinates coordinates
  (\c lti,\c ltj) and (\c rbi,\c rbj). The ellipse is drawn in the color i@(color)
*/

inline void capd::krak::Frame::ellipse(int lti,int ltj,int rbi,int rbj,int color)
{
   arc(lti,ltj,rbi,rbj,lti,ltj,lti,ltj,color);
}

/**
  Draws an ellipse  bound by the rectangle of world coordinates coordinates
  (\c lti,\c ltj) and (\c rbi,\c rbj). The ellipse is filled with the
  pattern i@(pattern) in the color i@(color)
*/

inline void capd::krak::Frame::ellipseFill(int lti,int ltj,int rbi,int rbj,int color,int pattern)
{
   arcFill(lti,ltj,rbi,rbj,lti,ltj,lti,ltj,color,pattern);
}

/**
  Draws an ellipse  bound by the rectangle of world coordinates coordinates
  (\c swx,\c swy) and (\c nex,\c ney). The ellipse is drawn in the color i@(color)
*/

inline void capd::krak::Frame::ellipse(double swx,double swy, double nex, double ney,int color)
{
   arc(swx,swy,nex,ney,swx,swy,swx,swy,color);
}

/**
  Draws an ellipse  bound by the rectangle of world coordinates coordinates
  (\c swx,\c swy) and (\c nex,\c ney). The ellipse is drawn in the color i@(color)
  pattern i@(pattern) in the color i@(color)
*/

inline void capd::krak::Frame::ellipseFill(double swx,double swy, double nex, double ney,
      int color,int pattern)
{
   // This function is defined ONLY for the integer coordinates!
// arcFill(swx,swy,nex,ney,swx,swy,swx,swy,color,pattern);
   arcFill (x2i(swx), y2j(swy), x2i(nex), y2j(ney), x2i(swx), y2j(swy),
      x2i(swx), y2j(swy), color, pattern);
}
#endif

/**
  The same as operator<< - changes the current position to \c at
*/
inline capd::krak::Frame &capd::krak::Frame::operator>>(const capd::krak::At &at)
{
   cCol=at.col;
   cRow=at.row;
   return *this;
}

/**
  The same as operator<< - changes the foreground color to \c c
*/
inline capd::krak::Frame &capd::krak::Frame::operator>>(const capd::krak::FgColor &c)
{
   SetFgColor(c.color);
   return *this;
}

/**
  The same as operator<< - changes the background color to \c c
*/
inline capd::krak::Frame &capd::krak::Frame::operator>>(const capd::krak::BgColor &c)
{
   SetBgColor(c.color);
   return *this;
}

/**
  Checks if the pixel \c pxl (in device coordinates) is inside the frame
*/
inline int capd::krak::Frame::isInside(capd::krak::Pxl &p)
{
   return ( p.i >= lti && p.i < rbi && p.j >= ltj && p.j< rbj );
}

/*___________________________________________________________________________*/

namespace capd{
namespace krak{
inline capd::krak::frstring &operator<<(capd::krak::frstring &out, capd::krak::Tab &t)
{
   capd::krak::frstring f;
   f.resize(t.col);
   for(int i=0;i<t.col;i++) f.str[i]=' ';
   out=out^f;
   return out;
}

}} // the end of the namespace capd::krak

#endif // _CAPD_KRAK_FRAME_H_
