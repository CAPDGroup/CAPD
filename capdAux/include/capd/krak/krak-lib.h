
/////////////////////////////////////////////////////////////////////////////
/// @file krak-lib.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

/*
   This is the include file for the compiled krak library
   in a simplified verion for educational purposes
   Author: Marian Mrozek
*/

/* Basic color definitions */

#define WHITE       0
#define BLACK       1
#define RED         2
#define GREEN       3
#define BLUE        4
#define YELLOW      5
#define MAGENTA     6
#define CYAN        7
#define ORANGE      8
#define VIOLET      9
#define PINE       10
#define BROWN      11
#define OLIVE      12
#define DARKBLUE   13
#define ORANGERED  14
#define BLUEGREEN  15

#define FRAME_FG   -1
#define FRAME_BG   -2

/* Pattern definitions */

#define MAX_PATTERN_NO 17

#define EMPTY_P     0
#define SOLID_P     1
#define HLINE_P     2
#define VLINE_P     3
#define DHLINE_P    4
#define DVLINE_P    5
#define DOT_P       6
#define DDOT_P      7
#define DUST_P      8
#define DDUST_P     9
#define SLASH_P    10
#define ISLASH_P   11
#define WHDLINE_P  12
#define DWHDLINE_P 13
#define HASH_P     14
#define DHASH_P    15
#define WHASH_P    16
#define DWHASH_P   17
#define WVLINE_P   18
#define VDOT_P     19


#define NO_KEY 0x1ff
#define DRAG_KEY 0x1fd
#define BUTTON_KEY 0x1fe


namespace capd{
namespace krak{

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

class At{
public:
   int row,col;
   At(int r,int c);
};

class Tab{
public:
   int col;
   Tab(int c);
};

class FgColor{
public:
   int color;
   FgColor(int c);
};

class BgColor{
public:
   int color;
   BgColor(int c);
};

void SetBgCol(int col);
void SetFgCol(int col);

}} // the end of the namespace capd::krak

// ###############################  inline definitions ########################

inline capd::krak::At::At(int r,int c)
{
   row=r;col=c;
}

inline capd::krak::Tab::Tab(int c)
{
   col=c;
}

inline capd::krak::FgColor::FgColor(int c)
{
   color=c;
}

inline capd::krak::BgColor::BgColor(int c)
{
   color=c;
}

namespace capd{
namespace krak{
  
extern int fontHght,fontWdth;

/* primitive class Frame */

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

   Frame(void);
   Frame(
         int lti,int ltj,int rbi,int rbj,
         int bgc=WHITE,int fgc=BLACK,
         int im=fontWdth/2, int jm=fontHght/2
   );

   void setWorldCoord(double swx,double swy,double nex,double ney);

   void Clear(void);
   void Clear(int color);
   void Bound(int color=BLACK);

   void SetBgColor(int c);
   void SetFgColor(int c);

   void jump(int i,int j);
   void draw(int i,int j,int color=FRAME_FG);
   void drawText(const char *c,int i,int j,int color=FRAME_FG);
   void dot(int i,int j,int color=FRAME_FG);
   void circle(int i,int j,int r, int color=FRAME_FG);
   void line(int i1,int j1,int i2,int j2,int color=FRAME_FG);
   void box(int lti,int ltj,int rbi,int rbj,int color=FRAME_FG);
   void boxFill(int lti,int ltj,int rbi,int rbj,int col,int pattern=SOLID_P);

   void jump(double x,double y);
   void draw(double x,double y,int color=FRAME_FG);
   void drawText(const char *c,double x,double y,int color=FRAME_FG);
   void dot(double x,double y,int color=FRAME_FG);
   void circle(double x,double y,int r, int color=FRAME_FG);
   void line(double x1,double y1,double x2,double y2, int color=FRAME_FG);
   void box(double swx,double swy, double nex, double ney,int color=FRAME_FG);
   void boxFill(double swx,double swy, double nex, double ney,int col,int pattern=SOLID_P);

   void Xcrss(double x,double y,int size=1, int color=FRAME_FG);

   int precision(int p){prec = p; return p;}

   Frame &operator<<(char c);
   Frame &operator<<(int n);
   Frame &operator<<(long n);
   Frame &operator<<(double r);

   Frame &operator<<(const char *text);

   Frame &operator<<(const capd::krak::FgColor &c);
   Frame &operator<<(const capd::krak::BgColor &c);
   friend Frame &operator<<(Frame &f, const capd::krak::At &at);
   Frame &operator<<(capd::krak::Tab tab);


   Frame &operator>>(const capd::krak::At &at);
   Frame &operator>>(const capd::krak::FgColor &c);
   Frame &operator>>(const capd::krak::BgColor &c);

   Frame &operator>>(int &n);
   Frame &operator>>(long &n);
   Frame &operator>>(double &r);
};


// ###############################  inline definitions ########################

/**
Changes the foreground color of the frame, like:
<pre>
frm<< FgColor(RED)<<"red "<< FgColor(BLUE)<<"blue"
</pre>
*/
inline Frame &Frame::operator<<(const FgColor &c){
  SetFgColor(c.color);
  return *this;
};

/**
Changes the background color of the frame, like:
<pre>
frm<< BgColor(YELLOW)<<"yellow"<< BgColor(GREEN)<<"green"
</pre>
*/
inline Frame &Frame::operator<<(const BgColor &c){
  SetBgColor(c.color);
  return *this;
};

/**
  Moves the current position to the cell refered by \e at , like:
<pre>
frm<<At(30,30)<<"AAAAA";
</pre>
*/
inline Frame &operator<<(Frame &f,const  At &at){
  f.cCol=at.col;
  f.cRow=at.row;
  return f;
};

/**
  Moves the current position to the column refered by \e tab 
*/
inline Frame &Frame::operator<<(const Tab tab){
  cCol=tab.col;
  return *this;
};


/**
  Sets the background color to \e c .
*/
inline void Frame::SetBgColor(int c){
  bgColor=c;
};

/**
  Sets the foreground color to \e c .
*/
inline void Frame::SetFgColor(int c){
  fgColor=c;
};

int Button(void);

inline int  button(void){
  return Button();
};


double Clock(void);
void delay(double t);


enum FunctKeys{BSKey=8,TabKey,CRKey=13,PgUpKey,PgDnKey,EndKey,HomeKey,
    LeftKey,UpKey,RightKey,DownKey,
    InsKey,DelKey,EscKey=27,
          F1Key=89,F2Key,F3Key,F4Key,F5Key,F6Key,F7Key,F8Key,F9Key,
    DragKey=DRAG_KEY,ButtonKey=BUTTON_KEY,NoKey=NO_KEY};


/* Graphic window routines */

void closeGraphics(void);
void openGraphics(int hrs,int vrs,int bgcol,int fgcol,int ltx=100, int lty=100);
void clear(int color=WHITE);

/* Basic drawing routines in device coordinates */

void jump(int i,int j);
void draw(int i,int j,int color=FRAME_FG);
void drawText(const char *c,int i, int j, int  color=FRAME_FG);
void dot(int i,int j,int color=FRAME_FG);
void circ(int pixelRow,int pixelColumn,int r, int color=FRAME_FG);
void circFill(int pixelRow,int pixelColumn,int r, int color=FRAME_FG,int pattern=SOLID_P);
void box(int lti,int ltj,int rbi,int rbj,int color=FRAME_FG);
void boxFill(int lti,int ltj,int rbi,int rbj,int color=FRAME_FG,int pattern=SOLID_P);
void polygon(int coords[],int nPoints,int color=FRAME_FG);
void polygonFill(int coords[],int nPoints,int color=FRAME_FG,int pattern=SOLID_P);
void arc(int leftTopPixelRow,int leftTopPixelColumn,
         int rightBottomPixelRow,int rightBottomPixelColumn,
         int begPixelRow,int begPixelColumn,int endPixelRow,int endPixelColumn,
         int color=FRAME_FG);
void arcFill(int leftTopPixelRow,int leftTopPixelColumn,
         int rightBottomPixelRow,int rightBottomPixelColumn,
         int begPixelRow,int begPixelColumn,int endPixelRow,int endPixelColumn,
         int color=FRAME_FG,int pattern=SOLID_P);
void ellipse(int leftTopPixelRow,int leftTopPixelColumn,
         int rightBottomPixelRow,int rightBottomPixelColumn,int color=FRAME_FG);
void ellipseFill(int leftTopPixelRow,int leftTopPixelColumn,
         int rightBottomPixelRow,int rightBottomPixelColumn,int color=FRAME_FG,int pattern=SOLID_P);

void keyAndMouse(int &key, int &row, int &col);
void mouse(int &row, int &col);

/* Basic drawing routines in world coordinates */

void setWorldCoord(double leftBottomX,double leftBottomY,
         double rightTopX,double rightTopY);
void jump(double x,double y);
void draw(double x,double y,int color=FRAME_FG);
void drawText(const char *c,double x, double y, int  color=FRAME_FG);
void dot(double x,double y,int color=FRAME_FG);
void circ(double x,double y,double r, int color=FRAME_FG);
void circFill(double x,double y,double r, int color=FRAME_FG,int pattern=SOLID_P);
void box(double leftBottomX,double leftBottomY,
         double rightTopX,double rightTopY,int color=FRAME_FG);
void boxFill(double leftBottomX,double leftBottomY,
             double rightTopX,double rightTopY,int color=FRAME_FG,int pattern=SOLID_P);
void polygon(double coords[],int nPoints,int color=FRAME_FG);
void polygonFill(double coords[],int nPoints,int color=FRAME_FG,int pattern=SOLID_P);
void arc(double leftBottomX,double leftBottomY,
         double rightTopX,double rightTopY,
         double begX,double begY,double endX,double endY,
         int color=FRAME_FG);
void arcFill(double leftBottomX,double leftBottomY,
         double rightTopX,double rightTopY,
         double begX,double begY,double endX,double endY,
         int color=FRAME_FG,int pattern=SOLID_P);
void ellipse(double swx,double swy, double nex, double ney,int color=FRAME_FG);
void ellipseFill(double swx,double swy, double nex, double ney,int color=FRAME_FG,int pattern=SOLID_P);

void keyAndMouse(int& key, double& x,double& y);
void mouse(double x,double y);

/* Basic mouse and keybord routines */

int  button(void);
void waitButton(void);

void setTextSize(int size);
int getTextSize(void);
void setBackgroundColor(int color);
void setForegroundColor(int color);

extern Frame rootFrame;
extern Frame *rootFrm;
#define gout rootFrame
#define main() mainEntry(int,char**)

}// end of namespace krak
}// end of namespace capd

using namespace capd::krak;
