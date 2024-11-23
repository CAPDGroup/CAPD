/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.5.4)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2014 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* CVS $Id: cimath.cpp,v 1.31 2014/01/30 17:23:43 cxsc Exp $ */

/*
**
**  COmplex interval STandard functions LibrarY, CoStLy Version 1.0.3
**
**  Copyright (C) Markus Neher,        markus.neher@math.uni-karlsruhe.de
**                Ingo Eble,           ingoeble@web.de
**                Frithjof Blomquist,  Blomquist@math.uni-wuppertal.de
**
**  The complex interval elementary functions in C-XSC are based on the 
**  CoStLy library written by Markus Neher. Additional improvements have
**  been done by Frithjof Blomquist.
**  
**  References:
**  - Neher, M: "Complex Standard Functions and their Implementation in
**    the CoStLy Library", Preprint Nr. 04/18, Fakultaet fuer Mathematik,
**    Universitaet Karlsruhe, 2004.
**  - Blomquist, F.; Hofschuster, W.; Kraemer, W.: "Complex Interval Functions 
**    in C-XSC", Preprint BUW-WRSWT 2005/2, Universitaet Wuppertal, 2005.
**  - Neher, M: "Complex Standard Functions and their Implementation in
**    the CoStLy Library", to appear in ACM Transactions on Mathematical 
**    Software (TOMS).
**
*/


//Include header files
#include "cimath.hpp"  // Declaration of the complex functions
#include "rmath.hpp"   // "real" standard functions
#include "imath.hpp"   // "interval" standard functions
#include "dot.hpp"     // "dotprecision" standard functions

namespace cxsc{

const real Inf_Pi = 7074237752028440.0 / 2251799813685248.0; // Inf_Pi < Pi
// succ(Inf_Pi) > Pi is guaranteed!

const real pi2 = 4967757600021510.0 / 40564819207303340847894502572032.0;
const interval Pi2 = interval(pi2,succ(pi2)); 
// Inf_Pi+pi2 < pi;  Inf_Pi+succ(pi2) > pi;


  inline const interval& ZERO_INTERVAL()
    {
      static const interval zero = interval(0.0);
      return zero;
    }

  inline const interval& ONE_INTERVAL()
    {
      static const interval one = interval(1.0);
      return one;
    }

  inline const interval& PI()  // Blomquist;
    {
      static const interval pi = interval(Inf_Pi,succ(Inf_Pi));
      return pi;
    }

  inline const interval& HALFPI()  // Blomquist;
    {
      interval hp = PI();
      times2pown(hp,-1); // division by 2
      static const interval hpi = hp;
      return hpi;
    }

bool disjoint(const interval& x,const interval& y)
{
  real ix( Inf(x) ),iy( Inf(y) ),sx( Sup(x) ), sy( Sup(y) );
  real inf( ( ix > iy )? ix : iy );
  real sup( ( sx < sy )? sx : sy );

  return ( inf > sup );
}

//-- SQRT_2 --- INV_SQRT_2 --- LN(4) ------------------------------------------
//
const real Inf_rev_sqrt_2 = 6369051672525772.0 / 9007199254740992.0; 
// Inf_rev_sqrt_2 < 1/sqrt(2); 

inline const interval& INV_SQRT_2()
{ // optimal inclusion of 1/sqrt(2);
    static const interval t = interval(Inf_rev_sqrt_2,succ(Inf_rev_sqrt_2));
    return t;
}

const real Inf_U_atan = 6243314768165359.0 / 4503599627370496.0;
const interval U_atan(Inf_U_atan,succ(Inf_U_atan)); // optimal inclusion of 
// ln(4);


//
//-- end SQRT_2 and INV_SQRT_2 ------------------------------------------------


/* ***************************************************************************/
/* ***************************************************************************/
/* ***                      Single-valued functions                      *** */
/* ***************************************************************************/
/* ***************************************************************************/


/* ***************************************************************************/
/* *** Power operator  pow  is not listed here, since it relies on the    ****/
/* *** (multi-valued) logarithm                                           ****/
/* ***************************************************************************/


/* ***************************************************************************/
/* *** The hyperbolic functions exp, sin, cos, sinh, cosh are separable:  ****/
/* *** Their real and imaginary parts are products of real functions      ****/
/* ***************************************************************************/
/* ***   With Re(z)=x, Im(z)=y :                                          ****/
/* ***                                                                    ****/
/* ***        exp   :   Re(exp(z)) = exp(x) * cos(y)                      ****/
/* ***                  Im(exp(z)) = exp(x) * sin(y)                      ****/
/* ***                                                                    ****/
/* ***        sin   :   Re(sin(z)) = sin(x) * cosh(y)                     ****/
/* ***                  Im(sin(x)) = cos(x) * sinh(y)                     ****/
/* ***                                                                    ****/
/* ***        cos   :   Re(cos(z)) = cos(x) * cosh(y)                     ****/
/* ***                  Im(sin(x)) = -sin(x) * sinh(y)                    ****/
/* ***                                                                    ****/
/* ***        sinh  :   Re(sinh(z)) = sinh(x) * cos(y)                    ****/
/* ***                  Im(sinh(z)) = cosh(x) * sin(y)                    ****/
/* ***                                                                    ****/
/* ***        cosh  :   Re(cosh(z)) = cosh(x) * cos(y)                    ****/
/* ***                  Im(cosh(z)) = sinh(x) * sin(y)                    ****/
/* ***                                                                    ****/
/* ***************************************************************************/


cinterval exp(const cinterval& z) throw()
{
  interval
    A( exp( Re(z) ) ),
    B(      Im(z)   );
  return cinterval( A*cos( B ) , A*sin( B ) );
}

cinterval exp2(const cinterval& z) throw()
{
	return exp(z*Ln2_interval);
}

cinterval exp10(const cinterval& z) throw()
{
	return exp(z*Ln10_interval);
}

cinterval expm1(const cinterval& z) throw()
// exp(z) - 1;
// Calculates nearly optimal inclusions for not too wide intervals z.
// Blomquist, 03.12.2008;
{
	const interval cancl_test = interval(0.9,1.1);
	interval rez(Re(z)), imz(Im(z));
	interval exp_x, sin_y, cos_y, h, Xt;
	cinterval res;
	
	exp_x = exp(rez);
	sin_y = sin(imz);
	cos_y = cos(imz);
	
	h = exp_x*cos_y;

	if (h < cancl_test && cos_y < cancl_test)
	{
		h = lnp1(-sqr(sin_y));
		times2pown(h,-1);
		// h = 0.5 * ln(1-sqr( sin(y) ));
		h = expm1(rez+h); // Cancellation also possible here!
	}
	else h = h - 1; // Cancellation possible here (real part)!
	// h: Real part;
	imz = exp_x * sin_y;
	res = cinterval(h,imz);
	return res;
}

cinterval cos(const cinterval& z) throw()
{
  interval
    A( Re(z) ),
    B( Im(z) );
  return cinterval( cos( A )*cosh( B ) , -sin( A )*sinh( B ) );
}

cinterval sin(const cinterval& z) throw()
{
  interval
    A( Re(z) ),
    B( Im(z) );
  return cinterval( sin( A )*cosh( B ) , cos( A )*sinh( B ) );
}

cinterval cosh(const cinterval& z) throw()
{
  interval
    A( Re(z) ),
    B( Im(z) );
  return cinterval( cos( B )*cosh( A ) , sin( B )*sinh( A ) );
}

cinterval sinh(const cinterval& z) throw()
{
  interval
    A( Re(z) ),
    B( Im(z) );
  return cinterval( cos( B )*sinh( A ) , sin( B )*cosh( A ) );
}



//-----------------------------------------------------------------------------
//
//  Section 2: tan, cot, tanh, coth
//
//  The implementation of cot, tanh, and coth is based on tan:
//
//  cot( z )  = tan( pi/2 - z )
//  tanh( z ) = transp( i * tan( transp( i * z ) )
//  coth( z ) = i * cot( i * z ) = i * tan( pi/2 - i * z )
//
//-----------------------------------------------------------------------------

//-- tan ------------------------------------------------------------ 040827 --
//
//  Complex tangent function
//
void horizontal_check( //-------------------------------------------- 040726 --
     const interval& hy, const interval& cosh_2y, real irez, real srez,
     const interval& hxl, const interval& hxu, real& resxl, real& resxu )
//
//  Subroutine of tangent function.
//  Check intersections with extremal curves of tan on a horizontal boundary.
//  This subroutine is only called if an intersection occurs.
//  In this case, sinh( 2 * hy ) <> 0.0 (poles are handled before).
//
//  There may be 1 or 2 intersections.
//  If intersections lie next to a boundary of rez, then it is impossible to
//  decide if there are 1 or 2 intersections.
//  In this case, 2 intersections are assumed
//  (valid enclosure, at the expense of a potential slight overestimation).
//
{
  bool both = false, left = false, right = false;

  if (srez - irez > Inf( PI() ))
    //  2 intersections
    both = true;
  else
    {
      interval
	res_l = cos( 2 * hxl ) - cosh_2y,
	res_u = cos( 2 * hxu ) - cosh_2y;

      if( Inf( res_l * res_u ) > 0.0 )
	//  2 intersections
	both = true;
      else
	{
	  if( Sup( res_l * res_u ) < 0.0 )
	    {
	      //  1 intersection (3 intersections are PI() apart)
	      //  neither of the intervals res_l and res_u contains zero
	      if( Inf( res_l ) > 0.0 )
		//  "left" intersection
		left = true;
	      else
		//  "right" intersection
		right = true;
	    }
	  else
	    //
	    //  1 (or both) intersections lie next to a boundary point
	    //  here, the task is to decide if 2 intersections occurs
	    //  if only 1 intersection takes place, then which one?
	    //
	    {
	      interval
		sin_2xl = sin( 2 * hxl ),
		sin_2xu = sin( 2 * hxu );

	      if( !disjoint( ZERO_INTERVAL(), res_l ) )
		//  intersection on the left boundary
		{
		  if( Inf( sin_2xl ) >= 0.0 )
		    // "left" intersection
		    {
		      left = true;
		      //  remove the intersection by changing res_l!
		      res_l = - ONE_INTERVAL();
		    }
		  else
		    {
		      if( Sup( sin_2xl ) <= 0.0 )
			// "right" intersection
			{
			  right = true;
			  //  remove the intersection by changing res_l!
			  res_l =  ONE_INTERVAL();
			}
		      else
			//  zero is interior point of sin_2xl
			//  if the real sine function has optimal precision,
			//  this case should never happen
			both = true;
		    }
		}

	      if( !disjoint( ZERO_INTERVAL(), res_u ) )
		//  intersection on the right boundary
		{
		  if( Inf( sin_2xu ) >= 0.0 )
		    // "left" intersection
		    {
		      left = true;
		      //  remove the intersection by changing res_u!
		      res_u = ONE_INTERVAL();
		    }
		  else
		    {
		      if( Sup( sin_2xu ) <= 0.0 )
			// "right" intersection
			{
			  right = true;
			  //  remove the intersection by changing res_u!
			  res_u = - ONE_INTERVAL();
			}
		      else
			//  zero is interior point of sin_2xu
			//  if the real sine function has optimal precision,
			//  this case should never happen
			both = true;
		    }
		}
	      //
	      //  finally, check if there is a second intersection
	      //
	      if( Inf( res_l * res_u ) < 0.0 )
		both = true;
	    }
	}
    }
  //
  //  Calculate extremal values
  //
  interval re_tan = 1 / sinh( 2 * abs( hy ) );

  //  "left" intersection, sin( 2x ) > 0
  if( left || both )
    {
      resxl = min( resxl, Inf( re_tan ) );
      resxu = max( resxu, Sup( re_tan ) );
    }

  //  "right" intersection, sin( 2x ) < 0
  if( right || both )
    {
      resxl = min( resxl, - Sup( re_tan ) );
      resxu = max( resxu, - Inf( re_tan ) );
    }
} // end horizontal_check


cinterval tan( const cinterval& z ) throw() //----------------------- 040827 --
{
  interval
      rez = Re(z),   // rez = z.re(),
      imz = Im(z);   // imz = z.im();

  real
    irez = Inf(rez),
    srez = Sup(rez),
    iimz = Inf(imz),
    simz = Sup(imz);

  interval
    hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

  real
    resxl, resxu, resyl, resyu;
  //
  //  1st: check for poles
  //
  if( ( !disjoint( ZERO_INTERVAL(), imz ) ) && ( !disjoint( ZERO_INTERVAL(), cos( rez ) ) ) )
    cxscthrow (STD_FKT_OUT_OF_DEF("cinterval tan( const cinterval& Z); Pole(s) in Z"));
  //
  //  2nd: real part
  //
  //  evaluate tan on vertical boundaries
  //
  interval
    cos_2rez   = cos( 2 * rez ),
    sinh_imz_2 = sqr( sinh( imz ) );

  interval
    re_tan_l = sin( 2 * hxl ) / ( 2 * ( sqr( cos( hxl ) ) + sinh_imz_2 ) ),
    re_tan_u = sin( 2 * hxu ) / ( 2 * ( sqr( cos( hxu ) ) + sinh_imz_2 ) );

  resxl = min( Inf( re_tan_l ), Inf( re_tan_u ) );
  resxu = max( Sup( re_tan_l ), Sup( re_tan_u ) );
  //
  //  look for extremal values on horizontal boundaries
  //  if one of the horizontal boundaries is the x-axes,
  //  then the complex tangent is the real tangent on this
  //  boundary, and due to monotonicity, its range
  //  is already included in the ranges of the vertical boundaries
  //
  if( irez < srez )
    {
      interval
	cosh_2yl = - 1 / cosh( 2 * hyl ),
	cosh_2yu = - 1 / cosh( 2 * hyu );

      if( !disjoint( cos_2rez, cosh_2yl ) && iimz != 0.0 )
	//extremal curve intersects lower boundary
	horizontal_check( hyl, cosh_2yl, irez, srez, hxl, hxu, resxl, resxu );

      if( !disjoint( cos_2rez, cosh_2yu ) && simz != 0.0 )
	//extremal curve intersects upper boundary
	horizontal_check( hyu, cosh_2yu, irez, srez, hxl, hxu, resxl, resxu );
    }
  //
  //  3rd: imaginary part
  //
  //  evaluate tan on horizontal boundaries
  //
  interval
    cos_rez_2 = sqr( cos( rez ) );

  interval
    im_tan_l = sinh( 2 * hyl ) / ( 2 * ( cos_rez_2 + sqr( sinh( hyl ) ) ) ),
    im_tan_u = sinh( 2 * hyu ) / ( 2 * ( cos_rez_2 + sqr( sinh( hyu ) ) ) );

  resyl = min( Inf( im_tan_l ), Inf( im_tan_u ) );
  resyu = max( Sup( im_tan_l ), Sup( im_tan_u ) );

  //
  //  look for extremal values on vertical boundaries
  //  here, the situation is simpler than for the real part
  //  if 2 intersections with extremal curves appear ,
  //  one lies in the lower half plane, the other in the upper half plane
  //
  interval
    cos_2xl = cos( 2 * hxl ),
    cos_2xu = cos( 2 * hxu );
  interval im_tan;

  if( iimz < 0.0 )
    //  intersection in lower half plane?
    {
      interval
	imz_h = interval( iimz, min( simz, 0.0 ) ),
	cosh_2imz = - 1 / cosh( 2 * imz_h );

      if( ( !disjoint( cosh_2imz, cos_2xl ) ) )
	//extremal curve intersects left boundary
	//in this case, sin( 2 * xl ) <> 0.0 (no poles here!)
	{
	  im_tan = - 1 / abs( sin( 2 * hxl ) );
	  resyl = min( resyl, Inf( im_tan ) );
	  resyu = max( resyu, Sup( im_tan ) );
	}
      if( ( !disjoint( cosh_2imz, cos_2xu ) ) )
	//extremal curve intersects right boundary
	//in this case, sin( 2 * xu ) <> 0.0 (no poles here!)
	{
	  im_tan = - 1 / abs( sin( 2 * hxu ) );
	  resyl = min( resyl, Inf( im_tan ) );
	  resyu = max( resyu, Sup( im_tan ) );
	}
    }
  if( simz > 0.0 )
    //  intersection in upper half plane?
    {
      interval
	imz_h = interval( max( iimz, 0.0 ), simz ),
	cosh_2imz = - 1 / cosh( 2 * imz_h );

      if( ( !disjoint( cosh_2imz, cos_2xl ) ) )
	//extremal curve intersects left boundary
	//in this case, sin( 2 * xl ) <> 0.0 (no poles here!)
	{
	  im_tan = + 1 / abs( sin( 2 * hxl ) );
	  resyl = min( resyl, Inf( im_tan ) );
	  resyu = max( resyu, Sup( im_tan ) );
	}
      if( ( !disjoint( cosh_2imz, cos_2xu ) ) )
	//extremal curve intersects right boundary
	//in this case, sin( 2 * xu ) <> 0.0 (no poles here!)
	{
	  im_tan = + 1 / abs( sin( 2 * hxu ) );
	  resyl = min( resyl, Inf( im_tan ) );
	  resyu = max( resyu, Sup( im_tan ) );
	}
    }

  return cinterval( interval( resxl, resxu ), interval( resyl, resyu ) );
  // return cinterval( interval(  resxu-resxl, resxu-resxl ), interval( resyu, resyu ) );
  //return cinterval( im_tan_l, im_tan_u );

} // end tan


//-- cot ------------------------------------------------------------ 040730 --
//
//  cot( z )  = tan( pi/2 - z )
//

cinterval cot( const cinterval& z ) throw()
{ // Improved cot-function; Blomquist,05.07.2005;
  // Improvements near the zeros z_(s,k) = Pi*(k+0.5), k=-4,-3,-2,-1,0,1,2,3,4;
  // where X = Re(z) must be a point interval near the zeros z_(s,k);
  //  Blomquist,08.07.2005;
    interval Rez(Re(z)),re;
    real irez = Inf(Rez); // irez = Inf( Re(z) )
    irez = irez/Inf_Pi - 0.5;
    double dbr = _double(irez);
    int p,s(sign(irez));
    p = s>=0 ? CUTINT(dbr+0.5) : CUTINT(dbr-0.5);
    if (-5<p && p<5) {
	re = (Inf_Pi)*interval(p+0.5) - Rez;
	re = diam(re)==0 ? re + (Pi2)*(p+0.5) : HALFPI() - Rez;
    }
    else re = HALFPI() - Rez;
    cinterval res(cinterval(re,-Im(z)));
    return tan(res);
}

//
//-- end cot ------------------------------------------------------------------

//-- tanh ----------------------------------------------------------- 040904 --
//
//  tanh( z ) = transp( i * tan( transp( i * z ) )
//
cinterval tanh( const cinterval& z ) throw()
{
  cinterval res = tan( cinterval( Im(z), Re(z) ) );
  return cinterval( Im(res), Re(res) );
}
//
//-- end tanh -----------------------------------------------------------------

//-- coth ----------------------------------------------------------- 040904 --
//
//   coth( z ) = i * cot( i * z );
//

cinterval coth( const cinterval& z ) throw()
{ // coth( z ) = i * cot( i * z );
    cinterval zh = cinterval( -Im(z), Re(z) ); //  zh = i*z;
    cinterval res = cot(zh);
    return( cinterval( -Im(res),Re(res) ) );
}
//
//-- end coth -----------------------------------------------------------------


//-----------------------------------------------------------------------------
//
//  Part II: Multi-valued functions
//
//-----------------------------------------------------------------------------
//
//
interval Atan(const interval& y, const interval& x) throw()
// Calculating an inclusion of atan(y/x) with x<>[0,0];
// This help function must only be used for POINT intervals y,x !!
// This function avoids internal overflow by calculating y/x. 
{
    const int c = 54;
    real Infx(Inf(x)),
	 Infy(Inf(y));
    int ex_x(expo(Infx)),
	ex_y(expo(Infy)),
	signx(sign(Infx)),
	signy(sign(Infy)),
	signq;
    interval res(0);

    if (signy!=0) {
    signq = signx * signy;
    if (ex_y-ex_x > c) res = signq>0 ? HALFPI() : -HALFPI();
    else res = atan(y/x);
    }

    return res;
}

interval Atan(const interval& y, const real& x) throw()
// Calculating an inclusion of atan(y/x) with x<>0.0;
// This help function must only be used for POINT intervals y !!
// This function avoids internal overflow by calculating y/x. 
{
    interval xi(x);
    return Atan(y,xi);
}

//
//  For the Multi-valued functions, two different procedures are 
//  implemented.
//
//  First, there is a subroutine for computing the principal value
//  of the respective function. The principal value is usually
//  analytic in a subset of the complex plane and undefined otherwise.
//
//  Second, there are procedures for computing interval supersets
//  of all function values of the respective function, where feasible.
//
//-----------------------------------------------------------------------------
//
//  Section 1: Argument functions
//
//-----------------------------------------------------------------------------

//-- Arg: analytic argument function -------------------------------- 040916 --
//
//  (i)   Arg(Z) subset (-pi,pi).
//  (ii)  Arg([0,0]) = 0.
//  (iii) Arg(Z) is undefined if Z contains negative real numbers.
//  (iv)  Otherwise, Arg(Z) is the interval hull of { Arg(z) | z in Z, z<>0 }.
//
//  atan is the inverse function of tan(t), t in (-pi/2,pi/2).
//
interval Arg( const cinterval& z ) throw()
{
  real 
    srez = Sup( Re(z) ),
    irez = Inf( Re(z) ), 
    simz = Sup( Im(z) ),
    iimz = Inf( Im(z) );

  interval
    hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

  real resl, resu;

  if( iimz > 0.0 )
    //  case I: Im(z) > 0
    {
      resl = ( srez > 0.0 ? Inf( Atan( hyl,hxu ) ) : ( srez < 0.0 ? Inf( Atan( hyu,hxu ) + PI() ) : Inf( HALFPI() ) ) );
      resu = ( irez > 0.0 ? Sup( Atan( hyu,hxl ) ) : ( irez < 0.0 ? Sup( Atan( hyl,hxl ) + PI() ) : Sup( HALFPI() ) ) );
      return interval( resl, resu );
    }
  else
    {
      if( simz < 0.0 )
	//  case II: Im(z) < 0
	{
	  resl = ( irez < 0.0 ? Inf( Atan( hyu,hxl ) - PI() ) : ( irez > 0.0 ? Inf( Atan( hyl,hxl ) ) : - Sup( HALFPI() ) ) );
	  resu = ( srez < 0.0 ? Sup( Atan( hyl,hxu ) - PI() ) : ( srez > 0.0 ? Sup( Atan( hyu,hxu ) ) : - Inf( HALFPI() ) ) );
	  return interval( resl, resu );
	}
      else
	// 0 in Im(z)
	{
	  if( irez > 0.0 )
	    //  case III: Re(z) > 0
	    //  z contains positive real values
	    {
	      resl = ( iimz < 0.0 ? Inf( Atan( hyl,hxl ) ) : 0.0 );
	      return interval( resl, Sup( Atan( hyu,hxl ) ) );
	    }
	  else
	    //  z contains nonpositive real numbers
	    {
	      if( irez < 0.0 )
	        {
		  //  case IV: z contains negative real numbers
                  cxscthrow (STD_FKT_OUT_OF_DEF("interval Arg( const cinterval& z ); z contains negative real numbers"));
		  return interval(0.0);
		}
	      else
		//  case V: 0 in z, but z doesn't contain negative real numbers
		{
		  if( srez > 0.0 )
		    //  diam( Re(z) > 0.0 )
		    {
		      resl = ( iimz < 0.0 ? - Sup( HALFPI() ) : 0.0 );
		      resu = ( simz > 0.0 ? Sup( HALFPI() ) : 0.0 );
		      return interval( resl, resu );
		    }
		  else
		    //  Re(z) == 0.0
		    {
		      if( iimz == 0.0 && simz == 0.0 )
			//  Z == 0
			return ZERO_INTERVAL();
		      else
			{
			  resl = ( iimz < 0.0 ? - Sup( HALFPI() ) : Inf( HALFPI() ) );
			  resu = ( simz > 0.0 ? Sup( HALFPI() ) : - Inf( HALFPI() ) );
			  return interval( resl, resu );
			}
		    }
		}
	    }
	}
    }
}
//
//-- end Arg ------------------------------------------------------------------

//-- arg: non-analytic argument function ---------------------------- 040916 --
//
//  (i)   arg(Z) is defined for all Z in IC.
//  (ii)  arg(Z) subset [-pi,3*pi/2].
//  (iii) arg(Z) == Arg(Z) if Arg(Z) is well-defined.
//
//  atan is the inverse function of tan(t), t in (-pi/2,pi/2).
//
interval arg( const cinterval& z ) throw()
{
  real
    srez = Sup( Re(z) ),
    irez = Inf( Re(z) ),
    simz = Sup( Im(z) ),
    iimz = Inf( Im(z) );

  real resl, resu;

  if( irez < 0.0 && iimz <= 0.0 && simz >= 0.0  )
    //  z contains negative real values
    {
      if( srez > 0.0 )
	//  0 in z and 0 interior point of Re(z)
	{
	  resl = ( iimz < 0.0 ? - Sup( PI() ) : 0.0 );
	  resu = ( ( iimz < 0.0 && simz == 0.0 ) ? 0.0 : Sup( PI() ) );
	  return interval( resl, resu );
	}
      else
      { // srez <= 0.0
	  if( iimz == simz )
	    //  z is real interval containing no positive values
	    return PI();
	  else
	    // sup( Re(z) ) <= 0, diam( Im(z) ) > 0
	    {
	      if( srez == 0.0 )
		{
		  resl = ( simz > 0.0 ? Inf( HALFPI() ) : - Sup( PI() ) );
		  resu = ( iimz < 0.0 ? ( simz > 0.0 ? Sup( 3 * HALFPI() ) : - Inf( HALFPI() ) ) : Sup( PI() ) );
		  return interval( resl, resu );
		}
	      else
		//   sup( Re(z) ) < 0, diam( Im(z) ) > 0
		{
		  interval hyl(iimz), hyu(simz);
		  resl = ( simz > 0.0 ? Inf( Atan( hyu,srez ) + PI() ) : - Sup( PI() ) );
		  resu = ( iimz < 0.0 ? ( simz > 0.0 ? Sup( Atan( hyl,srez ) + PI() ) : 
                                                       Sup( Atan( hyl,srez ) - PI() ) ) : Sup( PI() ) );
		  return interval( resl, resu );
		}
	    }
	}
    }
  else
    //  Arg(z) is well-defined
    return  Arg( z );
}
//
//-- end arg ------------------------------------------------------------------


/* ***************************************************************************/
/* ***************************************************************************/
/* ***                      Multi-valued functions                        ****/
/* ***************************************************************************/
/* ***************************************************************************/


//-- arg_inclmon: non-analytic inclusion-monotone argument function - 040617 --
//
//  (i)  arg_inclmon(Z) is defined for all Z in IC.
//  (ii) arg_inclmon(Z) = [-pi,pi] if Arg(Z) is not defined.
//
interval arg_inclmon( const cinterval& z ) throw()
{
  if( Inf( Re(z) ) < 0.0 && Inf( Im(z) ) <= 0.0 && Sup( Im(z) ) >= 0.0  )
    return  interval( -Sup( PI() ),Sup( PI() ) );
  else 
    return Arg(z);
}
//
//-- end arg_inclmon ----------------------------------------------------------


//-----------------------------------------------------------------------------
//
//  Section 2: Logarithms
//
//-----------------------------------------------------------------------------

//-- Ln: analytic natural logarithm --------------------------------- 040625 --
//
//  Ln(z) is undefined if z contains zero; z must not touch the negative real
//  axis from below;
//
cinterval Ln( const cinterval& z ) throw()
{ // Blomquist;
  real
    srez = Sup( Re(z) ),
    simz = Sup( Im(z) ),
    iimz = Inf( Im(z) );
  interval a1( abs(Re(z)) ),
           a2( abs(Im(z)) );
  if ( Inf(a1) == 0.0 && Inf(a2) == 0.0 )
      //  z contains 0
      cxscthrow(STD_FKT_OUT_OF_DEF("cinterval LN( const cinterval& z ); z contains 0"));
  if ( srez<0 && iimz<0 && simz>=0 ) 
      cxscthrow(STD_FKT_OUT_OF_DEF("cinterval LN( const cinterval& z ); z not allowed")); 
  return cinterval( ln_sqrtx2y2(Re(z),Im(z)),arg(z) );
}
//
//-- end Ln -------------------------------------------------------------------

//-- ln: non-analytic natural logarithm ----------------------------- 040923 --
//
//  ln(z) is undefined if z contains zero.
//
cinterval ln( const cinterval& z ) throw()
{
  interval a1( abs(Re(z)) ),
           a2( abs(Im(z)) );
  if ( Inf(a1) == 0.0 && Inf(a2) == 0.0 ) {
    //  z contains 0
    cxscthrow(STD_FKT_OUT_OF_DEF("cinterval ln( const cinterval& z ); z contains 0"));
    return z;
  }
  else
//    return cinterval( ln( abs_z_2 ) / 2.0 , arg( z ) );
    return cinterval( ln_sqrtx2y2(Re(z),Im(z)), arg(z) );
}
//
//-- end ln -------------------------------------------------------------------

cinterval lnp1(const cinterval& z) throw()
{  // ln(1+z);
	// Calculates nearly optimal inclusions for not too wide intervals z.
   // Blomquist, 03.12.2008;
	const real c1 = 1.0;
	cinterval y;
	interval abs_z(abs(z));
	real
		srez = Sup( Re(z) ),
		simz = Sup( Im(z) ),
		iimz = Inf( Im(z) );

	if (-1 <= z) //  z contains -1
		cxscthrow(STD_FKT_OUT_OF_DEF(
			"cinterval lnp1(const cinterval& z); z contains -1"));
	if ( srez<-1 && iimz<0 && simz>=0 ) 
		cxscthrow(STD_FKT_OUT_OF_DEF(
			"cinterval lnp1(const cinterval& z); z not allowed"));
	
	if (Sup(abs_z) < c1)
	{
		abs_z = Re(z);
		abs_z = lnp1( abs_z*(2+abs_z) + sqr(Im(z)) );
		times2pown(abs_z,-1);
		y = cinterval( abs_z, arg(1+z) );
	}
	else
		y = Ln(1+z);
	return y;
}

cinterval log2( const cinterval& z ) throw()
{
	return Ln(z) / Ln2_interval;
}

cinterval log10( const cinterval& z ) throw()
{
	return Ln(z) / Ln10_interval;
}

//-----------------------------------------------------------------------------
//
//  Section 3: Root functions
//
//-----------------------------------------------------------------------------

interval Sqrt_zpx( const interval& x, const interval& y )
// Calculating sqrt(|z|+|x|) without any internal overflow;
// Notice:  |z| = sqrt(x^2+y^2);   sqrt(|z|+|x|) < Maxreal 
{
    const int c1 = 1021;
    real Infx(Inf(x)), Infy(Inf(y));
    int ex_x(expo(Infx)), ex_y(expo(Infy));
    interval xc(abs(x)),yc(y),res;
    bool yeq0(Infy==0);

    if ((ex_x>=c1) || (ex_y>=c1)) 
    {  // To avoid overflow, scaling is necessary: 
	times2pown(xc,-2);
	if (yeq0) {
	    times2pown(xc,1); 
	    res = sqrt(xc);
	} else {
	times2pown(yc,-2);
	res = sqrt( sqrtx2y2(xc,yc)+xc);
	}
	times2pown(res,1);  // Backscaling
    } else   // Normal calculation without scaling:
	if (yeq0) {
	    times2pown(xc,1); 
	    res = sqrt(xc);
	} else	res = sqrt( sqrtx2y2(xc,yc)+xc);
    return res;
}

//-- sqrt: analytic square root ------------------------------------- 040626 --
//
interval Re_Sqrt_point( const interval& rez, const interval& imz ) // 040626 --
//
//  Real part of analytic square root of A POINT INTERVAL ONLY.
//  Do not use this as a general function
//  - it's only a subroutine for sqrt and sqrt_all.
//  The calculation is void if (rez+I*imz) is not a complex number.
//
{
  real
    irez = Inf( rez ),
    iimz = Inf( imz );

  if( iimz == 0.0 )
  {
      if( irez >= 0.0 )
	return sqrt( rez );
      else
	return ZERO_INTERVAL();
  }
  else
  { // iimz <> 0 
      if (irez >= 0.0) 
	  return INV_SQRT_2() * Sqrt_zpx(rez,imz);
      else 
	  return INV_SQRT_2() * abs(iimz) / Sqrt_zpx(rez,imz);
  }
}

interval Im_Sqrt_point( const interval& rez, const interval& imz ) // 040626 --
//
//  Imaginary part of analytic square root of A POINT INTERVAL ONLY
//  Do not use this as a general function
//  - it's only a subroutine for sqrt and sqrt_all
//  The calculation is void if (rez+I*imz) is not a complex number
//
{
  real
    irez = Inf( rez ),
    iimz = Inf( imz );

  if( iimz == 0.0 )
    {
      if( irez >= 0.0 )
	return ZERO_INTERVAL();
      else
	return sqrt( -rez );
    }
   else
    {
      if( Inf( rez ) >= 0.0 )
	return ( iimz * INV_SQRT_2() ) / Sqrt_zpx(rez,imz);
      else
	{
	  if( iimz > 0.0 )
	    return INV_SQRT_2() * Sqrt_zpx(rez,imz);
	  else
	    return -INV_SQRT_2() * Sqrt_zpx(rez,imz);
	}
    }
}


cinterval sqrt( const cinterval& z ) throw() // -------------------------------
//
//  Analytic square root function with z in the principal branch.
//  The branch cut is the negative real axis. On the branch cut 
//  the function values are defined by comming from the II quadrant.
//  Blomquist, 23.06.2005;
//
{
  real
    irez = Inf( Re(z) ),
    srez = Sup( Re(z) ),
    iimz = Inf( Im(z) ),
    simz = Sup( Im(z) );
  interval
    hxl( irez ), hxu( srez ), hyl( iimz ), hyu( simz );
  real
    resxl,  resxu, resyl, resyu;

  if( irez < 0.0 && iimz < 0.0 && simz >= 0.0 ) 
  //  if z is not in the principal branch then the inclusion monotony is violated!
      cxscthrow(STD_FKT_OUT_OF_DEF("cinterval sqrt(const cinterval& z); z not in the principal branch."));

  if (iimz>=0.0)
  {
      resxl = Inf( Re_Sqrt_point( hxl,hyl ) );
      resxu = Sup( Re_Sqrt_point( hxu,hyu ) );

      resyl = Inf( Im_Sqrt_point( hxu,hyl ) );
      resyu = Sup( Im_Sqrt_point( hxl,hyu ) );
  } else 
      if (simz<=0.0) {
	  resxl = Inf( Re_Sqrt_point( hxl,hyu ) );
	  resxu = Sup( Re_Sqrt_point( hxu,hyl ) );

	  resyl = Inf( Im_Sqrt_point( hxl,hyl ) );
	  resyu = Sup( Im_Sqrt_point( hxu,hyu ) );
      } else {
	  resxl = Inf( sqrt( hxl ) );
	  resxu = ( -iimz > simz ? Sup( Re_Sqrt_point( hxu, hyl ) ) 
                                 : Sup( Re_Sqrt_point( hxu, hyu ) ) );
	  resyl = Inf( Im_Sqrt_point( hxl,hyl ) );
	  resyu = Sup( Im_Sqrt_point( hxl,hyu ) );
      }

  return cinterval( interval( resxl,resxu ), interval( resyl,resyu ) );
}

cinterval sqrtp1m1(const cinterval& z) throw()
// sqrt(1+z)-1;
// Calculates nearly optimal inclusions for not too wide intervals z.
// Blomquist, 08.07.2008;
{
	const real c = 1.5;
	cinterval res;
	interval absz(abs(z));
	real Sup_absz(Sup(absz));
	
	if (Sup_absz < c)
		res = z / (sqrt(z+1) + 1);
	else 
		res = sqrt(z+1) - 1;
	return res;
}

cinterval sqrt1px2(const cinterval& z) throw()
// sqrt(1+z^2);
// Blomquist, 03.07.2008;
{
	const real c = 5e8;
	cinterval res;
	interval absz(abs(z));
	real Inf_absz(Inf(absz));
	
	if (Inf_absz > c)
	{
		absz = 1 / interval(Inf_absz);
		Inf_absz = Sup(absz);
		res = cinterval( interval(-Inf_absz,Inf_absz),
							  interval(-Inf_absz,Inf_absz) );
		 // res is the correcture interval, i.e.
		 // z + res or -z + res is the inclusionof sqrt{1+z^2}
		res = (Inf(Re(z))>=0)? z + res : -z + res;
	}
	else 
	{
		res = cinterval( interval(0), interval(1) ); // res = i
		if ( Sup(abs(z-res))<0.5 || Sup(abs(z+res))<0.5 )
		{
			res = cinterval(-Im(z),Re(z)); // Res = i*z;
		   // (1 - i*z)*(1 + i*z) = 1+z^2;
			res = sqrt( (1-res)*(1+res) );
		}
		else
			res = sqrt(1+sqr(z));
	}
	if (Inf(Re(res))<0)
		res = cinterval( interval(0.0,Sup(Re(res))) , Im(res) );

	return res;
}
// -- end sqrt1px2 ------------------------------------------------------------

cinterval sqrtx2m1(const cinterval& z) throw()
// sqrt(z^2-1);
// Calculates nearly optimal inclusions for not too wide intervals z.
// Blomquist, 04.12.2008;
{
	const real c = 5e8;
	cinterval res,u;
	interval absz(abs(z));
	real Inf_absz(Inf(absz));
	
	if (Inf_absz > c)
	{
		absz = 1 / interval(Inf_absz);
		Inf_absz = Sup(absz);
		res = cinterval(interval(-Inf_absz,Inf_absz),
							 interval(-Inf_absz,Inf_absz)); // res = Delta
	   // res is the correcture interval, i.e.
		res = (Inf(Re(z))>=0)? z + res : -z + res;
	}
	else 
	{
		res = z-1;  u = z+1;
		res = (Sup(abs(res))<0.5 || Sup(abs(u))<0.5)? sqrt(res*u) : sqrt(sqr(z)-1);
	}
	
	if (Inf(Re(res))<0)
		res = cinterval( interval(0.0,Sup(Re(res))) , Im(res) );

	return res;
}

cinterval sqrt1mx2(const cinterval& z) throw()
// sqrt(1-z^2);
// Calculates nearly optimal inclusions for not too wide intervals z.
// Blomquist, 04.12.2008;
{
	const real c = 5e8;
	cinterval res,u;
	interval absz(abs(z));
	real Inf_absz(Inf(absz));
	
	if (Inf_absz > c)
	{
		absz = 1 / interval(Inf_absz);
		Inf_absz = Sup(absz);
		res = cinterval(interval(-Inf_absz,Inf_absz),
							 interval(-Inf_absz,Inf_absz)); // res = Delta
		u = cinterval(-Im(z),Re(z)); // u = i*z;
		// res is the correcture interval, i.e.
		// i*z + res or -i*z + res is the inclusion of sqrt{1-z^2}
		res = (Inf(Im(z))>=0)? -u + res : u + res;
	}
	else 
	{
		res = 1-z;  u = 1+z;
		res = (Sup(abs(res))<0.5 || Sup(abs(u))<0.5)? sqrt(res*u) : sqrt(1-sqr(z));
	}
	if (Inf(Re(res))<0)
		res = cinterval( interval(0.0,Sup(Re(res))) , Im(res) );

	return res;
}

//-- sqrt_all ------------------------------------------------------- 040621 --
//
//  sqrt_all(z) computes a list of 2 intervals containing all square roots of z
//
std::list<cinterval> sqrt_all( const cinterval& z )
{
  real
    irez = Inf( Re(z) ),
    srez = Sup( Re(z) ),
    iimz = Inf( Im(z) ),
    simz = Sup( Im(z) );
  interval
    hxl( irez ), hxu( srez ), hyl( iimz ), hyu( simz );
  real
    resxl,  resxu, resyl, resyu;
  cinterval w;

  if( irez < 0.0 && iimz <= 0.0 && simz >= 0.0 )
    //  z contains negative real values
    {
      if( iimz == 0.0 )
	//  z in upper half plane
	//  principal values can be used
	{
	  //  min( Re ( sqrt( z ) ) ) in lower left  corner
	  //  max( Re ( sqrt( z ) ) ) in upper right corner
	  resxl = Inf( Re_Sqrt_point( hxl, hyl ) );
	  resxu = Sup( Re_Sqrt_point( hxu, hyu ) );
	  //  min( Im ( sqrt( z ) ) ) in lower right corner
	  //  max( Im ( sqrt( z ) ) ) in upper left  corner
	  resyl = Inf( Im_Sqrt_point( hxu, hyl ) );
	  resyu = Sup( Im_Sqrt_point( hxl, hyu ) );
	}
      else
	{
	  if( simz == 0.0 )
	    //  z in lower half plane
	    //  principal values can be used in lower corners
	    {
	      //  z in lower half plane
	      //  min( Re ( sqrt( z ) ) ) in upper left  corner
	      //  max( Re ( sqrt( z ) ) ) in lower right corner
	      resxl = 0;                              ;
	      resxu = Sup( Re_Sqrt_point( hxu, hyl ) );
	      //  min( Im ( sqrt( z ) ) ) in lower left  corner
	      //  max( Im ( sqrt( z ) ) ) in upper right corner
	      resyl = Inf( Im_Sqrt_point( hxl, hyl ) );
	      resyu = ( srez > 0.0 ? 0.0 : -Inf( sqrt( -hxu ) ) );
	    }
	  else
	    //  0 is interior point of Im( z )
	    {
	      if( srez <= 0.0 )
		{
		  //  0 is no interior point of Re (z )
		  //  in quadrant III,    Re( z ) = Im_Sqrt_point(|x|,y),
		  //	                  Im( z ) = Re_Sqrt_point(|x|,y)
		  //  min( Re ( sqrt( z ) ) ) in lower right corner = quadrant II/III
		  //  max( Re ( sqrt( z ) ) ) in upper right corner = quadrant II
		  resxl = Inf( Im_Sqrt_point( -hxu, hyl ) );
		  resxu = Sup( Re_Sqrt_point( hxu, hyu ) );
		  //  min( Im ( sqrt( z ) ) ) on real line
		  //  max( Im ( sqrt( z ) ) ) in lower or upper left corner
		  resyl = Inf( sqrt( -hxu ) );
		  resyu = ( - iimz > simz ? Sup( Re_Sqrt_point( -hxl, hyl ) ) : Sup( Im_Sqrt_point( hxl, hyu ) ) );
		}
	      else
		//  0 is interior point of z
		//  here, the principal values apply in all corner points
		{
		  resxl = 0;
		  //  max( Re ( sqrt( z ) ) ) in lower or upper right corner
		  resxu = ( - iimz > simz ? Sup( Re_Sqrt_point( hxu, hyl ) ) : Sup( Re_Sqrt_point( hxu, hyu ) ) );
		  //  min( Im ( sqrt( z ) ) ) in lower left corner
		  //  max( Im ( sqrt( z ) ) ) in upper left corner
		  resyl = Inf( Im_Sqrt_point( hxl, hyl ) );
		  resyu = Sup( Im_Sqrt_point( hxl, hyu ) );
		}
	    }
	}
      w = cinterval( interval( resxl, resxu ), interval( resyl, resyu ) );
    }
  else
    //  sqrt( z ) is well-defined
    w = sqrt( z );

  std::list<cinterval> res;
  res.push_back( w );
  res.push_back( -w );

  return res;
}
//
//-- end sqrt_all -------------------------------------------------------------


//-- sqrt(z,n): analytic n-th root ---------------------------------- 040624 --
//
interval Re_Sqrt_point( const interval& rez, const interval& imz,
			int n ) // before: unsigned int n  ---------- 040624 --
//
//  Real part of analytic n-th root.
//  Do not use this as a general function
//  - it's only a subroutine for sqrt(z,n) and sqrt_all(z,n).
//  The calculation is validated but largely overestimating
//  if (rez+I*imz) is not a complex number.
//
{
  interval abs_z_2 = sqr( rez ) + sqr( imz );
  if( Sup( abs_z_2 ) == 0.0 )
    //  z == 0
    return ZERO_INTERVAL();
  else
    return sqrt( abs_z_2, 2 * n ) * cos( Arg( cinterval( rez, imz ) ) / n );
}

interval Im_Sqrt_point( const interval& rez, const interval& imz,
			int n ) // before: unsigned int n  --- 040624 --
//
//  Imaginary part of analytic n-th root.
//  Do not use this as a general function
//  - it's only a subroutine for sqrt(z,n) and sqrt_all(z,n).
//  The calculation is validated but largely overestimating
//  if (rez+I*imz) is not a complex number.
//
{
  interval abs_z_2 = sqr( rez ) + sqr( imz );
  if( Sup( abs_z_2 ) == 0.0 )
    //  z == 0
    return ZERO_INTERVAL();
  else
    return sqrt( abs_z_2, 2 * n ) * sin( Arg( cinterval( rez, imz ) ) / n );
}

cinterval sqrt( const cinterval& z, int n ) throw() // ----- 040915 --
//
//  Analytic n-th root function
//  sqrt(z,n) is undefined if z contains negative real numbers.
//
{
  if( n == 0 )
    return cinterval( ONE_INTERVAL() );
  if( n == 1 )
    return z;
  if( n == 2 )
    return sqrt( z );
  else
    {
      real
	irez = Inf( Re(z) ),
	srez = Sup( Re(z) ),
	iimz = Inf( Im(z) ),
	simz = Sup( Im(z) );
      interval
	hxl( irez ), hxu( srez ), hyl( iimz ), hyu( simz );
      real
	resxl,  resxu, resyl, resyu;

      if( irez < 0.0 && iimz <= 0.0 && simz >= 0.0 )
        {
	  //  z contains negative real values
	  cxscthrow(STD_FKT_OUT_OF_DEF("cinterval sqrt(const cinterval& z, int n ); z contains negative real values."));
          return z;
        }
      else
	{
	  if( simz < 0.0 )
	    {
	      //  z in lower half plane
              cinterval hres = sqrt( cinterval( Re(z), -Im(z) ), n );
	      return cinterval( Re(hres), -Im(hres) );
	    }
	  else
	    {
	      if( iimz > 0.0 ) 
		{
		  //  z in upper half plane
		  interval tangle = tan( ( PI() * n ) / ( 2 * ( n-1 ) ) );
		  real tanglel = Inf( tangle ), 
                       tangleu = Sup( tangle );
		  //
		  //  min( Re( Root( z ) ) )
		  //
		  if ( irez >= 0.0  ||  Sup( hyl / irez ) <= tanglel )
		    //  lower boundary right of phi = n*Pi/(2n-2)
		    //  min( Re( Root( z ) ) ) in lower left corner
		    resxl = Inf( Re_Sqrt_point( hxl, hyl, n ) );
		  else
		    {
		      if( srez < 0.0 && Inf( hyl / srez ) >= tangleu )
			//  lower boundary left of phi = n*Pi/(2n-2)
			//  min( Re( Root( z ) ) ) in lower right corner
			resxl = Inf( Re_Sqrt_point( hxu, hyl, n ) );
		      else
			//  lower boundary intersects phi = n*Pi/(2n-2)
			//  min( Re( Root( z ) ) ) on intersection
			resxl = Inf( Re_Sqrt_point( iimz / tangle , hyl, n ) );
		    }
		  //
		  //  max( Re( Root( z ) ) )
		  //
		  if ( irez >= 0.0 || Sup( hyu / irez ) <= tanglel )
		    //  upper boundary right of phi = n*Pi/(2n-2)
		    //  max( Re( Root( z ) ) ) in upper right corner
		    resxu = Sup( Re_Sqrt_point( interval(srez), interval(simz), n ) );
		  else
		    {
		      if ( srez < 0.0 && Inf( hyu / srez ) >= tangleu )
			//  upper boundary left of phi = n*Pi/(2n-2)
			//  max( Re( Root( z ) ) ) in upper left corner
			resxu = Sup( Re_Sqrt_point( hxl, hyu, n ) );
		      else
			//  upper boundary intersects phi = n*Pi/(2n-2)
			//  max( Re( Root( z ) ) ) on upper left or right corner
			resxu = max( Sup( Re_Sqrt_point( hxl, hyu, n ) ), Sup( Re_Sqrt_point( hxu, hyu, n ) ) );
		    }
		  //
		  //  min( Im( Root( z ) ) )
		  //
		  if ( srez >= 0.0  || Sup( hyl / srez ) <= tanglel )
		    //  right boundary right of phi = n*Pi/(2n-2)
		    //  min( Im( Root( z ) ) ) in lower right corner
		    resyl = Inf( Im_Sqrt_point( hxu, hyl, n ) );
		  else
		    {
		      if( Inf( hyu / srez ) >= tangleu )
			//  right boundary left of phi = n*Pi/(2n-2)
			//  min( Im( Root( z ) ) ) in upper right corner
			resyl = Inf( Im_Sqrt_point( hxu, hyu, n ) );
		      else
			//  right boundary intersects phi = n*Pi/(2n-2)
			//  min( Im( Root( z ) ) ) on intersection
			resyl = Inf( Im_Sqrt_point( hxu, srez * tangle, n ) );
		    }
		  //
		  //  max( Im( Root( z ) ) )
		  //
		  if( irez >= 0.0 || Sup( hyl / irez ) <= tanglel )
		    //  left boundary right of phi = n*Pi/(2n-2)
		    //  max( Im( Root( z ) ) ) in upper left corner
		    resyu = Sup( Im_Sqrt_point( hxl, hyu, n ) );
		  else
		    {
		      if( Inf( hyu / irez ) >= tangleu )
			//  left boundary left of phi = n*Pi/(2n-2)
			//  max( Im( Root( z ) ) ) in lower left corner
			resyu = Sup( Im_Sqrt_point( hxl, hyl, n ) );
		      else
			//  left boundary intersects phi = n*Pi/(2n-2)
			//  max( Im( Root( z ) ) ) on lower or upper left corner
			resyu = max( Sup( Im_Sqrt_point( hxl, hyl, n ) ), Sup( Im_Sqrt_point( hxl, hyu, n ) ) );
		    }
		}
	      else
		{
		  //  z intersects positive real axis
		  //  min( Re( Root( z ) ) ) on real line
		  //  max( Re( Root( z ) ) ) in lower or upper right corner
		    resxl = ( irez == 0.0 ? 0.0 : Inf( sqrt( hxl, n ) ) );
		  resxu = ( - iimz > simz ? Sup( Re_Sqrt_point( hxu, hyl, n ) ) : Sup( Re_Sqrt_point( hxu, hyu, n ) ) );
		  //  min( Im ( sqrt( z ) ) ) in lower left corner
		  //  max( Im ( sqrt( z ) ) ) in upper left corner
		  resyl = Inf( Im_Sqrt_point( hxl, hyl, n ) );
		  resyu = Sup( Im_Sqrt_point( hxl, hyu, n ) );
		}
	      return cinterval( interval( resxl, resxu ), interval( resyl, resyu ) );
	    }
	}
    }
}
//
//-- end sqrt -----------------------------------------------------------------


//-- sqrt_all ------------------------------------------------------- 040628 --
//
std::list<cinterval> sqrt_all( const cinterval& z, int n )
//
//  sqrt_all(z,n) computes a list of n intervals containing all n-th roots of z
//
//  For n >=3, computing the optimal interval bounds is very expensive
//  and thus not considered cost-effective.
//
//  Hence, the polar form is used to calculate sqrt_all(z,n)
//  (observed overestimation of Re() and Im() in test cases: 5-15%).
//
//  z is enclosed into an interval in polar coordinates
//  (i.e. a sector of an annulus), from which the roots
//  (also polar intervals) are computed error-free
//  (save roundoff, which is enclosed).
//
//  The the final inclusion of the roots into a rectangular complex interval
//  involves some additional overestimation.
//
{
  std::list<cinterval> res;

  if( n == 0 )
    {
      res.push_back( cinterval( ONE_INTERVAL(), ZERO_INTERVAL() ) );
      return res;
    }
  else if( n == 1 )
    {
      res.push_back(z);
      return res;
    }
  else if( n == 2 ) return sqrt_all( z );
  else
    {
      interval
	arg_z = arg( z ), root_abs_z = sqrt( abs( z ), n );

      for(int k = 0; k < n; k++)
	{
	  interval arg_k = ( arg_z + 2 * k * PI() ) / n;

	  res.push_back( cinterval( root_abs_z * cos( arg_k ),
				    root_abs_z * sin( arg_k ) ) );
	}
      return res;
    }
}
//
//-- end sqrt_all -------------------------------------------------------------



//-----------------------------------------------------------------------------
//
//  Section 4: Power functions
//
//-----------------------------------------------------------------------------

//-- power_fast ----------------------------------------------------- 040714 --
//
//  Fast, validated power function for integer powers, based on exp and ln.
//  Medium amount of overestimation.
//
/*!
\param z The complex interval
\param n The integer exponent
\return The computed complex interval

Fast, validated power function for integer powers, based on \f$ \exp \f$ and \f$ \ln \f$. Medium amount of overestimation.

\sa power(const cinterval&,int)
*/
cinterval power_fast( const cinterval& z, int n ) throw()
{
  if( n == 0 )
    return cinterval( ONE_INTERVAL() );
  else if( n == 1 )
    return z;
  else if( n == -1 )
    return 1 / z;
  else if( n == 2 )
    return sqr(z); 
  else
    //  n >= 3  or  n <= -2
    //  If z is a large interval, z^n is a distorted set, for which the
    //  inclusion into a complex rectangle contains large overestimation.
    //  For example, if  n * Arg( z )  intersects the cartesian axes
    //  more than once then  0  is contained in the rectangular enclosure
    //  of z^n.
    //  For n <= -2, also inversion of z or z^n is required;
    //  both inversions lead to large overestimation of the resulting interval.
    //
    //  The computation of an optimal rectangular enclosure is implemented
    //  in power(z,n). In power_fast(z,n), z is enclosed into a sector S of
    //  an annulus. S^n is again some sector S' of a (different) annulus.
    //  S' is calculated exactly (apart from roundoff), and then enclosed
    //  into a rectangle. There is a certain amount of overestimation
    //  compared with the optimal rectangular enclosure of z^n, but the
    //  calculations are much cheaper here.
    {
      interval abs_z = abs( z );

      if( n < 0 && Inf( abs_z ) == 0.0 )
	//  z contains 0
        cxscthrow (STD_FKT_OUT_OF_DEF("cinterval power_fast(const cinterval& z, int n ); z contains 0."));
      if( Sup( abs_z ) == 0.0 )
        return cinterval( ZERO_INTERVAL(), ZERO_INTERVAL() );
      else
        {
          interval arg_z = arg( z );
	  interval abs_z_n = exp( n * ln( abs_z ) );

	  return cinterval( abs_z_n * cos( n * arg_z ),
                            abs_z_n * sin( n * arg_z ) );
	}
    }
}
//
//-- end power_fast -----------------------------------------------------------


//-- power ---------------------------------------------------------- 040720 --
//
//  Power function for integer powers with optimal (save roundoff) accuracy.
//
//  The extremal curves of Re(z^n) and Im(z^n) are straight lines defined
//  by the rays  arg(z) = k pi/ ( 2n - 2 ), k = 0,...4n-5.
//
//  Intersections of these rays with the boundary of z are called
//  "ray intersections" in the following.
//
cinterval power_point( const cinterval& z, int n ) //---------------- 040715 --
//
//  z^n for A POINT INTERVAL z ONLY.
//  Do not use this as a general function
//  - it's only a subroutine for power.
//  The calculation may break down otherwise.
//  The case 0 in z for negative n is handled in power(z,n).
//
{
  if( Inf( Re(z) ) == 0.0 && Inf( Im(z) ) == 0.0 )
    //  if z is a valid point interval, it must be 0
    return cinterval( ZERO_INTERVAL(), ZERO_INTERVAL() );
  else
    {
      interval abs_z = abs( z );
      interval arg_z = arg( z );
      interval abs_z_n = exp( n * ln( abs_z ) );

      return cinterval( abs_z_n * cos( n * arg_z ),
			abs_z_n * sin( n * arg_z ) );
    }
}

void update_res( const cinterval& res, //---------------------------- 040716 --
                 real& resxl, real& resxu, real& resyl, real& resyu )
//  Subroutine of power(z,n).
//  Update extremal values of power function.
{
  resxl = min( resxl, Inf( Re(res) ) );
  resxu = max( resxu, Sup( Re(res) ) );
  resyl = min( resyl, Inf( Im(res) ) );
  resyu = max( resyu, Sup( Im(res) ) );
}

void horizontal_check( //-------------------------------------------- 040720 --
     const interval& hy, real hyy, const interval& arg_h,
     real irez, real srez,
     real& resxl, real& resxu, real& resyl, real& resyu, int n )
//  Subroutine of power(z,n).
//  Check all relevant ray intersections on a horizontal boundary.
{
  double r_int;
  int
    il1, il2, iu1, iu2;
  cinterval res;
  //
  //  Here, the most difficult task is
  //  to determine the relevant ray intersections
  //  Both arg_hl and n can be negative, which makes it very complicated
  //  to decide the correct indexes for the "rightmost" intersections, etc.
  //
  //  To simplify the distinction of cases, we introduce the variable
  //  nofrays (number of rays) = abs( n-1 )
  //
  //  Note that we still have to pass n to power_point
  //
  iu1 = n-1;  if (iu1<0) iu1 = -iu1;
  int nofrays = iu1;
  real arg_hlR = Inf( 2 * nofrays * arg_h / PI() );
  double arg_hl = _double(arg_hlR);
  if( arg_hl >= 0.0 )
    modf( arg_hl, &r_int );
  else
    modf( arg_hl-1, &r_int );
  il1 = int( r_int + 1.0 );
  real arg_huR = Sup( 2 * nofrays * arg_h / PI() );
  double arg_hu = _double(arg_huR);
  if( arg_hu >= 0.0 )
    modf( arg_hu, &r_int );
  else
    modf( arg_hu-1, &r_int );
  iu1 = int( r_int );

  if( iu1 >= il1 ) 
  {
   //  at least one ray intersection
   //  maybe more?
      if( iu1 > il1 + 3 )
	//
	//  we're in trouble:
	//  there are more than 4 ray intersections
	//  now we must decide which of these are relevant
	//  depending on il1, iu1, n and on the quadrants,
	//  4 to 7 intersections must be considered
	//
      { // Neu Blomquist 13.12.2010
	if( n > 0 )
	  //  outmost intersections are relevant
	{
	    if( irez >= 0.0 )
	      //  z in quadrants I or IV:
	      //  only 4 rightmost intersections are relevant
	      {
		if( hyy < 0.0 )
		  //  quadrant IV
		  il1 = iu1 - 3;
		else
		  //  quadrant I
		  iu1 = il1 + 3;
	      }
	    else if( srez <= 0.0 )
	      //  z in quadrants II or III:
	      //  only 4 leftmost intersections are relevant
	      {
		if( hyy > 0.0 )
		  //  quadrant II
		  il1 = iu1 - 3;
		else
		  //  quadrant III
		  iu1 = il1 + 3;
	      }
	    else
	      //  z intersects imaginary axes
	      //  we may have two lists of intersections!
	      {
		iu2 = iu1;
		il2 = iu2 - 3;
		iu1 = il1 + 3;
		//  remove intersection points that are contained in both lists
		if( il2 <= iu1 )
		  il2 = iu1 + 1;
		//
		//  here, list 2 is processed
		//  list 1 is processed below
		//
		for(int i = il2; i <= iu2; i++)
		  {
		    interval cotangle = cot( ( PI() * i ) / ( 2 * nofrays ) );
		    res = power_point( cinterval( hy * cotangle , hy ), n );
		    update_res( res, resxl, resxu, resyl, resyu );
		  }
	      }
	}
	else
	  //  n < 0:
	  //  innermost intersections are relevant
	{
	    if( irez >= 0.0 )
	      //  z in quadrants I or IV:
	      //  only 4 leftmost intersections are relevant
	      {
		if( hyy > 0.0 )
		  //  quadrant I
		  il1 = iu1 - 3;
		else
		  //  quadrant IV
		  iu1 = il1 + 3;
	      }
	    else if( srez <= 0.0 )
	      //  z in quadrants II or III:
	      //  only 4 rightmost intersections are relevant
	      {
		if( hyy < 0.0 )
		  //  quadrant III
		  il1 = iu1 - 3;
		else
		  //  quadrant II
		  iu1 = il1 + 3;
	      }
	    else
	      //  z intersects imaginary axes
	      //  we have one lists of 5 to 7 intersections
	      {
		if( hyy > 0.0 )
		  {
		    il2 = nofrays - 3;
		    iu2 = nofrays + 3;
		  }
		else
		  {
		    il2 = -nofrays - 3;
		    iu2 = -nofrays + 3;
		  }
		//  list 1 contains all intersections
		//  list 2 is a filter for all relevant intersections
		//  store intersection of lists 1 and 2 in list 1
		il1 = ( il1 > il2 ? il1 : il2 );
		iu1 = ( iu1 < iu2 ? iu1 : iu2 );
	      }
	}
      } // Neu Blomquist 13.12.2010
      //
      //  list 1 has been left for processing
      //
      for(int i = il1; i <= iu1; i++)
	{
	  interval cotangle = cot( ( PI() * i ) / ( 2 * nofrays ) );
	  res = power_point( cinterval( hy * cotangle , hy ), n );
	  update_res( res, resxl, resxu, resyl, resyu );
	}
  }  
}

void vertical_check( //---------------------------------------------- 040720 --
     const interval& hx, real hxx, const interval& arg_h,
     real iimz, real simz,
     real& resxl, real& resxu, real& resyl, real& resyu, int n )
//  Subroutine of power(z,n).
//  Check all relevant ray intersections on a vertical boundary.
{
  double r_int;
  int
    il1, il2, iu1, iu2;
  cinterval res;
  //
  //  Here, the most difficult task is
  //  to determine the relevant ray intersections
  //  Both arg_hl and n can be negative, which makes it very complicated
  //  to decide the correct indexes for the topmost intersections, etc.
  //
  //  To simplify the distinction of cases, we introduce the variable
  //  nofrays (number of rays) = abs( n-1 )
  //
  //  Note that we still have to pass n to power_point
  //
  iu1 = n-1;  if (iu1<0) iu1 = -iu1;
  int nofrays = iu1;

  real arg_hlR = Inf( 2 * nofrays * arg_h / PI() );
  double arg_hl = _double(arg_hlR);
  if( arg_hl >= 0.0 )
    modf( arg_hl, &r_int );
  else
    modf( arg_hl-1, &r_int );
  il1 = int( r_int + 1.0 );
  real arg_huR = Sup( 2 * nofrays * arg_h / PI() );
  double arg_hu = _double(arg_huR);
  if( arg_hu >= 0.0 )
    modf( arg_hu, &r_int );
  else
    modf( arg_hu-1, &r_int );
  iu1 = int( r_int );

  if( iu1 >= il1 ) 
  {
  //  at least one ray intersection
  //  maybe more?
      if( iu1 > il1 + 3 )
	//
	//  we're in trouble:
	//  there are more than 4 ray intersections
	//  now we must decide which of these are relevant
	//  depending on il1, iu1, n and on the quadrants,
	//  4 to 7 intersections must be considered
	//
      { // Neu Blomquist 13.12.2010
	if( n > 0 )
	  //  outmost intersections are relevant
	{
	    if( iimz >= 0.0 )
	      //  z in quadrants I or II:
	      //  only 4 topmost intersections are relevant
	      {
		if( hxx > 0.0 )
		  //  quadrant I
		  il1 = iu1 - 3;
		else
		  //  quadrant II
		  iu1 = il1 + 3;
	      }
	    else if( simz <= 0.0 )
	      //  z in quadrants III or IV:
	      //  only 4 lowest intersections are relevant
	      {
		if( hxx < 0.0 )
		  //  quadrant III
		  il1 = iu1 - 3;
		else
		  //  quadrant IV
		  iu1 = il1 + 3;
	      }
	    else
	      //  z intersects real axes
	      //  we may have two lists of intersections!
	      {
		iu2 = iu1;
		il2 = iu2 - 3;
		iu1 = il1 + 3;
		//  remove intersection points that are contained in both lists
		if( il2 <= iu1 )
		  il2 = iu1 + 1;
		//
		//  here, list 2 is processed
		//  list 1 is processed below
		//
		for(int i = il2; i <= iu2; i++)
		  {
		    interval tangle = tan( ( PI() * i ) / ( 2 * nofrays ) );
		    res = power_point( cinterval( hx, hx * tangle ), n );
		    update_res( res, resxl, resxu, resyl, resyu );
		  }
	      }
	}
	else
	  //  n < 0:
	  //  innermost intersections are relevant
	{
	    if( iimz >= 0.0 )
	      //  z in quadrants I or II:
	      //  only 4 lowest intersections are relevant
	      {
		if( hxx < 0.0 )
		  //  quadrant II
		  il1 = iu1 - 3;
		else
		  //  quadrant I
		  iu1 = il1 + 3;
	      }
	    else if( simz <= 0.0 )
	      //  z in quadrants III or IV:
	      //  only 4 topmost intersections are relevant
	      {
		if( hxx > 0.0 )
		  //  quadrant IV
		  il1 = iu1 - 3;
		else
		  //  quadrant III
		  iu1 = il1 + 3;
	      }
	    else
	      //  z intersects real axes
	      //  we have one lists of 5 to 7 intersections
	      {
		if( hxx > 0.0 )
		  {
		    il2 = -3;
		    iu2 =  3;
		  }
		else
		  {
		    il2 = 2 * nofrays - 3;
		    iu2 = 2 * nofrays + 3;
		  }
		//  list 1 contains all intersections
		//  list 2 is a filter for all relevant intersections
		//  store intersection of lists 1 and 2 in list 1
		il1 = ( il1 > il2 ? il1 : il2 );
		iu1 = ( iu1 < iu2 ? iu1 : iu2 );
	      }
	}
      } // Neu Blomquist 13.12.2010 
      //
      //  list 1 has been left for processing
      //
      for(int i = il1; i <= iu1; i++)
	{
	  interval tangle = tan( ( PI() * i ) / ( 2 * nofrays ) );
	  res = power_point( cinterval( hx, hx * tangle ), n );
	  update_res( res, resxl, resxu, resyl, resyu );
	}
  }
}

/*!
\param z The complex interval
\param n The integer exponent
\return The computed complex interval

Power function for integer powers with optimal (save roundoff) accuracy.

\sa power_fast(const cinterval&,int)
*/
cinterval power( const cinterval& z, int n ) throw() //---- 040720 --
//
//  Power function for integer powers with optimal (save roundoff) accuracy.
//
{
  if( n == 0 )
    return cinterval( ONE_INTERVAL() );
  else if( n == 1 )
    return z;
  else if( n == -1 )
    return 1 / z;
  else if( n == 2 )
    return sqr(z);
  else
    //  n >= 3  or  n <= -2
    //
    //  n may have a large absolute value and z may be a large interval.
    //  In this case, we calculate function values at specific points
    //  and look for min and max.
    //
    //  An implementation with fewer function evaluations would be possible,
    //  at the expense of an unreadable source code.
    {
      interval abs_z = abs( z );

      if( n < 0 && Inf( abs_z ) == 0.0 )
	//  z contains 0
        cxscthrow (STD_FKT_OUT_OF_DEF("cinterval power(const cinterval& z, int n ); z contains 0."));
      if( Sup( abs_z ) == 0.0 )
        return cinterval( ZERO_INTERVAL(), ZERO_INTERVAL() );
      else
        {
	  real
	    irez = Inf(Re(z)),
	    srez = Sup(Re(z)),
	    iimz = Inf(Im(z)),
	    simz = Sup(Im(z));
	  interval
	    hxl( irez ), hxu( srez ), hyl( iimz ), hyu( simz );
	  real
	    resxl,  resxu, resyl, resyu;
	  cinterval res;
	  //
	  //  extremal values lie in corner points or on intersections of the
	  //  boundary of z with any ray  arg( w ) = k*Pi/(2*(n-1))
	  //
	  //  first check the corners of z
	  //
	  res = power_point( cinterval( hxl, hyl ), n );
	  resxl = Inf( Re(res) );
	  resxu = Sup( Re(res) );
	  resyl = Inf( Im(res) );
	  resyu = Sup( Im(res) );
	  res = power_point( cinterval( hxu, hyl ), n );
	  update_res( res, resxl, resxu, resyl, resyu );
	  res = power_point( cinterval( hxl, hyu ), n );
	  update_res( res, resxl, resxu, resyl, resyu );
	  res = power_point( cinterval( hxu, hyu ), n );
	  update_res( res, resxl, resxu, resyl, resyu );
	  //
	  //  consider the origin, if it is a boundary point
	  //  (negative n have been taken care of above)
	  //
          //	  if( irez * srez <= 0.0 && iimz * simz <= 0.0 && 
          //               irez * srez * iimz * simz == 0.0 )
	  if ( 0<=Re(z) && 0<=Im(z) && 
	       (irez==0 || srez==0 || iimz==0 || simz==0 ) )
	    update_res( cinterval( ZERO_INTERVAL(), ZERO_INTERVAL()), 
                                         resxl, resxu, resyl, resyu );
	  //
	  //  now check for ray intersections
	  //  for each quadrant and each boundary, we must consider
	  //  the function values at at most 4 consecutive intersections
	  //
	  //  Ia: lower boundary
	  //
	  if( iimz != 0.0 )
	    {
	      interval arg_h = arg( cinterval( Re(z), hyl ) );

	      //  check if there are ray intersections
	      horizontal_check( hyl, iimz, arg_h, irez, srez,
				resxl, resxu, resyl, resyu, n );
	    }
	  //  Ib: upper boundary
	  if( simz != 0.0 )
	    {
	      interval arg_h = arg( cinterval( Re(z), hyu ) );

	      //  check if there are ray intersections
	      horizontal_check( hyu, simz, arg_h, irez, srez,
				resxl, resxu, resyl, resyu, n );
	    }
	  //  IIa: left boundary
	  if( irez != 0.0 )
	    {
	      interval arg_h = arg( cinterval( hxl, Im(z) ) );

	      //  check if there are ray intersections
	      vertical_check( hxl, irez, arg_h, iimz, simz,
			      resxl, resxu, resyl, resyu, n );
	    }
	  //  IIb: right boundary
	  if( srez != 0.0 )
	    {
	      interval arg_h = arg( cinterval(  hxu, Im(z) ) );

	      //  check if there are ray intersections
	      vertical_check( hxu, srez, arg_h, iimz, simz,
			      resxl, resxu, resyl, resyu, n );
	    }
	  return cinterval( interval( resxl, resxu ), interval( resyl, resyu ) );
	}
    }
}
//
//-- end power ----------------------------------------------------------------

void times2pown(cinterval& x, int n) throw()
// Fast multiplication with 2^n
{
    interval R(Re(x)),I(Im(x));
    times2pown(R,n);
    times2pown(I,n);
    x = cinterval(R,I);
}


//-- pow ------------------------------------------------------------ 040627 --
//
//  Analytic power function for real interval exponent, based on Ln.
//

cinterval pow( const cinterval& z, const interval& p ) throw()
// Neue Version von pow(...) von 040627
{
    return exp( p*Ln(z) );
}

//
//-- end pow ------------------------------------------------------------------

//-- pow ------------------------------------------------------------ 040627 --
//
//  Analytic power function for complex interval exponent, based on Ln.
//

cinterval pow( const cinterval& z, const cinterval& p ) throw()
{
    return exp( p*Ln(z) );
}

//
//-- end pow ------------------------------------------------------------------


//-- pow_all -------------------------------------------------------- 041013 --
//
//  Non-analytic power function for real interval exponent.
//
//  If 0 \not\in z, then compute four rectangular intervals that comprehend
//  an annulus which contains all values  zeta^phi, zeta in z, phi in p.
//
//  If 0 in z and negative reals in p, then abort execution
//  (potential modification: return the entire complex plane).
//
std::list<cinterval> pow_all( const cinterval& z, const interval& p ) throw()
{
  interval abs_z = abs( z );

  if( 0.0 < Inf( abs_z ) )
    {
      interval abs_z_p = exp( p * ln( abs_z ) );

      //  Inner and outer radii of the annulus are inf/sup( abs_z_n )
      //  Inscribed square has side length sqrt( 2 ) * rad_1
      interval rad_1 = INV_SQRT_2() * Inf( abs_z_p );
      interval rad_2 = interval(Sup( abs_z_p ));

      std::list<cinterval> res;

      //  4 intervals covering the annulus, counter-clockwise
      res.push_back( cinterval( interval(  Inf( rad_1 ),  Sup( rad_2 ) ),
				interval( -Sup( rad_1 ),  Sup( rad_2 ) ) ) );
      res.push_back( cinterval( interval( -Sup( rad_2 ),  Sup( rad_1 ) ),
				interval(  Inf( rad_1 ),  Sup( rad_2 ) ) ) );
      res.push_back( cinterval( interval( -Sup( rad_2 ), -Inf( rad_1 ) ),
				interval( -Sup( rad_2 ),  Sup( rad_1 ) ) ) );
      res.push_back( cinterval( interval( -Sup( rad_1 ),  Sup( rad_2 ) ),
				interval( -Sup( rad_2 ), -Inf( rad_1 ) ) ) );

      return res;
    }
  else
    //  0 in z
    {
      if( Inf( p ) > 0.0 )
	//
	//  All values   zeta^phi, zeta in z, phi in p   lie in a disk
	//
	{
	  interval abs_z_p = exp( p * ln( interval( Sup( abs_z ) ) ) );
	  real rad_p = Sup( abs_z_p );

	  std::list<cinterval> res;

	  res.push_back( cinterval( interval( -rad_p, rad_p ),
				    interval( -rad_p, rad_p ) ) );

	  return res;
	}
      else
        {
	//
	//  The set   zeta^phi, zeta in z, phi in p   is unbounded
	//  if inf( p ) < 0.  0^p is undefined for p <= 0.
	//
	  cxscthrow(STD_FKT_OUT_OF_DEF("pow_all(cinterval, interval); 0^p is undefined for p <= 0."));
	  std::list<cinterval> res;
          return res;
        } 
    }
}
//
//-- end pow_all --------------------------------------------------------------



//-----------------------------------------------------------------------------
//
//  Section 5: asin, acos, asinh, acosh
//
//  The implementation of acos, asinh, and acosh is based on asin:
//
//  acos( z )  = pi / 2 - asin( z )
//  asinh( z ) = i * asin( - i * z )
//  acosh( z ) = i * acos( z ) = +/- i * ( pi / 2 - asin( z ) )
//
//  Only principal values are computed.
//
//-----------------------------------------------------------------------------

//-- asin ----------------------------------------------------------- 040905 --
//
//  Analytic inverse sine function
//
interval f_aux_asin( const interval& x, const interval& y ) //------- 040905 --
//
//  auxiliary function for evaluating asin
//  f_aux_asin(z) = ( |z+1| + |z-1| ) / 2
//
{
    interval res;
//  interval sqr_y = sqr( y );
//  interval res = ( sqrt( sqr( x + 1.0 ) + sqr_y ) + sqrt( sqr( x - 1.0 ) + sqr_y ) ) / 2.0;
    if (y == 0.0 && abs(x) == 1.0) res = 1.0; 
    else res = (sqrtx2y2(x+1.0,y) + sqrtx2y2(x-1.0,y)) / 2.0;  // Blomquist;

    if ( Sup(res)==Infinity ) // Blomquist: program stop, if overflow occurs.
	cxscthrow (STD_FKT_OUT_OF_DEF("cinterval asin( const cinterval& z); z out of range"));

  //  correct overestimation of lower boundary
  //  (x is a point interval)
    real hlb = max( 1.0, abs( Sup( x ) ) );
    if( Inf( res ) < hlb )
      //  invalid overestimation!
	res = interval( hlb, Sup(res) );
  return res;
}


interval f_aux_asin_Vn(const interval& x, const interval& y)
{ // normal calculation of V;
    interval V,f1,f2;
    f1 = x+1.0;  f2 = x-1;
    V = abs(f1)*sqrtp1m1(sqr(y/f1)) + abs(f2)*sqrtp1m1(sqr(y/f2));
    times2pown(V,-1);
    return V;
} // f_aux_asin_Vn

interval ACOSH_p1(const interval& x, const interval& y)
{ // Calculating an inclusion for acosh(1+V/2) if |x|<1;
    const int p = -80;
    real r1(Inf(y)),t;
    int ex(expo(r1));
    interval res(0.0),u,V;

    if (ex>-2000 && ex <= -80) {
	u = abs(y)/sqrt1mx2(x);
	t = pred(pred(Inf(u)));
	res = interval(t,Sup(u));
    } else 
	if (ex>p) {
	    V = f_aux_asin_Vn(x,y); // usual calculation
	    res = acoshp1(V);
	}
    return res;
} // ACOSH_p1

interval ACOSH_f_aux( const interval& x, const interval& y )
// Calculating acosh( f_aux_asin(x,y) ); x,y: point intervals !!
// Blomquist, 22.04.2005;
{
    interval res,delta;
    real rx(abs(Inf(x))), ry(abs(Inf(y)));

    if (rx>2.0 || ry>2.0) res = acosh( f_aux_asin(x,y) ); // as before!
    else { // Now with improvements !!
	if (rx == 1.0) {
	    delta = abs(y);
	    if (expo(Inf(delta))<=-50) {
		res = sqrt(delta);
		rx = Sup(res);
		if (rx>0) res = interval(pred(Inf(res)),succ(rx));
	    } else {
	    times2pown(delta,-1); // delta: |y|/2;
	    delta = sqrtp1m1(sqr(delta)) + delta;
	    res = acoshp1(delta);
	    }
	}
	else 
	    if (rx<1.0) res = ACOSH_p1(x,y);
	    else res = acoshp1( (abs(x)-1.0) + f_aux_asin_Vn(x,y) );
    }
    return res;
} // ACOSH_f_aux


interval Asin_beta( const interval& x, const interval& y )
// Calculating the improved real part of asin(z); Blomquist 19.06.2005;
// Re(asin(z)) = asin[ 2x/(sqrt[(x+1)^2+y^2] + sqrt[(x-1)^2+y^2]) ] = asin[beta]
// Improvements for beta --> +1  and  beta --> -1  are necessary because of nearly
// vertical tangents of the real asin(t)-function for |t|-->+1;   
{
    const real c1 = 0.75;
    bool neg_b;
    real Infxa;
    interval res,beta,abs_beta,delta,tm,tp,u,v,xa;
    beta = x / ( (sqrtx2y2(1+x,y) + sqrtx2y2(1-x,y))/2 );
    if (Inf(beta)<-1) Inf(beta)=-1;
    if (Sup(beta)> 1) Sup(beta)=1; 
    abs_beta = abs(beta);
    if (Inf(abs_beta)<c1) res = asin(beta); // Normal calculation 
    else { // Inf(abs_beta)>=c1; Calculation now with improvements: 
	xa = x;
	neg_b = Inf(x)<0;
	if (neg_b) xa = -xa; // Inf(xa) >0 :
	Infxa = Inf(xa);
	if (Infxa > 1) {
	    tm = y/(xa-1); 
	    tp = y/(xa+1);
	    u = sqrtp1m1(sqr(tm));
	    v = sqrtp1m1(sqr(tp));
	    delta = (tm*tp - u*v) / (2+u+v);
	} else
	    if (Infxa == 1) {
		u = abs(y);
		times2pown(u,-1); // u = |y|/2
		delta = u - sqrtp1m1(sqr(u));
	    } else {
		tp = 1+xa;  tm = 1-xa;
		delta = tm*(sqrt(1+sqr(y/tm))+1) - tp*sqrtp1m1(sqr(y/tp));
		times2pown(delta,-1);
	    }
	res = HALFPI() - asin( sqrt(delta*(2-delta)) );
	if (neg_b) res = -res;
    }
    return res;
}

cinterval asin( const cinterval& z ) throw() //---------------------- 040730 --
{
    const real gr = 6.355804e307; // upper bound for abs(rez),abs(imz)
    interval
	rez = Re(z),
	imz = Im(z);

    real irez = Inf(rez),
	 srez = Sup(rez),
	 iimz = Inf(imz),
	 simz = Sup(imz);

    interval hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

    real resxl, resxu, resyl, resyu;

    bool bl    = (iimz< 0.0) && (simz>0.0),
         raxis = (iimz==0.0) && (simz==0);

  //
  //  1st: check for singularities
  //
    if( (irez<-1 && (bl || (iimz<0  && simz==0))) || 
	(srez >1 && (bl || (iimz==0 && simz>0))) )
    cxscthrow(STD_FKT_OUT_OF_DEF("cinterval asin( const cinterval& z ); z contains singularities."));
  //
  //  check for too large bounds of abs(rez) and abs(imz) to prevent 
  //  overflow by calculating f_aux_asin(...)
  //
  resxl = max(abs(irez),abs(srez));  
  resxu = max(abs(iimz),abs(simz));
  if (resxl>gr || resxu>gr) 
      cxscthrow(STD_FKT_OUT_OF_DEF("cinterval asin( const cinterval& z ); z with too large bounds."));
  //
  //  2nd: real part
  //
  if( iimz < 0.0 && simz > 0.0 )
    //  z intersects [-1,1]
    {
      if( irez <= 0.0 )
	resxl = Inf( asin( hxl ) );
      else
//	resxl = Inf( asin( hxl / f_aux_asin( hxl, interval( max( - iimz, simz ) ) ) ) );
	resxl = Inf( Asin_beta(hxl,interval( max(-iimz,simz) )) ); // Blomquist, 19.06.2005;
      if( srez < 0.0 )
//	resxu = Sup( asin( hxu / f_aux_asin( hxu, interval( max( - iimz, simz ) ) ) ) );
	resxu = Sup( Asin_beta(hxu,interval( max(-iimz,simz) )) ); // Blomquist, 19.06.2005;
      else
	resxu = Sup( asin( hxu ) );
     }
  else
    {
      if( ( iimz >= 0.0 && irez >= 0.0  ) || ( simz <= 0.0 && irez <= 0.0 ) )
	//  left boundary in quadrants I or III
	//  min( Re( z ) ) in upper left corner
	//  resxl = Inf( asin( hxl / f_aux_asin( hxl, hyu ) ) );
	  resxl = Inf( Asin_beta(hxl,hyu) ); // Blomquist, 19.06.2005;
      else
	//  left boundary in quadrants II or IV
	//  min( Re( z ) ) in lower left corner
	//  resxl = Inf( asin( hxl / f_aux_asin( hxl, hyl ) ) );
	  resxl = Inf( Asin_beta(hxl,hyl) ); // Blomquist, 19.06.2005;
      if( ( iimz >= 0.0 && srez >= 0.0  ) || ( simz <= 0.0 && srez <= 0.0 ) )
	//  right boundary in quadrants I or III
	//  max( Re( z ) ) in lower right corner
	//  resxu = Sup( asin( hxu / f_aux_asin( hxu, hyl ) ) );
	  resxu = Sup( Asin_beta(hxu,hyl) ); // Blomquist, 19.06.2005;
      else
	//  right boundary in quadrants II or IV
	//  max( Re( z ) ) in upper right corner
	//  resxu = Sup( asin( hxu / f_aux_asin( hxu, hyu ) ) );
	  resxu = Sup( Asin_beta(hxu,hyu) ); // Blomquist, 19.06.2005;
    }
  //
  //  3rd: imaginary part
  //
  if (raxis) { // Interval argument is now a subset of the real axis.
               // Blomquist, 16.06.2005;
      if (srez<0.0) resyl =  Inf( ACOSH_f_aux( hxu, hyu )); 
      else          resyl = -Sup( ACOSH_f_aux( hxu, hyu )); 
      if (irez>0.0) resyu = -Inf( ACOSH_f_aux( hxl, hyu ));
      else          resyu =  Sup( ACOSH_f_aux( hxl, hyu ));
  } else 
  if( simz <= 0.0 )
    //  z in lower half plane
    //  min( Im( z ) ) in point with max |z|
    //  max( Im( z ) ) in point with min |z|
    {
	if( irez < -srez )
	//  most of z in quadrant III
	{
	  resyl = - Sup( ACOSH_f_aux( hxl, hyl ) );
	  if( srez < 0.0 )
	    resyu = - Inf( ACOSH_f_aux( hxu, hyu ) );
	  else
	    resyu = - Inf( ACOSH_f_aux( ZERO_INTERVAL(), hyu ) );
	}
      else
	//  most of z in quadrant IV
	{
	  resyl = - Sup( ACOSH_f_aux( hxu, hyl ) );
	  if( irez > 0.0 )
	    resyu = - Inf( ACOSH_f_aux( hxl, hyu ) );
	  else
	    resyu = - Inf( ACOSH_f_aux( ZERO_INTERVAL(), hyu ) );
	}
    }
  else if( iimz >= 0.0 )
    //  z in upper half plane
    //  min( Im( z ) ) in point with min |z|
    //  max( Im( z ) ) in point with max |z|
    {   
	if( irez < -srez )  // if( irez + srez < 0.0 )
	//  most of z in quadrant II
	{
	  resyu = Sup( ACOSH_f_aux( hxl, hyu ) );
	  if( srez < 0.0 )
	    resyl = Inf( ACOSH_f_aux( hxu, hyl ) );
	  else
	    resyl = Inf( ACOSH_f_aux( ZERO_INTERVAL(), hyl ) );
	}
      else
	//  most of z in quadrant I
	{
	  resyu = Sup( ACOSH_f_aux( hxu, hyu ) );
	  if( irez > 0.0 )
	     resyl = Inf( ACOSH_f_aux( hxl, hyl ) );
	  else
	    resyl = Inf( ACOSH_f_aux( ZERO_INTERVAL(), hyl ) );
	}
    }
  else
    //  z intersects imaginary axes
    //  min( Im( z ) ) in point in lower half plane with max |z|
    //  max( Im( z ) ) in point in upper half plane with max |z|
    {
      if( irez < -srez ) // if( irez + srez < 0.0 )
	//  most of z in quadrants II and IV
	{
	  resyl = - Sup( ACOSH_f_aux( hxl, hyl ) );
	  resyu =   Sup( ACOSH_f_aux( hxl, hyu ) );
	}
      else
	{
	  resyl = - Sup( ACOSH_f_aux( hxu, hyl ) );
	  resyu =   Sup( ACOSH_f_aux( hxu, hyu ) );
	}
    }

  return cinterval( interval( resxl, resxu ), interval( resyl, resyu ) );

}
//
//-- end asin -----------------------------------------------------------------



interval Acos_beta( const interval& x, const interval& y )
// Calculating the improved real part of acos(z); Blomquist 05.06.2005;
// Re(acos(z)) = acos[ 2x/(sqrt[(x+1)^2+y^2] + sqrt[(x-1)^2+y^2]) ]
{
    const real c1 = 0.75;
    interval res(0),beta,delta,tm,tp,u,v,xa;
    real Infy(Inf(y)),Infx(Inf(x));
    beta = x / ( (sqrtx2y2(1+x,y) + sqrtx2y2(1-x,y))/2 );
    if (Inf(beta)<-1) Inf(beta)=-1;
    if (Sup(beta)> 1) Sup(beta)= 1; 

    if (Sup(beta)<c1) 
	if (Sup(beta)<-c1) { // Improvement for beta --> -1
	    xa = -x;       // Inf(xa)>0:
	    Infx = -Infx;  // Infx > 0:
	    if (Infx > 1) {
		tm = y/(xa-1);  tp = y/(xa+1);
		u = sqrtp1m1(sqr(tm));
		v = sqrtp1m1(sqr(tp));
		delta = (tm*tp - u*v) / (2+u+v);
	    } else
		if (Infx == 1) {
		    u = abs(y);
		    times2pown(u,-1);  // u = |y|/2
		    delta = u - sqrtp1m1(sqr(u));
		} else {
		    tp = 1+xa;  tm = 1-xa;
		    delta = tm*(sqrt(1+sqr(y/tm))+1) - tp*sqrtp1m1(sqr(y/tp));
		    times2pown(delta,-1);
		}
	    res = PI() - asin( sqrt(delta*(2-delta)) );
	} else res = acos(beta); // Normal calculation
    else  // Sup(beta)>=c1
	if (Infx>1) 
	{
	    tm = y/(x-1);
	    if (expo(Sup(tm))<=-27) {
		if (Infy!=0) {
		    u = abs(Infy) / sqrtx2m1(x);
		    res = interval( pred(Inf(u)),succ(succ(Sup(u))) );
		}
	    } else {
		tp = y/(x+1);
		u = sqrtp1m1(sqr(tm));
		v = sqrtp1m1(sqr(tp));
		delta = (tm*tp - u*v) / (2+u+v);
		res = asin(sqrt(delta*(2-delta)));
	    }
	} else
	    if (Infx==1) {
		if (expo(Inf(y))<=-52) {
		u = sqrt(abs(y));
		if (Sup(u)==0) res = 0;
		else res = interval(pred(Inf(u)),succ(succ(Sup(u))));
		} else {
		    u = abs(y);
		    times2pown(u,-1); // u = |y|/2
		    delta = u - sqrtp1m1(sqr(u));
		    res = asin( sqrt(delta*(2-delta)) );
		}
	    } else {
		tp = 1+x;  tm = 1-x;
		delta = tm*(sqrt(1+sqr(y/tm))+1) - tp*sqrtp1m1(sqr(y/tp));
		times2pown(delta,-1);
		res = asin( sqrt(delta*(2-delta)) );
	    }

    return res;
} 


//-- acos ----------------------------------------------------------- 040730 --
//
// cinterval acos( const cinterval& z )
// {
//   w := acos(z);
//   Re(w) in a new Version,
//   Im(w) = -Im(asin(z));  Blomquist, 14.06.2005;
// }
//
//-- acos -----------------------------------------------------------------


//-- acos ----------------------------------------------------------- 040730 --
//
cinterval acos( const cinterval& z ) throw()  //--------------------- 040730 --
{
    const real gr = 6.355804e307; // upper bound for abs(rez),abs(imz)
    interval 
	rez = Re(z),
	imz = Im(z);

    real
	irez = Inf(rez),
	srez = Sup(rez),
	iimz = Inf(imz),
	simz = Sup(imz);

    interval
	hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

    bool bl    = (iimz< 0.0) && (simz>0.0),
         raxis = (iimz==0.0) && (simz==0);
  real
    resxl, resxu, resyl, resyu;
  //
  //  1st: check for singularities
  //
    if( (irez<-1 && (bl || (iimz<0  && simz==0))) || 
	(srez >1 && (bl || (iimz==0 && simz>0))) )
    cxscthrow(STD_FKT_OUT_OF_DEF("cinterval acos( const cinterval& z ); z contains singularities."));
  //
  //  check for too large bounds of abs(rez) and abs(imz) to prevent 
  //  overflow by calculating f_aux_asin(...)
  //
  resxl = max(abs(irez),abs(srez));  
  resxu = max(abs(iimz),abs(simz));
  if (resxl>gr || resxu>gr) 
      cxscthrow(STD_FKT_OUT_OF_DEF("cinterval acos( const cinterval& z ); z with too large bounds."));
  //
  //  2nd: real part
  //
  //  Blomquist, 06.06.2005;
      if( iimz < 0.0 && simz > 0.0 )
      //  z intersects [-1,1] on the x-axis
      {
	  if( irez <= 0.0 ) resxu = Sup( acos( hxl ) );
	  else resxu = Sup( Acos_beta(hxl,interval( max(-iimz,simz) )) );

	  if( srez < 0.0 ) 
	       resxl = Inf( Acos_beta(hxu,interval(max(-iimz,simz))) );
	  else resxl = Inf( acos( hxu ) ); 
      }
      else
      {
	  if (irez<0 && srez>0) 
	      // z intersects the posizive or negative y-axis
	      if (iimz >= 0) {
		  resxl = Inf( Acos_beta(hxu,hyl) );
		  resxu = Sup( Acos_beta(hxl,hyl) );
	      } else {
		  resxl = Inf( Acos_beta(hxu,hyu) );
		  resxu = Sup( Acos_beta(hxl,hyu) );
	      }
	  else 
	  {
	  if( ( iimz >= 0.0 && irez >= 0.0  ) || ( simz <= 0.0 && irez < 0.0 ) )
	  //  left boundary in quadrants I or III
	  //  min( Re( z ) ) in lower right corner
	      resxl = Inf( Acos_beta(hxu,hyl) );
	  else
	  //  left boundary in quadrants II or IV
	  //  min( Re( z ) ) in upper right corner
	      resxl = Inf( Acos_beta(hxu,hyu) );

	  if( ( iimz >= 0.0 && srez > 0.0  ) || ( simz <= 0.0 && srez <= 0.0 ) )
	  //  right boundary in quadrants I or III
	  //  max( Re( z ) ) in upper left corner
	      resxu = Sup( Acos_beta(hxl,hyu) );
	  else
	  //  right boundary in quadrants II or IV
	  //  max( Re( z ) ) in lower left corner
	      resxu = Sup( Acos_beta(hxl,hyl) );
	  }
      }
  //
  //  3rd: imaginary part
  //
  if (raxis) { // Interval argument is now a subset of the real axis.
               // Blomquist, 16.06.2005;
      if (srez<0.0) resyl =  Inf( ACOSH_f_aux( hxu, hyu )); 
      else          resyl = -Sup( ACOSH_f_aux( hxu, hyu )); 
      if (irez>0.0) resyu = -Inf( ACOSH_f_aux( hxl, hyu ));
      else          resyu =  Sup( ACOSH_f_aux( hxl, hyu ));
  } else 
  if( simz <= 0.0 )
    //  z in lower half plane
    //  min( Im( z ) ) in point with max |z|
    //  max( Im( z ) ) in point with min |z|
    {
      if( irez + srez < 0.0 )
	//  most of z in quadrant III
	{
	  resyl = - Sup( ACOSH_f_aux( hxl, hyl ) );
	  if( srez < 0.0 )
	    resyu = - Inf( ACOSH_f_aux( hxu, hyu ) );
	  else
	    resyu = - Inf( ACOSH_f_aux( ZERO_INTERVAL(), hyu ) );
	}
      else
	//  most of z in quadrant IV
	{
	  resyl = - Sup( ACOSH_f_aux( hxu, hyl ) );
	  if( irez > 0.0 )
	    resyu = - Inf( ACOSH_f_aux( hxl, hyu ) );
	  else
	    resyu = - Inf( ACOSH_f_aux( ZERO_INTERVAL(), hyu ) );
	}
    }
  else if( iimz >= 0.0 )
    //  z in upper half plane
    //  min( Im( z ) ) in point with min |z|
    //  max( Im( z ) ) in point with max |z|
    {   
	if( irez < -srez )  // if( irez + srez < 0.0 )
	//  most of z in quadrant II
	{
	  resyu = Sup( ACOSH_f_aux( hxl, hyu ) );
	  if( srez < 0.0 )
	    resyl = Inf( ACOSH_f_aux( hxu, hyl ) );
	  else
	    resyl = Inf( ACOSH_f_aux( ZERO_INTERVAL(), hyl ) );
	}
      else
	//  most of z in quadrant I
	{
	  resyu = Sup( ACOSH_f_aux( hxu, hyu ) );
	  if( irez > 0.0 )
	     resyl = Inf( ACOSH_f_aux( hxl, hyl ) );
	  else
	    resyl = Inf( ACOSH_f_aux( ZERO_INTERVAL(), hyl ) );
	}
    }
  else
    //  z intersects imaginary axes
    //  min( Im( z ) ) in point in lower half plane with max |z|
    //  max( Im( z ) ) in point in upper half plane with max |z|
    {
      if( irez < -srez ) // if( irez + srez < 0.0 )
	//  most of z in quadrants II and IV
	{
	  resyl = - Sup( ACOSH_f_aux( hxl, hyl ) );
	  resyu = Sup( ACOSH_f_aux( hxl, hyu ) );
	}
      else
	{
	  resyl = - Sup( ACOSH_f_aux( hxu, hyl ) );
	  resyu = Sup( ACOSH_f_aux( hxu, hyu ) );
	}
    }

  return cinterval( interval( resxl, resxu ), -interval( resyl, resyu ) );

}
//
//-- end acos -----------------------------------------------------------------


//-- asinh ---------------------------------------------------------- 040730 --
//
cinterval asinh( const cinterval& z ) throw()
//
//  asinh( Z ) = i * asin( -i * z )
//
{
  cinterval res = asin( cinterval( Im(z), -Re(z) ) );
  return cinterval( -Im(res), Re(res) );
}
//
//-- end asinh ----------------------------------------------------------------


//-- acosh ---------------------------------------------------------- 040908 --
//
cinterval acosh( const cinterval& z ) throw()
//
//  acosh( z ) = i * acos( z ) = +/- i * ( pi / 2 - asin( z ) )
//
{
  interval
    rez = Re(z),
    imz = Im(z);

  real
    irez = Inf(rez),
    srez = Sup(rez),
    iimz = Inf(imz),
    simz = Sup(imz);

  interval
    hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

  real
    resxl, resxu, resyl, resyu;

  // cinterval res;

  //
  //  1st: check for singularities
  //
  if( ( iimz <= 0.0 && simz >= 0.0 ) && ( irez < 1.0 ) )
      cxscthrow(STD_FKT_OUT_OF_DEF("cinterval acosh( const cinterval& z ); z contains singularities."));
  //  With this restriction the complex interval argument and the real axis must not have any common 
  //  point, if irez < +1;
  //  So for example the negative real axis must not be touched from above if irez<1, although this
  //  should be possible if the principal branch is considered! So the above restriction is too widely in
  //  some cases;  Blomquist, 21.06.2005;
  //
  //  2nd: z in upper half plane (or on the real axis)
  //  acosh( z ) =  + i * ( pi / 2 - asin( z ) )
  //
  if( iimz > 0.0 )
  {
      cinterval res = acos(z);
      return cinterval( -Im(res),Re(res) );
  }
  //
  //  3rd: z in lower half plane
  //  acosh( z ) =  - i * ( pi / 2 - asin( z ) )
  //
  if( simz < 0.0 )
    {
//      cinterval res = HALFPI() - asin( z );
	cinterval res = acos(z);  // Blomquist, 14.06.2005
	return cinterval( Im(res), -Re(res) );
    }
  //
  //  z intersects [1,infinity)
  //
  //  real part
  //  minimum on the left on real axes, maximum in lower or upper right corner
  //
  resxl = Inf( acosh( hxl ) );
  interval ytilde( max( -iimz, simz ) );
//  resxu = Sup( acosh( f_aux_asin( hxu, ytilde ) ) );
  resxu = Sup( ACOSH_f_aux(hxu,ytilde) ); // Blomquist, 14.06.2005;
  //
  //  imaginary part
  //  minimum in lower left corner, maximum in upper left corner
  //
  //   resyl = -Sup( acos( hxl / f_aux_asin( hxl, hyl ) ) );
  resyl = -Sup( Acos_beta(hxl,hyl) ); // Blomquist, 14.06.2005;
  //  resyu =  Sup( acos( hxl / f_aux_asin( hxl, hyu ) ) );
  resyu = Sup( Acos_beta(hxl,hyu) );

  return cinterval( interval( resxl, resxu ), interval( resyl, resyu ) );

}
//
//-- end acosh ----------------------------------------------------------------


//-----------------------------------------------------------------------------
//
//  Section 6: atan, acot, atanh, acoth
//
//  The implementation of acot, atanh, and acoth is based on atan:
//
//  acot( z )  =  atan( 1 / z )
//  atanh( z ) = - i * atan( i * z )
//  acoth( z ) = i * acot( i * z )
//
//  Only principal values are computed.
//
//-----------------------------------------------------------------------------

//-- atan ----------------------------------------------------------- 040912 --
//
//  Analytic inverse tangent function
//

void re_vert( const real& x, const interval& hx,
	      const real& rew_inf, const real& rew_sup,
 	      real& resxl, real& resxu ) //---------------------- 040729 --
//
//  Subroutine of analytic inverse tangent function.
//  Evaluate real part on a vertical boundary.
//
{
  if( x == 0.0 )
    //  singularities have been handled before, hence Re( w ) > 0
    {
      resxl = 0.0;
      resxu = 0.0;
    }
  else
    {
      if( x > 0.0 )
	//  w in quadrants I and/or II
	//  atan is the inverse function of tan(t), t in (-pi/2,pi/2).
	{
	  resxl = rew_sup > 0.0 ? Inf( Atan( 2 * hx,rew_sup )/2.0 )
		    : ( rew_sup < 0.0 ? Inf( (Atan( 2*hx,rew_sup ) + PI() )/2.0 )
		                          : Inf( HALFPI()/2.0 ) );

	  resxu = rew_inf > 0.0 ? Sup( Atan( 2*hx,rew_inf )/2.0 )
		        : ( rew_inf < 0.0 ? Sup( (Atan( 2*hx,rew_inf ) + PI())/2.0 )
			                  : Sup( HALFPI()/2.0 ) );
	}
      else
	//  w in quadrants III and/or IV
	{
	  resxl = rew_inf < 0.0 ? Inf( (Atan( 2*hx,rew_inf ) - PI())/2.0 )
                            : ( rew_inf > 0.0 ? Inf( Atan( 2*hx,rew_inf )/2.0 )
			                      : -Sup( HALFPI()/2.0 ) );
	  resxu = rew_sup < 0.0 ? Sup( (Atan( 2*hx,rew_sup ) - PI())/2.0 )
                            : ( rew_sup > 0.0 ? Sup( Atan( 2*hx,rew_sup )/2.0 )
			                      : -Inf( HALFPI()/2.0 ) );
	}
    }
} //  re_vert

interval Aux_1_atan(const real& x)
// x>=0;
// Calculating: ln[ 1+2/(sqrt(1+x^2)-1) ], [x] = x,
// [x] is a point interval !
// Blomquist; 19.02.05;
{
const int exOv = +54;
const int exUn = -26;

interval res,
    ix(x),  // ix is point interval with x>=0;
    r,t;
 int ex(expo(x));

 if (ex>=exOv) { // preventing overflow
     r = 2/ix;
     t = r*pred(1.0);
     r = r*succ(1.0);
     res = interval(Inf(t),Sup(r));
 } else
     if (ex<=exUn) { // x < 2^(-27)
         res = U_atan - 2*ln(ix);
     } else { // normal calculation
	 t = sqrtp1m1( sqr(ix) ); // t = sqrt(1+x^2)-1
	 res = lnp1(2/t); // res = ln[1 + 2/(sqrt(1+x^2)-1) ]
     }
 return res;
} // Aux_1_atan

interval Q_atan_UPSIGN(const interval& x, const interval& y)
{
// x: abs(Re(z));  x is real interval
// y: Inf(Im(z));  y is point interval
// Q_atan_UPSIGN: ln[ 1 + 4y/(x^2+(1-y)^2) ]

    const int n = 511; 
    interval res,t,t1,t2;
    int ex_x,ex,s;
    if (y==1.0) {
	if (Inf(x)>1.0) {
	    t = 2/x;
	    res = lnp1(sqr(t));
	} else
	    if(Sup(x)<1) res = ln(4+sqr(x)) - 2*ln(x);
	    else { // Punkt 3.:
		t = interval(Sup(x));  t = 2/t;
		t1 = lnp1(sqr(t)); 
		t = interval(Inf(x));
		t2 = ln(4.0+sqr(t)) - 2*ln(t);
		res = interval(Inf(t1),Sup(t2));  
	    }
    } else { // y <> [1,1]
	ex_x = expo( Sup(abs(x)) );
	ex   = expo( Sup(abs(y)) );
	if (ex_x>ex) ex = ex_x;  // Maximum 
	if (ex>n) { // scaling:
	    s = n-ex-1;
	    t = x;  times2pown(t,s);  // fast scaling with 2^(s)
	    t1 = y; times2pown(t1,s); // fast scaling with 2^(s)
	    t2 = sqr(t) + sqr(comp(0.5,s+1)-t1); // t2: denominator
	    t = y / t2; // scaled quotient 
	    times2pown(t,2*s+2); // back-scaling with 2^(s+2); '+2': factor 4 !!
	    res = lnp1(t);
	} else res = lnp1(4*y/(sqr(x)+sqr(1-y))); // normal calculation
    }
    return res;
} // Q_atan_UPSIGN

cinterval atan( const cinterval& z ) throw() //----- 040912 --
{
  interval
    rez = Re(z),
    imz = Im(z);

  real
    irez = Inf(rez),
    srez = Sup(rez),
    iimz = Inf(imz),
    simz = Sup(imz);

  const int n = 511; // For possible scaling

  interval
    hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

  real
    resxl, resxu, resyl, resyu;
  //
  //  1st: check for singularities
  //
  if( ( irez <= 0.0 && srez >= 0.0 ) && ( iimz <= -1.0 || simz >= 1.0 ) )
    cxscthrow(STD_FKT_OUT_OF_DEF("cinterval atan( const cinterval& z ); z contains singularities."));
  //
  //  2nd: real part
  //  Re( atan( z ) ) = Arg( w ) / 2, where w = 1 - x^2 - y^2 + i * 2x )
  //
  //  evaluate atan on vertical boundaries
  //
  interval
//      y_sqr = sqr( imz ),
//      rew_l = (1 - y_sqr) - sqr( hxl ),  // Blomquist; before: rew_l = 1 - sqr(hxl) - y_sqr, 
//      rew_u = (1 - y_sqr) - sqr( hxu );  // Blomquist; before: rew_u = 1 - sqr(hxu) - y_sqr; 
      rew_l, rew_u;

/*  ------------------------------ Blomquist ---------------------------------------------------  */
/*  ----------       Improvements for Im(z) = [1,1]  or  Im(z) = [-1,-1]        ------------------*/
  bool sqrImz_1 = (iimz==simz) && (iimz==1.0 || iimz==-1.0); // Test for Im(z) = [1,1] or [-1,-1]

  if (sqrImz_1) { 
      rew_l = -abs(hxl);  hxl = interval(sign(irez)); 
      rew_u = -abs(hxu);  hxu = interval(sign(srez));
  }
  else {
      int ex,s;
      interval imz_, scf;
      int ex1 = expo(iimz);  int ex2 = expo(simz); 
      if (ex2>ex1) ex1 = ex2;

      ex = expo(irez);
      if(ex1>ex) ex = ex1; // Maximum
      if (ex>n) { // Scaling necessary
	  s = n - ex - 1;
	  scf = interval(comp(0.5,s+1)); // scf: scaling factor 2^s
	  times2pown(scf,s); // scf = 2^(2*s);
	  times2pown(hxl,s); // hxl = hxl * 2^s
	  imz_ = imz;
	  times2pown(imz_,s); // imz_ = imz_ * 2^s
	  rew_l = (scf - sqr(imz_)) - sqr(hxl); // here now without overflow!!
	  times2pown(hxl,s); // hxl = hxl * 2^s
      } else rew_l = (1 - sqr( imz )) - sqr( hxl );

      ex = expo(srez);
      if(ex1>ex) ex = ex1; // Maximum
      if (ex>n) { // Scaling necessary
	  s = n - ex - 1;
	  scf = interval(comp(0.5,s+1)); // scf: scaling factor 2^s
	  times2pown(scf,s); // scf = 2^(2*s);
	  times2pown(hxu,s); // hxu = hxu * 2^s
	  imz_ = imz;
	  times2pown(imz_,s); // imz_ = imz_ * 2^s
	  rew_u = (scf - sqr(imz_)) - sqr(hxu); // here now without overflow!!
	  times2pown(hxu,s); // hxu = hxu * 2^s
      } else rew_u = (1 - sqr( imz )) - sqr( hxu );
  }
/*  ------------------------------ Blomquist; 22.02.05; ----------------------------------------  */

  //
  //  left boundary
  //
  real rew_inf = Inf( rew_l );
  real rew_sup = Sup( rew_l );
  re_vert( irez, hxl, rew_inf, rew_sup, resxl, resxu );

  //
  //  right boundary
  //
  rew_inf = Inf( rew_u );
  rew_sup = Sup( rew_u );
  real res_l, res_u;
  re_vert( srez, hxu, rew_inf, rew_sup, res_l, res_u );

  resxl = min( resxl, res_l );
  resxu = max( resxu, res_u );
  //
  //  look for extremal values on horizontal boundaries
  //  since atan( x+iy ) = atan( x-iy ),
  //  intersections can be considered in the upper half plane
  //
  real abs_y_min = Inf( abs( imz ) );

  if( abs_y_min > 1.0 )
    {
      interval
	abs_hyl = interval( abs_y_min ),
//      abs_hxl = sqrt( sqr( abs_hyl ) - 1.0 );
        abs_hxl = sqrtx2m1(abs_hyl);  // Blomquist;

      if( Sup( abs_hxl ) > irez && Inf( abs_hxl ) < srez )
	//  extremal curve intersects lower boundary of x+i|y| in quadrant I
	//  intersection in Q I or Q IV: update minimum
	//	resxl = inf( atan( abs_y_min / abs_hxl ) ) / 2.0;
	resxl = Inf( (PI() - atan( 1.0 / abs_hxl ))/2.0 );
      else if( -Inf( abs_hxl ) > irez && -Sup( abs_hxl ) < srez )
	//  extremal curve intersects lower boundary of x+i|y| in quadrant II
	//  intersection in Q II or Q III: update maximum
	resxu = Sup( (atan( 1.0 / abs_hxl ) - PI())/2.0 );
    }
  //  3rd: imaginary part
  //  Im( atan( z ) ) = +/- Ln( 1 +/- 4y/( x^2 + (1 -/+ y)^2 ) ) / 4
  //
  //  evaluate atan on horizontal boundaries
  interval
    abs_rez = abs(rez), 
    im_atan_l, im_atan_u;

  if( iimz < 0.0 )
//    im_atan_l = -ln( 1 - 4 * hyl / ( x_sqr + sqr( 1 + hyl ) ) ) / 4.0;
//    im_atan_l = -lnp1(-4 * hyl / ( x_sqr + sqr( 1 + hyl ) )) / 4.0;  // Blomquist
      im_atan_l = -Q_atan_UPSIGN(abs_rez,-hyl);  // Blomquist (Versuch)
  else
//    im_atan_l = ln( 1 + 4 * hyl / ( x_sqr + sqr( 1 - hyl ) ) ) / 4.0;
    im_atan_l = Q_atan_UPSIGN(abs_rez,hyl);  // Blomquist
  times2pown(im_atan_l,-2);

  if( simz < 0.0 )
//    im_atan_u = -ln( 1 - 4 * hyu / ( x_sqr + sqr( 1 + hyu ) ) ) / 4.0;
//    im_atan_u = -lnp1(-4 * hyu / ( x_sqr + sqr( 1 + hyu ) ) ) / 4.0; // Blomquist
      im_atan_u = -Q_atan_UPSIGN(abs_rez,-hyu);  // Blomquist
  else
//    im_atan_u = ln( 1 + 4 * hyu / ( x_sqr + sqr( 1 - hyu ) ) ) / 4.0;
      im_atan_u = Q_atan_UPSIGN(abs_rez,hyu);  // Blomquist
  times2pown(im_atan_u,-2);

  resyl = min( Inf( im_atan_l ), Inf( im_atan_u ) );
  resyu = max( Sup( im_atan_l ), Sup( im_atan_u ) );
  //
  //  look for extremal values on vertical boundaries,
  //  if vertical boundaries intersect extremal curves
  //
  real abs_x_min = Inf( abs( rez ) );
  interval
    x_extr = interval( abs_x_min ),
//    y_extr = sqrt( 1.0 + sqr( x_extr ) );
    y_extr = sqrt1px2(x_extr);                     // Blomquist;

  if( Inf( y_extr ) < simz && Sup( y_extr ) > iimz )
    //  extremal curve intersects left boundary of |x|+iy in quadrant I
    //  update maximum
    //  resyu = Sup( ln( 1 + 4 * y_extr / ( sqr( x_extr ) + sqr( 1 - y_extr ) ) ) ) / 4.0;
    //  resyu = Sup( Aux_1_atan(abs_x_min)/4.0 );    // Blomquist
  {
      rez = Aux_1_atan(abs_x_min);
      times2pown(rez,-2);
      resyu = Sup(rez);
  }

  if( -Sup( y_extr ) < simz && -Inf( y_extr ) > iimz )
    //  extremal curve intersects left boundary of |x|+iy in quadrant IV
    //  update minimum
    //  resyl = -Sup( ln( 1 + 4 * y_extr / ( sqr( x_extr ) + sqr( 1 - y_extr ) ) ) ) / 4.0;
    //  resyl = -Sup( Aux_1_atan(abs_x_min)/4.0 );  // Blomquist
  {
      rez = Aux_1_atan(abs_x_min);
      times2pown(rez,-2);
      resyl = -Sup(rez);
  }

  return cinterval( interval( resxl, resxu ), interval( resyl, resyu ) );

}
//
//-- end atan -----------------------------------------------------------------


//-- acot ----------------------------------------------------------- 040912 --
//
//  Analytic inverse cotangent function
//  acot( z ) = atan( 1/z )
//  The code of acot( z ) is almost identical to the code of atan( z )
//
cinterval acot( const cinterval& z ) throw()
{
  interval
    rez = Re(z),
    imz = Im(z);

  real
    irez = Inf(rez),
    srez = Sup(rez),
    iimz = Inf(imz),
    simz = Sup(imz);

  const int n = 511; // For possible scaling

  interval
    hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

  real
    resxl, resxu, resyl, resyu;
  //
  //  1st: check for singularities
  //
  if( ( (irez <= 0.0) && (srez >= 0.0) ) && ( (iimz <= 1.0) && (simz >= -1.0) ) )
    cxscthrow(STD_FKT_OUT_OF_DEF("cinterval acot( const cinterval& z ); z contains singularities."));
  //
  //  2nd: real part
  //  Re( atan(  z  ) )   = Arg( w ) / 2, where w = 1 - x^2 - y^2 + i * 2x )
  //  Re( atan( 1 / z ) ) = Arg( w ) / 2, where w = x^2 + y^2 - 1 + i * 2x )
  //
  //  evaluate acot on vertical boundaries
  //
  interval
//    y_sqr = sqr( imz ),
//    rew_l = (y_sqr - 1) + sqr(hxl),
//    rew_u = (y_sqr - 1) + sqr(hxu);
//    rew_l = (sqr( hxl )-1) + y_sqr, 
//    rew_u = (sqr( hxu )-1) + y_sqr; 
      rew_l, rew_u;
/*  ------------------------------ Blomquist ---------------------------------------------------  */
/*  ----------       Improvements for Im(z) = [1,1]  or  Im(z) = [-1,-1]        ------------------*/
  bool sqrImz_1 = (iimz==simz) && (iimz==1.0 || iimz==-1.0); // Test for Im(z) = [1,1] or [-1,-1]

  if (sqrImz_1) { 
      rew_l = abs(hxl);  hxl = interval(sign(irez)); 
      rew_u = abs(hxu);  hxu = interval(sign(srez));
 }
  else {
      int ex,s;
      interval imz_, scf;
      int ex1 = expo(iimz);  int ex2 = expo(simz); 
      if (ex2>ex1) ex1 = ex2;

      ex = expo(irez);
      if(ex1>ex) ex = ex1; // Maximum
      if (ex>n) { // Scaling necessary
	  s = n - ex - 1;
	  scf = interval(comp(0.5,s+1)); // scf: scaling factor 2^s
	  times2pown(scf,s); // scf = 2^(2*s);
	  times2pown(hxl,s); // hxl = hxl * 2^s
	  imz_ = imz;
	  times2pown(imz_,s); // imz_ = imz_ * 2^s
	  rew_l = (sqr(imz_) - scf) + sqr(hxl); // here now without overflow!!
	  times2pown(hxl,s); // hxl = hxl * 2^s
      } else rew_l = (sqr( imz ) - 1.0) + sqr( hxl );

      ex = expo(srez);
      if(ex1>ex) ex = ex1; // Maximum
      if (ex>n) { // Scaling necessary
	  s = n - ex - 1;
	  scf = interval(comp(0.5,s+1)); // scf: scaling factor 2^s
	  times2pown(scf,s); // scf = 2^(2*s);
	  times2pown(hxu,s); // hxu = hxu * 2^s
	  imz_ = imz;
	  times2pown(imz_,s); // imz_ = imz_ * 2^s
	  rew_u = (sqr(imz_) - scf) + sqr(hxu); // here now without overflow!!
	  times2pown(hxu,s); // hxu = hxu * 2^s
      } else rew_u = (sqr( imz )-1.0) + sqr( hxu );
  }
/*  ------------------------------ Blomquist; 22.02.05; ----------------------------------------  */

  //
  //  left boundary
  //
  real rew_inf = Inf( rew_l );
  real rew_sup = Sup( rew_l );
  re_vert( irez, hxl, rew_inf, rew_sup, resxl, resxu );
  //
  //  right boundary
  //
  rew_inf = Inf( rew_u );
  rew_sup = Sup( rew_u );
  real res_l, res_u;
  re_vert( srez, hxu, rew_inf, rew_sup, res_l, res_u );

  resxl = min( resxl, res_l );
  resxu = max( resxu, res_u );
  //
  //  look for extremal values on horizontal boundaries
  //  since acot( x+iy ) = acot( x-iy ),
  //  intersections can be considered in the upper half plane
  //
  real abs_y_min = Inf( abs( imz ) );

  if( abs_y_min > 1.0 )
    {
      interval
	abs_hyl = interval( abs_y_min ),
//	abs_hxl = sqrt( sqr( abs_hyl ) - 1.0 );
        abs_hxl = sqrtx2m1(abs_hyl);  // Blomquist;
      if( Sup( abs_hxl ) > irez && Inf( abs_hxl ) < srez )
	//  extremal curve intersects lower boundary of x+i|y| in quadrant I
	//  intersection in Q I or Q IV: update maximum
	resxu = Sup( atan( 1.0 / abs_hxl ) / 2.0 );
      if( -Inf( abs_hxl ) > irez && -Sup( abs_hxl ) < srez )
	//  extremal curve intersects lower boundary of x+i|y| in quadrant II
	//  intersection in Q II or Q III: update minimum
	resxl = -Sup( atan( 1.0 / abs_hxl ) / 2.0 );
     }
  //
  //  3rd: imaginary part
  //  Im( atan( z ) ) = +/- Ln( 1 +/- 4y/( x^2 + (1 -/+ y)^2 ) ) / 4
  //  Im( acot ) = - Im ( atan ): We calculate Im( atan ) and return "-"
  //
  //  evaluate atan on horizontal boundaries
  //
  interval
//  x_sqr = sqr( rez ), // overflow is avoided by calling Q_atan_UPSIGN(...)
    im_atan_l, im_atan_u,
    abs_rez = abs(rez);  // Blomquist;
  if( iimz < 0.0 )
//    im_atan_l = -ln( 1 - 4 * hyl / ( x_sqr + sqr( 1 + hyl ) ) ) / 4.0;
      im_atan_l = -Q_atan_UPSIGN(abs_rez,-hyl);  // Blomquist
  else
//    im_atan_l = ln( 1 + 4 * hyl / ( x_sqr + sqr( 1 - hyl ) ) ) / 4.0;
      im_atan_l = Q_atan_UPSIGN(abs_rez,hyl);  // Blomquist
  times2pown(im_atan_l,-2);

  if( simz < 0.0 )
//    im_atan_u = -ln( 1 - 4 * hyu / ( x_sqr + sqr( 1 + hyu ) ) ) / 4.0;
      im_atan_u = -Q_atan_UPSIGN(abs_rez,-hyu);  // Blomquist
  else
//    im_atan_u = ln( 1 + 4 * hyu / ( x_sqr + sqr( 1 - hyu ) ) ) / 4.0;
      im_atan_u = Q_atan_UPSIGN(abs_rez,hyu);  // Blomquist
  times2pown(im_atan_u,-2);

  resyl = min( Inf( im_atan_l ), Inf( im_atan_u ) );
  resyu = max( Sup( im_atan_l ), Sup( im_atan_u ) );
  //
  //  look for extremal values on vertical boundaries,
  //  if vertical boundaries intersect extremal curves
  //
  real abs_x_min = Inf( abs( rez ) );
  interval
    x_extr = interval( abs_x_min ),
//    y_extr = sqrt( 1.0 + sqr( x_extr ) );
    y_extr = sqrt1px2(x_extr);  // Blomquist
  if( Inf( y_extr ) < simz && Sup( y_extr ) > iimz )
    //  extremal curve intersects left boundary of |x|+iy in quadrant I
    //  update maximum
    //  resyu = Sup( ln( 1 + 4 * y_extr / ( sqr( x_extr ) + sqr( 1 - y_extr ) ) ) ) / 4.0;
    //  resyu = Sup( Aux_1_atan(abs_x_min)/4.0 );    // Blomquist
  {
      rez = Aux_1_atan(abs_x_min);
      times2pown(rez,-2);
      resyu = Sup(rez);
  }

  if( -Sup( y_extr ) < simz && -Inf( y_extr ) > iimz )
    //  extremal curve intersects left boundary of |x|+iy in quadrant IV
    //  update minimum
    //  resyl = -Sup( ln( 1 + 4 * y_extr / ( sqr( x_extr ) + sqr( 1 - y_extr ) ) ) ) / 4.0;
    //  resyl = -Sup( Aux_1_atan(abs_x_min)/4.0 );  // Blomquist
  {
      rez = Aux_1_atan(abs_x_min);
      times2pown(rez,-2);
      resyl = -Sup(rez);
  }

  return cinterval( interval( resxl, resxu ), interval( -resyu, -resyl ) );

}
//
//-- end acot -----------------------------------------------------------------


//-- atanh ---------------------------------------------------------- 040912 --
//
cinterval atanh( const cinterval& z ) throw()
//
//  atanh( z ) = - i * atan( i * z )
//
{
  cinterval res = atan( cinterval( -Im(z), Re(z) ) );
  return cinterval( Im(res), -Re(res) );
}
//
//-- end atanh ----------------------------------------------------------------

//-- acoth ---------------------------------------------------------- 040912 --
//
cinterval acoth( const cinterval& z ) throw()
//
//  acoth( z ) = i * acot( i * z )
//
{
  cinterval res = acot(  cinterval( -Im(z), Re(z) ) );
  return cinterval( -Im(res), Re(res) );
}
//
//-- end acoth ----------------------------------------------------------------


cinterval sqr(const cinterval& z) throw()
// Improvement of the sqr(z)-function; Blomquist, 24.06.2005;
{   
    interval rez(Re(z)), reza(abs(rez)),
	     imz(Im(z)), imza(abs(imz));
    real
	irez = Inf(reza),
	srez = Sup(reza),
	iimz = Inf(imza),
	simz = Sup(imza);
    interval
	hxl(irez), hxu(srez), hyl(iimz), hyu(simz);
    real
	resxl, resxu;

    resxl = Inf( (hxl-hyu)*(hxl+hyu) );
    resxu = Sup( (hxu-hyl)*(hxu+hyl));  

    hxl = rez * imz;
    times2pown(hxl,1); // fast multiplikation with 2;

    return cinterval( interval(resxl,resxu), hxl );
}

} // namespace cxsc

/*

  End of File: cimath.cpp

*/
