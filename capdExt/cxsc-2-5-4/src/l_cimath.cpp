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

/* CVS $Id: l_cimath.cpp,v 1.21 2014/01/30 17:23:46 cxsc Exp $ */

/*
**
**  File: l_cimath.cpp, 2007/03/05.
**
**  Copyright (C) Markus Neher,        markus.neher@math.uni-karlsruhe.de
**                Ingo Eble,           ingoeble@web.de
**                Frithjof Blomquist,  Blomquist@math.uni-wuppertal.de
*/

//Include header files
#include <l_real.hpp>
#include <l_cimath.hpp>  // Complex functions in staggered format
#include <l_imath.hpp>   // "l_interval" standard functions

namespace cxsc{

    inline l_real min(const l_real& a, const l_real& b)
    {
	return (a<b)? a : b;
    }

    inline l_real max(const l_real& a, const l_real& b)
    {
	return (a>b)? a : b;
    }

    inline const l_interval& ONE_INTERVAL()
    {
	static const l_interval one = l_interval(1.0);
	return one;
    }

bool disjoint(const l_interval& x, const l_interval& y)
{
  l_real ix( Inf(x) ),iy( Inf(y) ),sx( Sup(x) ), sy( Sup(y) );
  l_real inf( ( ix > iy )? ix : iy );
  l_real sup( ( sx < sy )? sx : sy );

  return ( inf > sup );
}


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


l_cinterval exp(const l_cinterval& z) throw()
{
    int stagsave = stagprec,
	stagmax = 19;
    l_interval lreal(Re(z)), limg(Im(z));
    cinterval dz(z);
    l_cinterval y;
    bool grob = ( (Sup(Re(dz)) > succ(succ(Inf(Re(dz))))) || 
                  (Sup(Im(dz)) > succ(succ(Inf(Im(dz))))) );

    if (stagprec==1 || grob) y = exp(dz); // schnelle Auswertung!
    else
    {
	if (stagprec < stagmax) stagprec++;
	else stagprec = stagmax;
	l_interval
	    A( exp(lreal) ),
	    B( limg );
	y = l_cinterval( A*cos( B ) , A*sin( B ) );
	stagprec = stagsave;
	y = adjust(y);
    }
    return y;
}

l_cinterval exp2(const l_cinterval& z) throw()
{
	return exp(z*Ln2_l_interval());
}

l_cinterval exp10(const l_cinterval& z) throw()
{
	return exp(z*Ln10_l_interval());
}

l_cinterval expm1(const l_cinterval& z) throw()
// exp(z) - 1;
// Calculates nearly optimal inclusions for not too wide intervals z.
// Blomquist, 03.12.2008;
{
	int stagsave = stagprec,
 stagmax = 30; // l_imath.cpp: sqrt() uses stagmax = 30;
 if (stagprec>stagmax) stagprec = stagmax;
	
 const l_interval cancl_test = l_interval(0.995,1.005);
 l_interval rez(Re(z)), imz(Im(z));
 l_interval exp_x, sin_y, cos_y, h, Xt;
 l_cinterval res;
	
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
 res = l_cinterval(h,imz);
	
 stagprec = stagsave;
 res = adjust(res);
		
 return res;
}

l_cinterval cos(const l_cinterval& z) throw()
{
    int stagsave = stagprec,
	stagmax = 19;
    l_interval lreal(Re(z)),limg(Im(z));
    cinterval dz(z);
    l_cinterval y;
    bool grob = ( (Sup(Re(dz)) > succ(succ(Inf(Re(dz))))) || 
                  (Sup(Im(dz)) > succ(succ(Inf(Im(dz))))) );

    if (stagprec==1 || grob) y = cos(dz); // schnelle Auswertung!
    else
    {
	if (stagprec < stagmax) stagprec++;
	else stagprec = stagmax;
	l_interval
	    A( lreal ),
	    B( limg );
	y = l_cinterval( cos( A )*cosh( B ) , -sin( A )*sinh( B ) );
	stagprec = stagsave;
	y = adjust(y);
    }
    return y;
}

l_cinterval sin(const l_cinterval& z) throw()
{
    int stagsave = stagprec,
	stagmax = 19;
    l_interval lreal(Re(z)),limg(Im(z));
    cinterval dz(z);
    l_cinterval y;
    bool grob = ( (Sup(Re(dz)) > succ(succ(Inf(Re(dz))))) || 
                  (Sup(Im(dz)) > succ(succ(Inf(Im(dz))))) );

    if (stagprec==1 || grob) y = sin(dz); // schnelle Auswertung!
    else
    {
	if (stagprec < stagmax) stagprec++;
	else stagprec = stagmax;
	l_interval
	    A( lreal ),
	    B( limg );
	y = l_cinterval( sin( A )*cosh( B ) , cos( A )*sinh( B ) );
	stagprec = stagsave;
	y = adjust(y);
    }
    return y;
}

l_cinterval cosh(const l_cinterval& z) throw()
{
    int stagsave = stagprec,
	stagmax = 19;
    l_interval lreal(Re(z)),limg(Im(z));
    cinterval dz(z);
    l_cinterval y;
    bool grob = ( (Sup(Re(dz)) > succ(succ(Inf(Re(dz))))) || 
                  (Sup(Im(dz)) > succ(succ(Inf(Im(dz))))) );

    if (stagprec==1 || grob) y = cosh(dz); // schnelle Auswertung!
    else
    {
	if (stagprec < stagmax) stagprec++;
	else stagprec = stagmax;
	l_interval
	    A( lreal ),
	    B( limg );
	y = l_cinterval( cos( B )*cosh( A ) , sin( B )*sinh( A ) );
	stagprec = stagsave;
	y = adjust(y);
    }
    return y;
}

l_cinterval sinh(const l_cinterval& z) throw()
{
    int stagsave = stagprec,
	stagmax = 19;
    l_interval lreal(Re(z)),limg(Im(z));
    cinterval dz(z);
    l_cinterval y;
    bool grob = ( (Sup(Re(dz)) > succ(succ(Inf(Re(dz))))) || 
                  (Sup(Im(dz)) > succ(succ(Inf(Im(dz))))) );

    if (stagprec==1 || grob) y = sinh(dz); // schnelle Auswertung!
    else
    {
	if (stagprec < stagmax) stagprec++;
	else stagprec = stagmax;
	l_interval
	    A( lreal ),
	    B( limg );
	y = l_cinterval( cos( B )*sinh( A ) , sin( B )*cosh( A ) );
	stagprec = stagsave;
	y = adjust(y);
    }
    return y;
}

l_cinterval sqr(const l_cinterval& z) throw()
// Blomquist, 21.11.2006;
{   
    dotprecision akku;
    l_interval rez(Re(z)), reza(abs(rez)),
	imz(Im(z)), imza(abs(imz));
    l_real
	irez = Inf(reza),
	srez = Sup(reza),
	iimz = Inf(imza),
	simz = Sup(imza);

    akku = 0;
    accumulate(akku,irez,irez);
    accumulate(akku,-simz,simz);
    irez = cxsc::rnd_down(akku);

    akku = 0;
    accumulate(akku,srez,srez);
    accumulate(akku,-iimz,iimz);
    srez = cxsc::rnd_up(akku); 

    rez = rez * imz;
    times2pown(rez,1); // fast multiplikation with 2;

    return l_cinterval( l_interval(irez,srez), rez );
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
     const l_interval& hy, const l_interval& cosh_2y, l_real irez, l_real srez,
     const l_interval& hxl, 
     const l_interval& hxu, l_real& resxl, l_real& resxu )
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

  if (srez - irez > Inf( Pid2_l_interval() ))
    //  2 intersections
    both = true;
  else
    {
      l_interval
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
	      l_interval
		sin_2xl = sin( 2 * hxl ),
		sin_2xu = sin( 2 * hxu );

	      if( !disjoint( l_interval(0), res_l ) )
		//  intersection on the left boundary
		{
		  if( Inf( sin_2xl ) >= 0.0 )
		    // "left" intersection
		    {
		      left = true;
		      //  remove the intersection by changing res_l!
		      res_l = -l_interval(1);
		    }
		  else
		    {
		      if( Sup( sin_2xl ) <= 0.0 )
			// "right" intersection
			{
			  right = true;
			  //  remove the intersection by changing res_l!
			  res_l =  l_interval(1);
			}
		      else
			//  zero is interior point of sin_2xl
			//  if the real sine function has optimal precision,
			//  this case should never happen
			both = true;
		    }
		}

	      if( !disjoint( l_interval(0), res_u ) )
		//  intersection on the right boundary
		{
		  if( Inf( sin_2xu ) >= 0.0 )
		    // "left" intersection
		    {
		      left = true;
		      //  remove the intersection by changing res_u!
		      res_u = l_interval(1);
		    }
		  else
		    {
		      if( Sup( sin_2xu ) <= 0.0 )
			// "right" intersection
			{
			  right = true;
			  //  remove the intersection by changing res_u!
			  res_u = -l_interval(1);
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
  l_interval re_tan = 1 / sinh( 2 * abs( hy ) );

  //  "left" intersection, sin( 2x ) > 0
  if( left || both )
    {
      resxl = min( resxl, Inf( re_tan ) );
      resxu = max( resxu, Sup( re_tan ) );
    }

  //  "right" intersection, sin( 2x ) < 0
  if( right || both )
    {
      resxl = min( resxl, -Sup( re_tan ) );
      resxu = max( resxu, -Inf( re_tan ) );
    }
} // end horizontal_check


l_cinterval Tan( const l_cinterval& z ) throw() // -------------- 040827 --
{
    l_cinterval y;

    l_interval
	rez = Re(z),   // rez = z.re(),
	imz = Im(z);   // imz = z.im();

    l_real
	irez = Inf(rez),
	srez = Sup(rez),
	iimz = Inf(imz),
	simz = Sup(imz);

    l_interval
	hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

    l_real
	resxl, resxu, resyl, resyu;
  //
  //  1st: check for poles
  //
    if( ( !disjoint( l_interval(0), imz ) ) && 
	( !disjoint( l_interval(0), cos( rez ) ) ) )
	cxscthrow (STD_FKT_OUT_OF_DEF("l_cinterval tan( const l_cinterval& z); Pole(s) in z"));
  //
  //  2nd: real part
  //
  //  evaluate tan on vertical boundaries
  //
    l_interval
	cos_2rez   = cos( 2 * rez ),
	sinh_imz_2 = sqr( sinh( imz ) );

    l_interval
	re_tan_l=sin( 2*hxl ) / ( 2*( sqr( cos( hxl ) ) + sinh_imz_2 ) ),
	re_tan_u=sin( 2*hxu ) / ( 2*( sqr( cos( hxu ) ) + sinh_imz_2 ) );

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
	l_interval
	    cosh_2yl = - 1 / cosh( 2 * hyl ),
	    cosh_2yu = - 1 / cosh( 2 * hyu );

	if( !disjoint( cos_2rez, cosh_2yl ) && iimz != 0.0 )
	//extremal curve intersects lower boundary
	    horizontal_check(hyl,cosh_2yl,irez,srez,hxl,hxu,resxl,resxu);

	if( !disjoint( cos_2rez, cosh_2yu ) && simz != 0.0 )
	//extremal curve intersects upper boundary
	    horizontal_check(hyu,cosh_2yu,irez,srez,hxl,hxu,resxl,resxu);
    }
  //
  //  3rd: imaginary part
  //
  //  evaluate tan on horizontal boundaries
  //
    l_interval
	cos_rez_2 = sqr( cos( rez ) );

    l_interval
	im_tan_l = sinh( 2*hyl ) / (2*(cos_rez_2 + sqr( sinh( hyl ) ))),
	im_tan_u = sinh( 2*hyu ) / (2*(cos_rez_2 + sqr( sinh( hyu ) )));

    resyl = min( Inf( im_tan_l ), Inf( im_tan_u ) );
    resyu = max( Sup( im_tan_l ), Sup( im_tan_u ) );

  //
  //  look for extremal values on vertical boundaries
  //  here, the situation is simpler than for the real part
  //  if 2 intersections with extremal curves appear ,
  //  one lies in the lower half plane, the other in the upper half plane
  //
    l_interval
	cos_2xl = cos( 2 * hxl ),
	cos_2xu = cos( 2 * hxu );
    l_interval im_tan;

    if( iimz < 0.0 )
    //  intersection in lower half plane?
    {
	l_interval
	    imz_h = l_interval( iimz, min( simz, l_real(0.0) ) ),
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
	l_interval
	    imz_h = l_interval( max( iimz, l_real(0.0) ), simz ),
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

    y = l_cinterval( l_interval(resxl,resxu ),l_interval(resyl,resyu ) );

    return y;

} // Tan

l_cinterval tan( const l_cinterval& z ) throw() 
{
// tan(z) has the poles z_k = pi*(1+2k)/2;   |k| in {0,1,2,3,...}.
// z = z_k + eps = pi*(1+2k)/2 + eps;  With  |eps|<<1  k can be calculated
// by:  k = Re(z)/pi - 0.5; With this k the komplex value eps is given
// by:  eps = z - pi*(1+2k)/2;  pi = 3.1415926... ;
// It holds: 
// tan(z) = tan(z_k+eps) = tan[pi*(1+2k)/2 + eps] 
//        = tan[pi/2 + pi*k + eps] = tan[pi/2 + eps] = -1 / tan(eps);
// Definitions:  u = Re(eps);  u = abs(u);  v = Im(eps);  v = abs(v);  
// if (Sup(u)<S && Sup(v)<S) tan(z) = -1 / Tan(eps);
// else tan(z) = Tan(z);     S = 1e-15;
// Thus, near the poles tan(z) is calculated in higher accuracy with
// -1 / Tan(eps);
// Blomquist, 28.09.2007;
    const real S = 1e-15;
    int stagsave = stagprec,
	stagmax = 19; 
    if (stagprec < stagmax) stagprec++;
    else stagprec = stagmax;

    l_cinterval y,eps;
    l_interval rez = Re(z);

    interval re(rez),u,v;
    real x(mid(re));
    double dbr = _double(x), pi(3.14159265358979323);
    int k,s;
    dbr = dbr/pi - 0.5;
    s = sign(dbr);
    k = (s>=0)? CUTINT(dbr+0.5) : CUTINT(dbr-0.5);
    if (k<-2147483647)
	cxscthrow (STD_FKT_OUT_OF_DEF(
                 "l_cinterval tan(const l_cinterval& z); z out of range"));
    eps = z - Pid2_l_interval()*(1+2*k);

    u = Re(eps);  u = abs(u);
    v = Im(eps);  v = abs(v);

    if (Sup(u)<S && Sup(v)<S)
	y = -l_cinterval(1) / Tan(eps);
    else y = Tan(z);

    stagprec = stagsave;
    y = adjust(y);

    return y;
} // tan()

l_cinterval cot( const l_cinterval& z ) throw() 
{
// cot(z) has the poles z_k = k*pi;   |k| in {0,1,2,3,...}.
// z = z_k + eps = k*pi + eps;  With  |eps|<<1  k can be calculated
// by:  k = Re(z)/pi;   With this k the komplex value eps is given
// by:  eps = z - k*pi = z - (pi/2)*2k;    pi = 3.1415926... ;
// It holds: 
// cot(z) = cot(z_k+eps) = cot(k*pi + eps) 
//        = cot(eps) = 1 / tan(eps);
// Definitions:  u = Re(eps);  u = abs(u);  v = Im(eps);  v = abs(v);  
// if (Sup(u)<S && Sup(v)<S) cot(z) = 1 / tan(eps);
// else cot(z) = tan(pi/2 - z);     S = 1e-15;
// Thus, near the poles cot(z) is calculated in higher accuracy with
// 1 / tan(eps);
// Blomquist, 29.09.2007;
    const real S = 1e-15;
    int stagsave = stagprec,
	stagmax = 19; 
    if (stagprec < stagmax) stagprec++;
    else stagprec = stagmax;

    l_cinterval y,eps;

    l_interval 	rez = Re(z);
    interval re(rez),u,v;
    real x(mid(re));

    double dbr = _double(x), pi(3.14159265358979323);
    int k,s;
    dbr = dbr/pi;
    s = sign(dbr);
    k = (s>=0)? CUTINT(dbr+0.5) : CUTINT(dbr-0.5);
    if (k<-2147483647)
	cxscthrow (STD_FKT_OUT_OF_DEF(
                 "l_cinterval cot(const l_cinterval& z); z out of range"));
    eps = z - Pid2_l_interval()*2*k;
    u = Re(eps);  u = abs(u);
    v = Im(eps);  v = abs(v);
    if (Sup(u)<S && Sup(v)<S)
	y = cinterval(1) / Tan(eps);
    else y = Tan(Pid2_l_interval() - z);

    stagprec = stagsave;
    y = adjust(y);

    return y;
} // cot



//-- tanh ---------------------------------------------------------
//
//  tanh( z ) = transp( i * tan( transp( i * z ) )
//
l_cinterval tanh( const l_cinterval& z ) throw()
{
    int stagsave = stagprec,
	stagmax = 19; 
    if (stagprec < stagmax) stagprec++;
    else stagprec = stagmax;

    l_cinterval res = tan( l_cinterval( Im(z), Re(z) ) ),y;
    y = l_cinterval( Im(res), Re(res) );

    stagprec = stagsave;
    y = adjust(y);

    return y;
}
//
//-- end tanh -------------------------------------------------------

//-- coth -----------------------------------------------------------
//
//   coth( z ) = i * cot( i * z );
//

l_cinterval coth( const l_cinterval& z ) throw()
{ // coth( z ) = i * cot( i * z );
    int stagsave = stagprec,
	stagmax = 19; 
    if (stagprec < stagmax) stagprec++;
    else stagprec = stagmax;

    l_cinterval zh = l_cinterval( -Im(z), Re(z) ); //  zh = i*z;
    l_cinterval res = cot(zh);
    zh = l_cinterval( -Im(res), Re(res) );

    stagprec = stagsave;
    zh = adjust(zh);

    return zh;
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

l_interval Atan(const l_interval& y, const l_interval& x) throw()
// Calculating an inclusion of atan(y/x) with x<>[0,0];
// This help function must only be used for POINT intervals y,x !!
// This function avoids internal overflow by calculating y/x.
// atan(t), t in R, is a point symmetrical function, so we only need 
// to consider t >= 0;
// For sufficiently great t the values atan(t) will be included by
// Pid2_l_interval(), declared in l_imath.hpp. Thus, for t > t0, we must
// fulfill   atan(t) > INF := Inf( Pid2_l_interval() ). For stagprec = 19
// it holds:  INF = pi/2 - d, d := 2.6817655261... * 10^-308;
// --->  atan(t) > INF = pi/2 - d  must be fulfilled for t > t0 = ?
// atan(t) is a monotonic function, so t0 is determined by
// atan(t0) = INF = pi/2 - d     --->  t0 = tan(pi/2 - d) = cot(d);
// With  cot(d) = 1/d - d/3 - (d^3)/45 - ...  < 1/d,  0 < d < pi,  we get:
//
// t > 1/d = 3.7288... * 10^+307 --> atan(t) > INF := Inf(Pid2_l_interval());
// 
// In the following algorithm with  ex_x,ex_y  it holds:
// Infy[1] = my * 2^ex_y;    Infx[1] = mx * 2^ex_x;   0.5 <= my,mx < 1;
// For  ex_y - ex_x >= N0:=1023  in this algorithm atan(y/x) is included by 
// Pid2_l_interval(); It is easy to show, that for ex_y - ex_x = N0:=1023 
// the minimal value of y/x is given by u := 2^1022 = 4.49423...*10^+307, 
// which is greater than 1/d = 3.7288... * 10^+307, so the inclusion of
// atan(y/x) with Pid2_l_interval() is correct!
// On the other side, for ex_y - ex_x <= 1022 the quotient y/x is smaller
// than 2^1023 < MaxReal = 0.9999...*2^1024, so that an overflow by
// calculating y/x is not possible!
// Blomquist, 04.12.2006;   
{
    const int c = 1022;
    l_interval res(0);
    l_real Infx(Inf(x)),
	Infy(Inf(y));
    int ex_x(expo_gr(Infx)),
	ex_y(expo_gr(Infy)),
	signx(sign(Infx)),
	signy(sign(Infy)),
	signq;
    if (signy!=0) {
	signq = signx * signy;
	if (ex_y-ex_x > c) 
	    res = signq>0 ? Pid2_l_interval() : -Pid2_l_interval();
	else res = atan(y/x);
    }

    return res;
}

l_interval Atan(const l_interval& y, const l_real& x) throw()
// Calculating an inclusion of atan(y/x) with x<>0.0;
// This help function must only be used for POINT intervals y !!
// This function avoids internal overflow by calculating y/x. 
{
    l_interval xi(x);
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
l_interval Arg( const l_cinterval& z ) throw()
{
    l_real 
	srez = Sup( Re(z) ),
	irez = Inf( Re(z) ), 
	simz = Sup( Im(z) ),
	iimz = Inf( Im(z) );

    l_interval
	hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

    l_real resl, resu;

    if( iimz > 0.0 )
    //  case I: Im(z) > 0
    {
	resl = ( srez > 0.0 ? Inf( Atan( hyl,hxu ) ) : 
	       ( srez < 0.0 ? Inf( Atan( hyu,hxu ) + Pi_l_interval() ) : 
		 Inf( Pid2_l_interval() ) ) );
	resu = ( irez > 0.0 ? Sup( Atan( hyu,hxl ) ) : 
		 ( irez < 0.0 ? Sup( Atan( hyl,hxl ) + Pi_l_interval() ) : 
		   Sup( Pid2_l_interval() ) ) );
	return l_interval( resl, resu );
    }
    else
    {
	if( simz < 0.0 )
	//  case II: Im(z) < 0
	{
	    resl = ( irez < 0.0 ? Inf( Atan( hyu,hxl ) - Pi_l_interval() ) : 
		     ( irez > 0.0 ? Inf( Atan( hyl,hxl ) ) : 
		       -Sup( Pid2_l_interval() ) ) );
	  resu = ( srez < 0.0 ? Sup( Atan( hyl,hxu ) - Pi_l_interval() ) : 
		   ( srez > 0.0 ? Sup( Atan( hyu,hxu ) ) : 
		     -Inf( Pid2_l_interval() ) ) );
	  return l_interval( resl, resu );
	}
	else
	// 0 in Im(z)
	{
	    if( irez > 0.0 )
	    //  case III: Re(z) > 0
	    //  z contains positive real values
	    {
		resl = iimz < 0.0 ? Inf( Atan( hyl,hxl ) ) : l_real(0.0);
		return l_interval( resl, Sup( Atan( hyu,hxl ) ) );
	    }
	    else
	    //  z contains nonpositive real numbers
	    {
		if( irez < 0.0 )
	        {
		  //  case IV: z contains negative real numbers
                  cxscthrow (STD_FKT_OUT_OF_DEF("l_interval Arg( const l_cinterval& z ); z contains negative real numbers"));
		  return l_interval(0.0);
		}
		else
		//  case V: 0 in z, but z doesn't contain negative real numbers
		{
		    if( srez > 0.0 )
		    //  diam( Re(z) > 0.0 )
		    {
			resl = iimz < 0.0 ? -Sup(Pid2_l_interval()) : 
				 l_real(0.0);
			resu = simz > 0.0 ?  Sup(Pid2_l_interval()) : 
				 l_real(0.0);
			return l_interval( resl, resu );
		    }
		    else
		    //  Re(z) == 0.0
		    {
			if( iimz == 0.0 && simz == 0.0 )
			//  Z == 0
			    return l_interval(0.0);
			else
			{
			   resl = ( iimz < 0.0 ? - Sup( Pid2_l_interval() ) : 
				    Inf( Pid2_l_interval() ) );
			   resu = ( simz > 0.0 ? Sup( Pid2_l_interval() ) : 
				    -Inf( Pid2_l_interval() ) );
			   return l_interval( resl, resu );
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

l_interval arg( const l_cinterval& z ) throw()
{
    l_real
	srez = Sup( Re(z) ),
	irez = Inf( Re(z) ),
	simz = Sup( Im(z) ),
	iimz = Inf( Im(z) );

    l_real resl, resu;

    if( irez < 0.0 && iimz <= 0.0 && simz >= 0.0  )
    //  z contains negative real values
    {
	if( srez > 0.0 )
	    //  0 in z and 0 interior point of Re(z)
	{
	    resl = ( iimz < 0.0 ? - Sup( Pi_l_interval() ) : l_real(0.0) );
	    resu = ( ( iimz < 0.0 && simz == 0.0 ) ? l_real(0.0) : 
		     Sup( Pi_l_interval() ) );
	    return l_interval( resl, resu );
	}
	else
	{ // srez <= 0.0
	    if( iimz == simz )
	    //  z is real interval containing no positive values
		return Pi_l_interval();
	    else
	    // sup( Re(z) ) <= 0, diam( Im(z) ) > 0
	    {
		if( srez == 0.0 )
		{
		    resl = ( simz > 0.0 ? Inf( Pid2_l_interval() ) : 
			     -Sup( Pi_l_interval() ) );
		    resu = ( iimz < 0.0 ? 
			     ( simz > 0.0 ? Sup( 3 * Pid2_l_interval() ) : 
			       -Inf( Pid2_l_interval() ) ) : 
			     Sup( Pi_l_interval() ) );
		    return l_interval( resl, resu );
		}
		else
		//   sup( Re(z) ) < 0, diam( Im(z) ) > 0
		{
		    l_interval hyl(iimz), hyu(simz);
		    resl = ( simz > 0.0 ? 
			     Inf( Atan( hyu,srez ) + Pi_l_interval() ) : 
			     -Sup( Pi_l_interval() ) );
		    resu = ( iimz < 0.0 ? ( simz > 0.0 ? 
					    Sup( Atan( hyl,srez ) + 
						 Pi_l_interval() ) : 
					    Sup( Atan( hyl,srez ) - 
						 Pi_l_interval() ) ) : 
			     Sup( Pi_l_interval() ) );
		    return l_interval( resl, resu );
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


// ***************************************************************************
// ***************************************************************************
// ***                      Multi-valued functions                        ****
// ***************************************************************************
// ***************************************************************************


//-- arg_inclmon: non-analytic inclusion-monotone argument function - 040617 --
//
//  (i)  arg_inclmon(Z) is defined for all Z in IC.
//  (ii) arg_inclmon(Z) = [-pi,pi] if Arg(Z) is not defined.
//

l_interval arg_inclmon( const l_cinterval& z ) throw()
{
    if( Inf( Re(z) ) < 0.0 && Inf( Im(z) ) <= 0.0 && Sup( Im(z) ) >= 0.0  )
	return l_interval( -Sup( Pi_l_interval() ),Sup( Pi_l_interval() ) );
    else return Arg(z);
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
l_cinterval Ln( const l_cinterval& z ) throw()
{ // Blomquist;
    int stagsave = stagprec,
	stagmax = 19; 
    if (stagprec < stagmax) stagprec++;
    else stagprec = stagmax;

    l_cinterval y;
    l_real
	srez = Sup( Re(z) ),
	simz = Sup( Im(z) ),
	iimz = Inf( Im(z) );
    l_interval a1( abs(Re(z)) ),
	a2( abs(Im(z)) );
    if ( Inf(a1) == 0.0 && Inf(a2) == 0.0 )
	//  z contains 0
	cxscthrow(STD_FKT_OUT_OF_DEF("l_cinterval Ln( const l_cinterval& z ); z contains 0"));
    if ( srez<0 && iimz<0 && simz>=0 ) 
	cxscthrow(STD_FKT_OUT_OF_DEF("l_cinterval Ln( const l_cinterval& z ); z not allowed")); 
    y = l_cinterval( ln_sqrtx2y2(Re(z),Im(z)), arg(z) );

    stagprec = stagsave;
    y = adjust(y);

   return y;
}
//
//-- end Ln -------------------------------------------------------------------

//-- ln: non-analytic natural logarithm ----------------------------- 040923 --
//
//  ln(z) is undefined if z contains zero.
//
l_cinterval ln( const l_cinterval& z ) throw()
{
    int stagsave = stagprec,
	stagmax = 19; 
    if (stagprec < stagmax) stagprec++;
    else stagprec = stagmax;

    l_cinterval y;
    l_interval a1( abs(Re(z)) ),
	a2( abs(Im(z)) );
    if ( Inf(a1) == 0.0 && Inf(a2) == 0.0 ) 
	//  z contains 0
	cxscthrow(STD_FKT_OUT_OF_DEF
              ("l_cinterval ln( const l_cinterval& z ); z contains 0"));
    y = l_cinterval( ln_sqrtx2y2(Re(z),Im(z)), arg(z) );

    stagprec = stagsave;
    y = adjust(y);

    return y; 
}
//
//-- end ln -------------------------------------------------------------------

l_cinterval lnp1(const l_cinterval& z) throw()
{  // ln(1+z);
	// Calculates nearly optimal inclusions for not too wide intervals z.
   // Blomquist, 03.12.2008;
	int stagsave = stagprec,
 	stagmax = 30; 
 	if (stagprec > stagmax) 
		stagprec = stagmax;
	
 	const real c1 = 1.0;
 	l_cinterval y;
 	l_interval abs_z(abs(z));
 	l_real
		srez = Sup( Re(z) ),
		simz = Sup( Im(z) ),
		iimz = Inf( Im(z) );
	
	if (-1 <= z) //  z contains -1
		cxscthrow(STD_FKT_OUT_OF_DEF(
			"l_cinterval lnp1(const l_cinterval& z); z contains -1"));
	if ( srez<-1 && iimz<0 && simz>=0 ) 
		cxscthrow(STD_FKT_OUT_OF_DEF(
			"l_cinterval lnp1(const l_cinterval& z); z not allowed"));
	
	if (Sup(abs_z) < c1)
	{
		abs_z = Re(z);
		abs_z = lnp1( abs_z*(2+abs_z) + sqr(Im(z)) );
		times2pown(abs_z,-1);
		y = l_cinterval( abs_z, arg(1+z) );
	}
	else
		y = Ln(1+z);
	stagprec = stagsave;
	y = adjust(y);

	return y;
}

l_cinterval log2( const l_cinterval& z ) throw()
{
	return Ln(z) / Ln2_l_interval();
}

l_cinterval log10( const l_cinterval& z ) throw()
{
	return Ln(z) / Ln10_l_interval();
}

//-----------------------------------------------------------------------------
//
//  Section 3: Root functions
//
//-----------------------------------------------------------------------------

l_interval Sqrt_zpx( const l_interval& x, const l_interval& y, int& d )
// Calculating:
//     res = sqrt(2^(-d)*|z|+2^(-d)*|x|) = 2^(-d/2) * sqrt(|z| + |x|)
//     without any internal overflow;
// Notice:  |z| = sqrt(x^2+y^2);   sqrt(|z|+|x|) < Maxreal; 
// Blomquist, 05.03.2007;
{
    const int c1 = 1021;
    l_real Infx(Inf(x)), Infy(Inf(y));
    int ex_x(expo_gr(Infx)), ex_y(expo_gr(Infy));
    l_interval xc(abs(x)),yc(y),res;
    bool yeq0(Infy==0);
    d = 0;

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
	if (yeq0) 
	{
	    times2pown(xc,1); 
	    res = sqrt(xc);
	} else 
	{
	    if (ex_y>ex_x) ex_x = ex_y; // ex_x = Max{ex_x,ex_y}
	    if (ex_x%2 != 0) ex_x--;
	    if (ex_x<-50)
	    {
		d = ex_x; // d is divisible by 2;
		Times2pown(xc,-d);
		Times2pown(yc,-d);
		res = sqrt( sqrtx2y2(xc,yc)+xc);
	    } else res = sqrt( sqrtx2y2(xc,yc)+xc);
	}
    return res;
} // Sqrt_zpx()

//-- sqrt: analytic square root ------------------------------------- 040626 --
//
l_interval Re_Sqrt_point( const l_interval& rez, const l_interval& imz ) 
//
//  Real part of analytic square root of A POINT INTERVAL ONLY.
//  Do not use this as a general function
//  - it's only a subroutine for sqrt and sqrt_all.
//  The calculation is void if (rez+I*imz) is not a complex number.
//  Blomquist, 05.03.2007;
{
    int d;
    l_interval res;
    l_real
	irez = Inf( rez ),
	iimz = Inf( imz );

    if( iimz == 0.0 )
    {
	if( irez >= 0.0 )
	    return sqrt( rez );
	else return l_interval(0.0);
    }
    else
    { // iimz <> 0 
	res = Sqrt_zpx(rez,imz,d);
	if (irez >= 0.0) 
	{
	    if (d!=0) times2pown(res,d/2);
	    return Sqrt2r_l_interval() * res;
	}
	else 
	{
	    iimz = abs(iimz);
	    if (d!=0) times2pown(iimz,-d/2);
	    return Sqrt2r_l_interval() * iimz / res;
	}
    }
}

l_interval Im_Sqrt_point( const l_interval& rez, const l_interval& imz ) 
//
//  Imaginary part of analytic square root of A POINT INTERVAL ONLY
//  Do not use this as a general function
//  - it's only a subroutine for sqrt and sqrt_all
//  The calculation is void if (rez+I*imz) is not a complex number
//  Blomquist, 05.03.2007;
{
    int d;
    l_interval res;
    l_real
	irez = Inf( rez ),
	iimz = Inf( imz );

    if( iimz == 0.0 )
    {
	if( irez >= 0.0 ) return l_interval(0.0);
	else return sqrt(-rez);
    }
    else
    {
	res = Sqrt_zpx(rez,imz,d);
	if( irez >= 0.0 )
	{
	    if (d!=0) times2pown(iimz,-d/2);
	    return Sqrt2r_l_interval() * iimz / res;
	}
	else
	{
	    if (d!=0) times2pown(res,d/2);
	    res = Sqrt2r_l_interval() * res;
	    if( iimz > 0.0 ) return res;
	    else return -res;
	}
    }
}

l_cinterval sqrt( const l_cinterval& z ) throw() // ----------------------
//
//  Analytic square root function with z in the principal branch.
//  The branch cut is the negative real axis. On the branch cut 
//  the function values are defined by comming from the II quadrant.
//  Blomquist, 23.06.2005;
//
{
    int stagsave = stagprec,
	stagmax = 19; 
    if (stagprec < stagmax) stagprec++;
    else stagprec = stagmax;

    l_cinterval y;
    l_real
	irez = Inf( Re(z) ),
	srez = Sup( Re(z) ),
	iimz = Inf( Im(z) ),
	simz = Sup( Im(z) );
    l_interval
	hxl( irez ), hxu( srez ), hyl( iimz ), hyu( simz );
    l_real
	resxl,  resxu, resyl, resyu;

    if( irez < 0.0 && iimz < 0.0 && simz >= 0.0 ) 
	//  if z is not in the principal branch then the inclusion monotony 
        // is violated!
	cxscthrow(STD_FKT_OUT_OF_DEF(
  "l_cinterval sqrt(const l_cinterval& z); z not in the principal branch."));

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
    y = l_cinterval( l_interval( resxl,resxu ), l_interval( resyl,resyu ) );

    stagprec = stagsave;
    y = adjust(y);

    return y;
}  

l_cinterval sqrtp1m1(const l_cinterval& z) throw()
// sqrt(1+z)-1;
// Calculates nearly optimal inclusions for not too wide intervals z.
// Blomquist, 08.07.2008;
{
	const real c = 0.125;
	int stagsave = stagprec,
 	stagmax = 30; // l_imath.cpp: sqrt() uses stagmax = 30;
 	if (stagprec>stagmax) stagprec = stagmax;
	
 	l_cinterval res;
 	l_interval absz(abs(z));
 	l_real Sup_absz(Sup(absz));
	
 	if (Sup_absz < c)
	 	res = z / (sqrt(z+1) + 1);
 	else 
	 	res = sqrt(z+1) - 1;
	
 	stagprec = stagsave;
 	res = adjust(res);	
 	return res;
}

l_cinterval sqrt1px2(const l_cinterval& z) throw()
// sqrt(1+z^2);
// Calculates nearly optimal inclusions for not too wide intervals z.
// Blomquist, 03.07.2008;
{
	const l_real c = l_real(1e152);
	int stagsave = stagprec,
 	stagmax = 30; // l_imath.cpp: sqrt() uses stagmax = 30;
 	if (stagprec>stagmax) stagprec = stagmax;
	
 	l_cinterval res;
 	l_interval absz(abs(z));
 	l_real Inf_absz(Inf(absz));
	
	if (Inf_absz > c)
 	{
	 	absz = 1 / l_interval(Inf_absz);
	 	Inf_absz = Sup(absz);
	 	res = l_cinterval( l_interval(-Inf_absz,Inf_absz),
						 	    l_interval(-Inf_absz,Inf_absz) );
		// res is the correcture interval, i.e.
		// z + res or -z + res is the inclusionof sqrt{1+z^2}
	 	res = (Inf(Re(z))>=0)? z + res : -z + res;
 	}
 	else 
 	{
	 	res = l_cinterval( l_interval(0), l_interval(1) ); // res = i
	 	if ( Sup(abs(z-res))<0.5 || Sup(abs(z+res))<0.5 )
	 	{
		 	res = l_cinterval(-Im(z),Re(z)); // Res = i*z;
		    // (1 - i*z)*(1 + i*z) = 1+z^2;
		 	res = sqrt( (1-res)*(1+res) );
	 	}
	 	else
		 	res = sqrt(1+sqr(z));
 	}
 	if (Inf(Re(res))<0)
	 	res = l_cinterval( l_interval(l_real(0),Sup(Re(res))) , Im(res) );
 	stagprec = stagsave;
 	res = adjust(res);	
 	return res;
}
// -- end sqrt1px2 ---------------------------------------------------------

l_cinterval sqrtx2m1(const l_cinterval& z) throw()
// sqrt(z^2-1);
// Calculates nearly optimal inclusions for not too wide intervals z.
// Blomquist, 04.12.2008;
{
	const l_real c = l_real(1e152);
	int stagsave = stagprec,
 	stagmax = 30; // l_imath.cpp: sqrt() uses stagmax = 30;
 	if (stagprec>stagmax) stagprec = stagmax;
	
 	l_cinterval res,u;
 	l_interval absz(abs(z));
 	l_real Inf_absz(Inf(absz));
	
 	if (Inf_absz > c)
 	{
	 	absz = 1 / l_interval(Inf_absz);
	 	Inf_absz = Sup(absz);
	 	res = l_cinterval(l_interval(-Inf_absz,Inf_absz),
							   l_interval(-Inf_absz,Inf_absz)); // res = Delta
	  // res is the correcture interval, i.e.
	   res = (Inf(Re(z))>=0)? z + res : -z + res;
 	}
 	else 
 	{
	 	res = z-1;  u = z+1;
	 	res = (Sup(abs(res))<0.5 || Sup(abs(u))<0.5)? sqrt(res*u) : sqrt(sqr(z)-1);
 	}
	
 	if (Inf(Re(res))<0)
	 	res = l_cinterval(l_interval(l_real(0),Sup(Re(res))),Im(res));
 
 	stagprec = stagsave;
 	res = adjust(res);	
 	return res;
}

l_cinterval sqrt1mx2(const l_cinterval& z) throw()
// sqrt(1-z^2);
// Calculates nearly optimal inclusions for not too wide intervals z.
// Blomquist, 04.12.2008;
{
	const l_real c = l_real(1e152);
	int stagsave = stagprec,
 	stagmax = 30; // l_imath.cpp: sqrt() uses stagmax = 30;
 	if (stagprec>stagmax) stagprec = stagmax;
	
 	l_cinterval res,u;
 	l_interval absz(abs(z));
 	l_real Inf_absz(Inf(absz));
	
 	if (Inf_absz > c)
 	{
	 	absz = 1 / l_interval(Inf_absz);
	 	Inf_absz = Sup(absz);
	 	res = l_cinterval(l_interval(-Inf_absz,Inf_absz),
							   l_interval(-Inf_absz,Inf_absz)); // res = Delta
	 	u = l_cinterval(-Im(z),Re(z)); // u = i*z;
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
	 	res = l_cinterval( l_interval(l_real(0),Sup(Re(res))) , Im(res) );
 	stagprec = stagsave;
 	res = adjust(res);
 	return res;
}


//-- sqrt_all ------------------------------------------------------- 040621 --
//
//  sqrt_all(z) computes a list of 2 intervals containing all square roots of z
//
std::list<l_cinterval> sqrt_all( const l_cinterval& z )
{
    l_real
	irez = Inf( Re(z) ),
	srez = Sup( Re(z) ),
	iimz = Inf( Im(z) ),
	simz = Sup( Im(z) );
    l_interval
	hxl( irez ), hxu( srez ), hyl( iimz ), hyu( simz );
    l_real
	resxl,  resxu, resyl, resyu;
    l_cinterval w;

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
		resyu = ( srez > 0.0 ? l_real(0.0) : -Inf( sqrt( -hxu ) ) );
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
	w = l_cinterval( l_interval(resxl,resxu), l_interval(resyl,resyu ) );
    }
    else
    //  sqrt( z ) is well-defined
	w = sqrt( z );

    std::list<l_cinterval> res;
    res.push_back( w );
    res.push_back( -w );

    return res;
}
//
//-- end sqrt_all -------------------------------------------------------------

//-- sqrt(z,n): analytic n-th root ---------------------------------- 040624 --
//
l_interval Re_Sqrt_point( const l_interval& rez, const l_interval& imz,
			int n ) // before: unsigned int n  ---------- 040624 --
//
//  Real part of analytic n-th root.
//  Do not use this as a general function
//  - it's only a subroutine for sqrt(z,n) and sqrt_all(z,n).
//  The calculation is validated but largely overestimating
//  if (rez+I*imz) is not a complex number.
//
{
    l_interval abs_z_2 = sqr( rez ) + sqr( imz );
    if( Sup( abs_z_2 ) == 0.0 )
    //  z == 0
	return l_interval(0);
    else
	return sqrt( abs_z_2, 2 * n ) * 
	    cos( Arg( l_cinterval( rez, imz ) ) / n );
}

l_interval Im_Sqrt_point( const l_interval& rez, const l_interval& imz,
			int n ) // before: unsigned int n  --- 040624 --
//
//  Imaginary part of analytic n-th root.
//  Do not use this as a general function
//  - it's only a subroutine for sqrt(z,n) and sqrt_all(z,n).
//  The calculation is validated but largely overestimating
//  if (rez+I*imz) is not a complex number.
//
{
    l_interval abs_z_2 = sqr( rez ) + sqr( imz );
    if( Sup( abs_z_2 ) == 0.0 )
    //  z == 0
	return l_interval(0);
    else
	return sqrt( abs_z_2, 2 * n ) * 
	    sin( Arg( l_cinterval( rez, imz ) ) / n );
}

l_cinterval sqrt( const l_cinterval& z, int n ) throw() // ----- 040915 --
//
//  Analytic n-th root function
//  sqrt(z,n) is undefined if z contains negative real numbers.
//
{
    if( n == 0 )
	return l_cinterval(l_interval(1));
    if( n == 1 )
	return z;
    if( n == 2 )
	return sqrt( z );
    else
    {
	l_real
	    irez = Inf( Re(z) ),
	    srez = Sup( Re(z) ),
	    iimz = Inf( Im(z) ),
	    simz = Sup( Im(z) );
	l_interval
	    hxl( irez ), hxu( srez ), hyl( iimz ), hyu( simz );
	l_real
	    resxl,  resxu, resyl, resyu;

	if( irez < 0.0 && iimz <= 0.0 && simz >= 0.0 )
        {
	  //  z contains negative real values
	    cxscthrow(STD_FKT_OUT_OF_DEF("l_cinterval sqrt(const l_cinterval& z, int n ); z contains negative real values."));
	    return z;
        }
	else
	{
	    if( simz < 0.0 )
	    {
	      //  z in lower half plane
		l_cinterval hres = sqrt( l_cinterval( Re(z), -Im(z) ), n );
		return l_cinterval( Re(hres), -Im(hres) );
	    }
	    else
	    {
		if( iimz > 0.0 ) 
		{
		  //  z in upper half plane
		    l_interval tangle = tan( ( Pi_l_interval() * n ) 
					     / ( 2 * ( n-1 ) ) );
		    l_real tanglel = Inf( tangle ), 
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
			    resxl = Inf( Re_Sqrt_point( iimz / tangle , 
							hyl, n ) );
		    }
		  //
		  //  max( Re( Root( z ) ) )
		  //
		    if ( irez >= 0.0 || Sup( hyu / irez ) <= tanglel )
		    //  upper boundary right of phi = n*Pi/(2n-2)
		    //  max( Re( Root( z ) ) ) in upper right corner
			resxu = Sup( Re_Sqrt_point( l_interval(srez), 
						    l_interval(simz), n ) );
		    else
		    {
			if ( srez < 0.0 && Inf( hyu / srez ) >= tangleu )
			//  upper boundary left of phi = n*Pi/(2n-2)
			//  max( Re( Root( z ) ) ) in upper left corner
			    resxu = Sup( Re_Sqrt_point( hxl, hyu, n ) );
			else
			//  upper boundary intersects phi = n*Pi/(2n-2)
			//  max( Re(Root( z )) ) on upper left or right corner
			    resxu = max( Sup( Re_Sqrt_point( hxl, hyu, n ) ), 
					 Sup( Re_Sqrt_point( hxu, hyu, n ) ) );
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
			  resyl = Inf(Im_Sqrt_point( hxu, srez * tangle, n ));
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
			//  max( Im(Root( z )) ) on lower or upper left corner
			  resyu = max( Sup( Im_Sqrt_point( hxl, hyl, n ) ), 
				       Sup( Im_Sqrt_point( hxl, hyu, n ) ) );
		  }
		}
		else
		{
		  //  z intersects positive real axis
		  //  min( Re( Root( z ) ) ) on real line
		  //  max( Re( Root( z ) ) ) in lower or upper right corner
		    resxl = ( irez == 0.0 ? l_real(0.0) : 
			      Inf( sqrt( hxl, n ) ) );
		    resxu = ( - iimz > simz ? 
			      Sup( Re_Sqrt_point( hxu, hyl, n ) ) : 
			      Sup( Re_Sqrt_point( hxu, hyu, n ) ) );
		  //  min( Im ( sqrt( z ) ) ) in lower left corner
		  //  max( Im ( sqrt( z ) ) ) in upper left corner
		    resyl = Inf( Im_Sqrt_point( hxl, hyl, n ) );
		    resyu = Sup( Im_Sqrt_point( hxl, hyu, n ) );
		}
		return l_cinterval( l_interval( resxl, resxu ), 
				    l_interval( resyl, resyu ) );
	    }
	}
    }
}
//
//-- end sqrt -----------------------------------------------------------------

//-- sqrt_all ------------------------------------------------------- 040628 --
//
std::list<l_cinterval> sqrt_all( const l_cinterval& z, int n )
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
    std::list<l_cinterval> res;

    if( n == 0 )
    {
	res.push_back( l_cinterval( l_interval(1), l_interval(0) ) );
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
	l_interval
	    arg_z = arg( z ), root_abs_z = sqrt( abs( z ), n );

	for(int k = 0; k < n; k++)
	{
	    l_interval arg_k = ( arg_z + 2 * k * Pi_l_interval() ) / n;

	    res.push_back( l_cinterval( root_abs_z * cos( arg_k ),
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
l_cinterval power_fast( const l_cinterval& z, int n ) throw()
{
    if( n == 0 )
	return l_cinterval( l_interval(1) );
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
	l_interval abs_z = abs(z);

	if( n < 0 && Inf( abs_z ) == 0.0 )
	//  z contains 0
	    cxscthrow (STD_FKT_OUT_OF_DEF("l_cinterval power_fast(const l_cinterval& z, int n ); z contains 0."));
	if( Sup( abs_z ) == 0.0 )
	    return l_cinterval( l_interval(0));
	else
        {
	    l_interval arg_z = arg( z );
	    l_interval abs_z_n = exp( n * ln( abs_z ) );

	    return l_cinterval( abs_z_n * cos( n * arg_z ),
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
l_cinterval power_point( const l_cinterval& z, int n ) //----------- 040715 --
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
	return l_cinterval( l_interval(0) );
    else
    {
	l_interval ln_absz = ln_sqrtx2y2(Re(z),Im(z));
	l_interval arg_z = arg(z);
	l_interval abs_z_n = exp( n * ln_absz );

	return l_cinterval( abs_z_n * cos( n * arg_z ),
			    abs_z_n * sin( n * arg_z ) );
    }
}

void update_res( const l_cinterval& res, //------------------------- 040716 --
                 l_real& resxl, l_real& resxu, l_real& resyl, l_real& resyu )
//  Subroutine of power(z,n).
//  Update extremal values of power function.
{
    resxl = min( resxl, Inf( Re(res) ) );
    resxu = max( resxu, Sup( Re(res) ) );
    resyl = min( resyl, Inf( Im(res) ) );
    resyu = max( resyu, Sup( Im(res) ) );
}

int trunc(const l_real& x)
{
    l_real y(x);  y += 0;
    bool neg(y[1]<0);
    if (neg) y = -y; // y >= 0;
    real r(y[1]);
    int res = int( _double(r) );
    y = y - res;
    if (y[1]<0.0 && res>0) res -= 1;
    if (neg) res = -res;
    return res;
}

void horizontal_check( //-------------------------------------------- 040720 --
     const l_interval& hy, l_real hyy, const l_interval& arg_h,
     l_real irez, l_real srez,
     l_real& resxl, l_real& resxu, l_real& resyl, l_real& resyu, int n )
//  Subroutine of power(z,n).
//  Check all relevant ray intersections on a horizontal boundary.
{
    int r_int, il1, il2, iu1, iu2;
    l_cinterval res;
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
    l_real arg_hlR = Inf( 2 * nofrays * arg_h / Pi_l_interval() );
    if( arg_hlR[1] >= 0.0 )
	r_int = trunc(arg_hlR);
    else
	r_int = trunc(arg_hlR - 1.0); 
    il1 = r_int + 1;

    l_real arg_huR = Sup( 2 * nofrays * arg_h / Pi_l_interval() );
    if( arg_huR[1] >= 0.0 )
	r_int = trunc(arg_huR);
    else
	r_int = trunc(arg_huR - 1.0);

    iu1 = r_int;

    if( iu1 >= il1 )
    //  at least one ray intersection
    //  maybe more?
    {
	if( iu1 > il1 + 3 ) {
	//
	//  we're in trouble:
	//  there are more than 4 ray intersections
	//  now we must decide which of these are relevant
	//  depending on il1, iu1, n and on the quadrants,
	//  4 to 7 intersections must be considered
	//
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
		    l_interval cotangle = cot( ( Pi_l_interval() * i ) / 
					       ( 2 * nofrays ) );
		    res = power_point( l_cinterval( hy * cotangle , hy ), n );
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
	}
      //
      //  list 1 has been left for processing
      //
	for(int i = il1; i <= iu1; i++)
	{
	    l_interval cotangle = cot( ( Pi_l_interval() * i ) / 
				       ( 2 * nofrays ) );
	    res = power_point( l_cinterval( hy * cotangle , hy ), n );
	    update_res( res, resxl, resxu, resyl, resyu );
	}
    }
}

void vertical_check( //---------------------------------------------- 040720 --
     const l_interval& hx, l_real hxx, const l_interval& arg_h,
     l_real iimz, l_real simz,
     l_real& resxl, l_real& resxu, l_real& resyl, l_real& resyu, int n )
//  Subroutine of power(z,n).
//  Check all relevant ray intersections on a vertical boundary.
{
    int r_int, il1, il2, iu1, iu2;
    l_cinterval res;
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

    l_real arg_hlR = Inf( 2 * nofrays * arg_h / Pi_l_interval() );
    if( arg_hlR[1] >= 0.0 )
	r_int = trunc(arg_hlR);
    else
	r_int = trunc(arg_hlR - 1.0);
    il1 = r_int + 1;

    l_real arg_huR = Sup( 2 * nofrays * arg_h / Pi_l_interval() );
    if( arg_huR[1] >= 0.0 )
	r_int = trunc(arg_huR);
    else
	r_int = trunc(arg_huR - 1.0);
    iu1 = r_int;

    if( iu1 >= il1 )
    //  at least one ray intersection
    //  maybe more?
    {
	if( iu1 > il1 + 3 ) {
	//
	//  we're in trouble:
	//  there are more than 4 ray intersections
	//  now we must decide which of these are relevant
	//  depending on il1, iu1, n and on the quadrants,
	//  4 to 7 intersections must be considered
	//
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
			l_interval tangle = tan( ( Pi_l_interval() * i ) / 
						 ( 2 * nofrays ) );
			res = power_point( l_cinterval( hx, hx * tangle ), n );
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
	  }  
      //
      //  list 1 has been left for processing
      //
	for(int i = il1; i <= iu1; i++)
	{
	    l_interval tangle = tan( ( Pi_l_interval() * i ) / 
				     ( 2 * nofrays ) );
	    res = power_point( l_cinterval( hx, hx * tangle ), n );
	    update_res( res, resxl, resxu, resyl, resyu );
	}
    }
}

l_cinterval power( const l_cinterval& z, int n ) throw() //---- 040720 --
//
//  Power function for integer powers with optimal (save roundoff) accuracy.
//
{
    if( n == 0 )
	return l_cinterval( l_interval(1) );
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
	l_interval abs_z = abs(z);

	if( n < 0 && Inf( abs_z ) == 0.0 )
	//  z contains 0
	    cxscthrow (STD_FKT_OUT_OF_DEF("l_cinterval power(const l_cinterval& z, int n ); z contains 0."));
	if( Sup( abs_z ) == 0.0 )
	    return l_cinterval( l_interval(0), l_interval(0) );
	else
        {
	    l_real
		irez = Inf(Re(z)),
		srez = Sup(Re(z)),
		iimz = Inf(Im(z)),
		simz = Sup(Im(z));
	    l_interval
		hxl( irez ), hxu( srez ), hyl( iimz ), hyu( simz );
	    l_real
		resxl,  resxu, resyl, resyu;
	    l_cinterval res;
	  //
	  //  extremal values lie in corner points or on intersections of the
	  //  boundary of z with any ray  arg( w ) = k*Pi/(2*(n-1))
	  //
	  //  first check the corners of z
	  //
	    res = power_point( l_cinterval( hxl, hyl ), n );
	    resxl = Inf( Re(res) );
	    resxu = Sup( Re(res) );
	    resyl = Inf( Im(res) );
	    resyu = Sup( Im(res) );
	    res = power_point( l_cinterval( hxu, hyl ), n );
	    update_res( res, resxl, resxu, resyl, resyu );
	    res = power_point( l_cinterval( hxl, hyu ), n );
	    update_res( res, resxl, resxu, resyl, resyu );
	    res = power_point( l_cinterval( hxu, hyu ), n );
	    update_res( res, resxl, resxu, resyl, resyu );
	  //
	  //  consider the origin, if it is a boundary point
	  //  (negative n have been taken care of above)
	  //
          //	    if( irez * srez <= 0.0 and iimz * simz <= 0.0 
          //		and irez * srez * iimz * simz == 0.0 )
		if ( 0<=Re(z) && 0<=Im(z) && 
		     (irez==0 || srez==0 || iimz==0 || simz==0 ) )
		update_res( l_cinterval( l_interval(0), l_interval(0)), 
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
		l_interval arg_h = arg( l_cinterval( Re(z), hyl ) );

	      //  check if there are ray intersections
		horizontal_check( hyl, iimz, arg_h, irez, srez,
				  resxl, resxu, resyl, resyu, n );
	    }
	  //  Ib: upper boundary
	    if( simz != 0.0 )
	    {
		l_interval arg_h = arg( l_cinterval( Re(z), hyu ) );

	      //  check if there are ray intersections
		horizontal_check( hyu, simz, arg_h, irez, srez,
				  resxl, resxu, resyl, resyu, n );
	    }
	  //  IIa: left boundary
	    if( irez != 0.0 )
	    {
		l_interval arg_h = arg( l_cinterval( hxl, Im(z) ) );

	      //  check if there are ray intersections
		vertical_check( hxl, irez, arg_h, iimz, simz,
				resxl, resxu, resyl, resyu, n );
	    }
	  //  IIb: right boundary
	    if( srez != 0.0 )
	    {
		l_interval arg_h = arg( l_cinterval( hxu, Im(z) ) );

	      //  check if there are ray intersections
		vertical_check( hxu, srez, arg_h, iimz, simz,
				resxl, resxu, resyl, resyu, n );
	    }
	    return l_cinterval( l_interval( resxl, resxu ), 
				l_interval( resyl, resyu ) );
	}
    }
}
//
//-- end power ----------------------------------------------------------------


//-- pow ------------------------------------------------------------ 040627 --
//
//  Analytic power function for real interval exponent, based on Ln.
//

l_cinterval pow( const l_cinterval& z, const l_interval& p ) throw()
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

l_cinterval pow( const l_cinterval& z, const l_cinterval& p ) throw()
{
    return exp( p*Ln(z) );
}

//
//-- end pow ------------------------------------------------------------------

//-- pow_all -------------------------------------------------------- 041013 --
//
//  Non-analytic power function for real l_interval exponent.
//
//  If 0 \not\in z, then compute four rectangular intervals that comprehend
//  an annulus which contains all values  zeta^phi, zeta in z, phi in p.
//
//  If 0 in z and negative reals in p, then abort execution
//  (potential modification: return the entire complex plane).
//
std::list<l_cinterval> pow_all( const l_cinterval& z, 
				const l_interval& p ) throw()
{
    l_interval abs_z = abs(z);

    if( 0.0 < Inf( abs_z ) )
    {
	l_interval abs_z_p = exp( p * ln( abs_z ) );

      //  Inner and outer radii of the annulus are inf/sup( abs_z_n )
      //  Inscribed square has side length sqrt( 2 ) * rad_1
	l_interval rad_1 = Sqrt2r_l_interval() * Inf( abs_z_p );
	l_interval rad_2 = l_interval(Sup( abs_z_p ));

	std::list<l_cinterval> res;

      //  4 intervals covering the annulus, counter-clockwise
	res.push_back(l_cinterval(l_interval( Inf( rad_1 ), Sup( rad_2 ) ),
				  l_interval( -Sup( rad_1 ), Sup( rad_2 ) )));
	res.push_back(l_cinterval(l_interval( -Sup( rad_2 ), Sup( rad_1 ) ),
				  l_interval( Inf( rad_1 ), Sup( rad_2 ) ) ));
	res.push_back(l_cinterval(l_interval( -Sup( rad_2 ), -Inf( rad_1 ) ),
				  l_interval( -Sup( rad_2 ), Sup(rad_1) ) ) );
	res.push_back(l_cinterval(l_interval( -Sup( rad_1 ), Sup( rad_2 ) ),
				  l_interval( -Sup( rad_2 ), -Inf(rad_1) ) ) );

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
	    l_interval abs_z_p = exp( p * ln( l_interval( Sup( abs_z ) ) ) );
	    l_real rad_p = Sup( abs_z_p );

	    std::list<l_cinterval> res;

	    res.push_back( l_cinterval( l_interval( -rad_p, rad_p ),
					l_interval( -rad_p, rad_p ) ) );

	    return res;
	}
	else
        {
	//
	//  The set   zeta^phi, zeta in z, phi in p   is unbounded
	//  if inf( p ) < 0.  0^p is undefined for p <= 0.
	//
	    cxscthrow(STD_FKT_OUT_OF_DEF("pow_all(l_cinterval, l_interval); 0^p is undefined for p <= 0."));
	    std::list<l_cinterval> res;
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
l_interval f_aux_asin( const l_interval& x, const l_interval& y ) 
//
//  auxiliary function for evaluating the imaginary part of asin(z);
//  f_aux_asin(z) = ( |z+1| + |z-1| ) / 2;  z = x + i*y;
//                = ( sqrt[(x+1)^2 + y^2] + sqrt[(x-1)^2 + y^2] )/2; 
//  Blomquist, 23.01.07;        
{
    l_interval res;
//  interval sqr_y = sqr( y );
//  interval res = ( sqrt( sqr( x + 1.0 ) + sqr_y ) + 
//                   sqrt( sqr( x - 1.0 ) + sqr_y ) ) / 2.0;
    res = abs(x);
    if (y != 0.0 || Inf(res) < 1.0) 
    { 
	res = sqrtx2y2(x+1.0,y) + sqrtx2y2(x-1.0,y);  // Blomquist;
	times2pown(res,-1);
    }

    if ( Sup(res)==Infinity ) // Blomquist: program stop, if overflow occurs.
	cxscthrow (STD_FKT_OUT_OF_DEF(
              "l_cinterval asin( const l_cinterval& z); z out of range"));

  //  Now we correct a possible overestimation of the lower bound
  //  of res.
  //  It holds: 
  //  (sqrt[(x+1)^2 + y^2] + sqrt[(x-1)^2 + y^2])/2 >= Max(1,|x|);
    l_real hlb = max( l_real(1.0), abs( Sup( x ) ) );
    if( Inf( res ) < hlb ) //  this is an invalid overestimation!
	res = l_interval( hlb, Sup(res) );
    return res;
}

l_interval ACOSH_f_aux( const l_interval& x, const l_interval& y )
// Hilfsfunktion zur Berechnung des Imagin�rteils von asin(z),
// z = x + i*y;
// Bezeichnungen:
// f_aux_asin(x,y) = alpha(x,y);
// alpha(x,y) := ( sqrt[(x+1)^2 + y^2] + sqrt[(x-1)^2 + y^2] )/2;
// R�ckgabewert:
// 1.    res = acosh( f_aux_asin(x,y) ), falls  |x|>2 oder |y|>2;
// 2.    res = ln[ alpha + sqrt(alpha^2 - 1) ]; d = alpha(x,y) -1
//           = ln[1 + sqrt(d)*(sqrt(d) + sqrt(2+d))], falls gilt:
//       |x|<=2 und |y|<=2;   x,y sind nur Punktintervalle!       
// Blomquist, 03.02.2007;

{
    const int c  = 1017,
              c1 = 1022;
    l_interval res,xa(abs(x)),ya(abs(y)),diff, diff1,S1,S2;
    l_real rx((Inf(xa))), ry(Inf(ya));
    int d,p,ex,ex1,ex2,s(0);

    if (rx>2.0 || ry>2.0) res = acosh( f_aux_asin(x,y) );
    else 
    { //  Jetzt mit Verbesserungen:!
	if (rx == 1.0) 
	{
	    S1 = 0.5 + (ya/4) / (sqrt1px2(ya/2) + 1);
	    S2 = sqrt(ya) * sqrt(S1);  // S2 = sqrt(delta);
	    res = lnp1(S2 * (S2 + sqrt( 2+sqr(S2) ))); 
	}
	else 
	    if (rx<1.0) // 0 <= x < +1;
	    {
		diff = 1 - xa;
		d  = expo_gr(ry);
		ex = expo_gr(diff);
		if (d > ex) ex = d;
		ex = -ex; 
		// ex = -Max( expo(ry[1]),expo(diff[1]) ) > 0;
		if (ex>LI_maxexpo1)
		{
		    times2pown(diff,LI_maxexpo1);
		    times2pown(diff,ex-LI_maxexpo1);
		    times2pown(ya,LI_maxexpo1);
		    times2pown(ya,ex-LI_maxexpo1);
		}
		else
		{
		    times2pown(diff,ex);
		    times2pown(ya,ex);
		}
		S2 = sqrtx2y2(diff,ya) + diff; // S2 '=' 1;
		if (ex>=LI_maxexpo1) 
		{
		    d = ex - c;  // c = 1017;
		    if (d%2 != 0) d += 1;
		} else d = 0;
		S2 = comp(0.5,ex-d+1) / S2;
		diff = 1+xa;
		S1 = comp(0.5,-d+1) / (sqrtx2y2(diff,y) + diff);
		S1 = S1 + S2;
		times2pown(S1,-1);  // S1 = { ... } / 2; 
		S1 = sqrt(S1); 
		d = d/2;
		if (d != 0) times2pown(S1,d);
		S1 = abs(y) * S1; // S1 = sqrt(alpha - 1)
		res = lnp1(S1 * (S1 + sqrt( 2+sqr(S1) ))); 
	    }
	    else  // 1 < x <= 2;
		if (y == 0)
		{
		    S1 = sqrt(xa-1);
		    res = lnp1( S1*(S1 + sqrt(xa+1)) );
		}
		else // 1 < x <= 2  and  0 < |y| <= 2
		{
		    diff = xa - 1;
		    ex1 = expo_gr(diff);
		    ex2 = expo_gr(ry);
		    ex = ex1;
		    if (ex2>ex) ex = ex2;
		    ex = -ex; // ex = -Max{expo(x-1),expo(y)}
		    p = -2*ex2 + 1 + ex1;
		    if (ex>p) p = ex; // p = Max{-2*ex2 + 1 + ex1,ex}
		    if (p>c1) 
		    {
			d = p - c;  // c = 1017;
			if (d%2 != 0) d += 1;
		    } else d = 0;
		    S1 = xa + 1;
		    if (d>c1) 
		    { 
 			S2 = l_interval(MinReal); 
			s = d - c1;  // s > 0;
			times2pown(S1,s);
			times2pown(ya,s);
		    }
		    else S2 = l_interval(comp(0.5,-d+1));
		    S1 = S2 / (sqrtx2y2(S1,ya) + S1); // S1 = u * 2^-d;

		    diff = xa - 1;  diff1 = diff;
		    S2 = l_interval(1);
		    Times2pown(S2,ex-d);
		    ya = abs(y);
		    if (ex>LI_maxexpo1)
		    {
			times2pown(diff,LI_maxexpo1);
			times2pown(diff,ex-LI_maxexpo1);
			times2pown(ya,LI_maxexpo1);
			times2pown(ya,ex-LI_maxexpo1);
		    } else
		    {
			times2pown(diff,ex);
			times2pown(ya,ex);
		    }
		    S2 = S2 / (sqrtx2y2(diff,ya) + diff);
                    // S2 = v * 2^-d;

		    s = -2*ex2-d+1;
		    Times2pown(diff1,s);
		    ya = abs(y);
		    Times2pown(ya,-ex2);
		    diff1 = diff1 / (ya*ya);
                    // diff1 = w * 2^-d;

		    res = diff1 + S2;
		    res = sqrt((res + S1)/2);
		    s = expo_gr(res);
		    if (d/2+s>c1)
		    { // Overflow vermeiden durch:
			ex = c1 - s;
			times2pown(res,ex);
			ya = abs(y);
			times2pown(ya,d/2-ex);
			S1 = ya * res;
		    } else
		    {
			times2pown(res,d/2);
			S1 = abs(y) * res;
		    }
		    res = lnp1(S1 * (S1 + sqrt( 2 + sqr(S1) )));	
		}
    }
	return res;
} // ACOSH_f_aux

l_interval Asin_beta( const l_interval& x, const l_interval& y )
// Calculating the improved real part of asin(z); Blomquist 22.01.2007;
// Re(asin(z)) = asin[ 2x/(sqrt[(x+1)^2+y^2] + sqrt[(x-1)^2+y^2]) ]=asin[beta];
// Improvements for beta --> +1  and  beta --> -1 are necessary because of 
// nearly vertical tangents of the real asin(t)-function for |t|-->+1;
// This function should only be used for POINT intervals x,y!  z = x + i*y;   
{
    const real c1 = 0.75;
    bool neg_x;
    l_real Infxa;
    int ex_d,ex_y,ex;
    l_interval res,beta,abs_beta,delta,w_delta,xa,Ne,Sqrt_Ne,Kla;
    Ne = sqrtx2y2(1+x,y) + sqrtx2y2(1-x,y);
    beta = x / ( Ne/2 );
    if (Inf(beta)<-1) Inf(beta) = -1;
    if (Sup(beta)> 1) Sup(beta) = 1; 
    abs_beta = abs(beta);
    if (Inf(abs_beta) < c1) res = asin(beta); // Normal calculation 
    else { // Inf(abs_beta)>=c1; Calculation now with improvements:
	Ne = Ne + 2.0; // Ne = 2 + sqrt[(1+x)^2 + y^2] + sqrt[(1-x)^2 + y^2];
	Sqrt_Ne = sqrt(Ne); 
	xa = x;
	neg_x = Inf(x)<0;
	if (neg_x) xa = -xa; // Inf(xa) >0 :
	Infxa = Inf(xa);

	if (Infxa > 1.0) 
	{ // x > +1;
	    if (y == 0.0) 
	    {
		w_delta = 0.0;
		delta = 0.0;
	    } 
	    else
	    {
		beta = xa - 1; // beta > 0;
		Infxa = Inf(y);
		ex_y = expo_gr(Infxa);
		ex_d = expo_gr(beta);
		ex = ex_y-ex_d-50;
		if (ex > 0) 
		{ 
		    times2pown(beta,ex);
		    res = l_interval(comp(0.5,-ex+1));
		    delta = abs(y)/beta;
		    Kla = sqrtx2y2(res,delta) + res;
		    w_delta = sqrt(2*beta)*delta/(Sqrt_Ne * sqrt(Kla));
		} 
		else
		{
		    delta = abs(y)/beta;
		    Kla = sqrt1px2(delta) + 1;
		    w_delta = sqrt(2*beta)*delta / (Sqrt_Ne * sqrt(Kla));
		}
		delta = sqr(w_delta);
	    }
	} 
	else
	    if (Infxa == 1.0) 
	    { // x = 1;
		delta = 2*abs(y) / Ne;
		w_delta = sqrt(2*abs(y)) / Sqrt_Ne;
	    } 
	    else 
	    { // 0.75 <= x < 1;
		if (y == 0.0) 
		{
		    delta = 1 - xa;  // delta is a point interval!
		    w_delta = sqrt(delta);
		} 
		else // 0.75 <= x < 1  and  y != 0:
		{
		    beta = 1 - xa; // beta > 0;
		    Infxa = Inf(y);
		    ex_y = expo_gr(Infxa);
		    ex_d = expo_gr(beta);
		    ex = ex_y-ex_d-50;  
		    if (ex > 0) 
		    { 
			times2pown(beta,ex);
			res = l_interval(comp(0.5,-ex+1));
			Kla = sqrtx2y2(res,y/beta) + res;
		    } 
		    else // Now without scaling
			Kla = sqrt1px2(y/beta) + 1;
		    times2pown(beta,1);  // beta * 2;
		    delta = beta * Kla / Ne;
		    w_delta = sqrt(beta) * sqrt(Kla) / Sqrt_Ne;
		}
	    }
	res = Pid2_l_interval() - asin( w_delta*sqrt(2-delta) );
	if (neg_x) res = -res;
    }
    return res;
} // Asin_beta(...)

l_cinterval asin( const l_cinterval& z ) throw() //------------- 040730 --
{
    const real gr = 6.355804e307; // upper bound for abs(rez),abs(imz)
    l_interval
	rez = Re(z),
	imz = Im(z);

    l_real irez = Inf(rez),
	 srez = Sup(rez),
	 iimz = Inf(imz),
	 simz = Sup(imz);

    l_interval hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

    l_real resxl, resxu, resyl, resyu;

    bool bl    = iimz< 0.0 && simz>0.0,
         raxis = iimz==0.0 && simz==0;

  //
  //  1st: check for singularities
  //
    if( (irez<-1 && (bl || (iimz<0  && simz==0))) || 
	(srez >1 && (bl || (iimz==0 && simz>0))) )
    cxscthrow(STD_FKT_OUT_OF_DEF("l_cinterval asin( const l_cinterval& z ); z contains singularities."));
  //
  //  check for too large bounds of abs(rez) and abs(imz) to prevent 
  //  overflow by calculating f_aux_asin(...)
  //
    resxl = max(abs(irez),abs(srez));  
    resxu = max(abs(iimz),abs(simz));
    if (resxl>gr || resxu>gr) 
	cxscthrow(STD_FKT_OUT_OF_DEF("l_cinterval asin( const l_cinterval& z ); z with too large bounds."));
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
	    resxl = Inf( Asin_beta(hxl,l_interval( max(-iimz,simz) )) ); // Blomquist, 19.06.2005;
	if( srez < 0.0 )
//	resxu = Sup( asin( hxu / f_aux_asin( hxu, interval( max( - iimz, simz ) ) ) ) );
	    resxu = Sup( Asin_beta(hxu,l_interval( max(-iimz,simz) )) ); // Blomquist, 19.06.2005;
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
		resyu = - Inf( ACOSH_f_aux( l_interval(0), hyu ) );
	}
	else
	//  most of z in quadrant IV
	{
	    resyl = - Sup( ACOSH_f_aux( hxu, hyl ) );
	    if( irez > 0.0 )
		resyu = - Inf( ACOSH_f_aux( hxl, hyu ) );
	    else
		resyu = - Inf( ACOSH_f_aux( l_interval(0), hyu ) );
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
	      resyl = Inf( ACOSH_f_aux( l_interval(0), hyl ) );
	}
	else
	//  most of z in quadrant I
	{
	    resyu = Sup( ACOSH_f_aux( hxu, hyu ) );
	    if( irez > 0.0 )
		resyl = Inf( ACOSH_f_aux( hxl, hyl ) );
	    else
		resyl = Inf( ACOSH_f_aux( l_interval(0), hyl ) );
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

  return l_cinterval( l_interval( resxl,resxu ), l_interval( resyl,resyu ) );

}
//
//-- end asin -----------------------------------------------------------------

l_interval Asin_arg( const l_interval& x, const l_interval& y )
// Asin_arg berechnet f�r Punktintervalle x,y mit 
// beta := 2*|x| / ( sqrt((x+1)^2+y^2) + sqrt((x-1)^2+y^2) )
// und delta := 1 - |beta|  eine Einschlie�ung von:  
//                arcsin[ sqrt(delta)*sqrt(2-delta) ].
// Das Punktintervall x darf auch negativ sein!
// Blomquist, 06.02.2007;
{
    l_interval res,Ne,Sqrt_Ne,xa,beta,delta,w_delta,Kla;
    l_real Infxa;
    int ex,ex_y,ex_d;
    bool neg_x;

    Ne = 2 + sqrtx2y2(1+x,y) + sqrtx2y2(1-x,y);
    Sqrt_Ne = sqrt(Ne); 
    xa = x;
    neg_x = Inf(x)<0;
    if (neg_x) xa = -xa; // Inf(xa) >0 :
    Infxa = Inf(xa);

    if (Infxa > 1.0) { // x > +1;
	if (y == 0.0) {
	    w_delta = 0.0;
	    delta = 0.0;
	} else
	{
	    beta = xa - 1; // beta > 0;
	    Infxa = Inf(y);
	    ex_y = expo_gr(Infxa);
	    ex_d = expo_gr(beta);
	    ex = ex_y-ex_d-50;
	    if (ex > 0) { 
		times2pown(beta,ex);
		res = l_interval(comp(0.5,-ex+1));
		delta = abs(y)/beta;
		Kla = sqrtx2y2(res,delta) + res;
		w_delta = sqrt(2*beta)*delta/(Sqrt_Ne * sqrt(Kla));
	    } else
	    {
		delta = abs(y)/beta;
		Kla = sqrt1px2(delta) + 1;
		w_delta = sqrt(2*beta)*delta / (Sqrt_Ne * sqrt(Kla));
	    }
	    delta = sqr(w_delta);
	}
    } else
	if (Infxa == 1.0) { // x = 1;
	    delta = 2*abs(y) / Ne;
	    w_delta = sqrt(2*abs(y)) / Sqrt_Ne;
	} else { // 0.75 <= x < 1;
	    if (y == 0.0){
		beta = 1 - xa; // beta is a point interval!
		delta = 4*beta / Ne;
		w_delta = 2*sqrt(beta) / Sqrt_Ne;
	    } else // 0.75 <= x < 1  and  y != 0:
	    {
		beta = 1 - xa; // beta > 0;
		Infxa = Inf(y);
		ex_y = expo_gr(Infxa);
		ex_d = expo_gr(beta);
		ex = ex_y-ex_d-50;  
		if (ex > 0) { 
		    times2pown(beta,ex);
		    res = l_interval(comp(0.5,-ex+1));
		    Kla = sqrtx2y2(res,y/beta) + res;
		} else // Now without scaling
		    Kla = sqrt1px2(y/beta) + 1;
		times2pown(beta,1);  // beta * 2;
		delta = beta * Kla / Ne;
		w_delta = sqrt(beta) * sqrt(Kla) / Sqrt_Ne;
	    }
	}
    res = asin( w_delta*sqrt(2-delta) );

    return res;
} // Asin_arg

l_interval Acos_beta( const l_interval& x, const l_interval& y )
// Calculating the improved real part of acos(z); Blomquist 05.06.2005;
// Re(acos(z)) = acos[ 2x / (sqrt[(x+1)^2+y^2] + sqrt[(x-1)^2+y^2]) ]
{
    const real c1 = 0.75;
    l_interval res(0),beta;
    beta = x / ( (sqrtx2y2(1+x,y) + sqrtx2y2(1-x,y))/2 );
    if (Inf(beta)<-1) Inf(beta)=-1;
    if (Sup(beta)> 1) Sup(beta)= 1; 

    if (Sup(beta)<c1) 
	if (Sup(beta)<-c1) { // Improvement for beta --> -1
	    res = Pi_l_interval() - Asin_arg(x,y);
	} else res = acos(beta); // Normal calculation
    else  // Sup(beta)>=c1
	res = Asin_arg(x,y);
    return res;
} 


//-- acos ----------------------------------------------------------- 040730 --
//
// l_cinterval acos( const l_cinterval& z )
// {
//   w := acos(z);
//   Re(w) in a new Version,
//   Im(w) = -Im(asin(z));  Blomquist, 14.06.2005;
// }
//
//-- acos -----------------------------------------------------------------


//-- acos ---------------------------------------------------
//
l_cinterval acos( const l_cinterval& z ) throw()
{
    const real gr = 6.355804e307; // upper bound for abs(rez),abs(imz)
    l_interval 
	rez = Re(z),
	imz = Im(z);

    l_real
	irez = Inf(rez),
	srez = Sup(rez),
	iimz = Inf(imz),
	simz = Sup(imz);

    l_interval
	hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

    bool bl    = iimz< 0.0 && simz>0.0,
         raxis = iimz==0.0 && simz==0;
    l_real
	resxl, resxu, resyl, resyu;
  //
  //  1st: check for singularities
  //
    if( (irez<-1 && (bl || (iimz<0  && simz==0))) || 
	(srez >1 && (bl || (iimz==0 && simz>0))) )
	cxscthrow(STD_FKT_OUT_OF_DEF("l_cinterval acos( const l_cinterval& z ); z contains singularities."));
  //
  //  check for too large bounds of abs(rez) and abs(imz) to prevent 
  //  overflow by calculating f_aux_asin(...)
  //
    resxl = max(abs(irez),abs(srez));  
    resxu = max(abs(iimz),abs(simz));
    if (resxl>gr || resxu>gr) 
	cxscthrow(STD_FKT_OUT_OF_DEF("l_cinterval acos( const l_cinterval& z ); z with too large bounds."));
  //
  //  2nd: real part
  //
  //  Blomquist, 05.02.2007;
      if( iimz < 0.0 && simz > 0.0 )
      //  z intersects [-1,1] on the x-axis
      {
	  if( irez <= 0.0 ) resxu = Sup( acos( hxl ) );
	  else resxu = Sup( Acos_beta(hxl,l_interval( max(-iimz,simz) )) );

	  if( srez < 0.0 ) 
	       resxl = Inf( Acos_beta(hxu,l_interval(max(-iimz,simz))) );
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
	  resyl = -Sup( ACOSH_f_aux( hxl, hyl ) );
	  if( srez < 0.0 )
	    resyu = -Inf( ACOSH_f_aux( hxu, hyu ) );
	  else
	    resyu = -Inf( ACOSH_f_aux( l_interval(0), hyu ) );
	}
      else
	//  most of z in quadrant IV
	{
	  resyl = -Sup( ACOSH_f_aux( hxu, hyl ) );
	  if( irez > 0.0 )
	    resyu = -Inf( ACOSH_f_aux( hxl, hyu ) );
	  else
	    resyu = -Inf( ACOSH_f_aux( l_interval(0), hyu ) );
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
	    resyl = Inf( ACOSH_f_aux( l_interval(0), hyl ) );
	}
      else
	//  most of z in quadrant I
	{
	  resyu = Sup( ACOSH_f_aux( hxu, hyu ) );
	  if( irez > 0.0 )
	     resyl = Inf( ACOSH_f_aux( hxl, hyl ) );
	  else
	    resyl = Inf( ACOSH_f_aux( l_interval(0), hyl ) );
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
	  resyl = -Sup( ACOSH_f_aux( hxl, hyl ) );
	  resyu = Sup( ACOSH_f_aux( hxl, hyu ) );
	}
      else
	{
	  resyl = -Sup( ACOSH_f_aux( hxu, hyl ) );
	  resyu =  Sup( ACOSH_f_aux( hxu, hyu ) );
	}
    }

  return l_cinterval( l_interval( resxl, resxu ),-l_interval( resyl, resyu ) );

}
//
//-- end acos -----------------------------------------------------------------


//-- asinh ---------------------------------------------------------- 040730 --
//
l_cinterval asinh( const l_cinterval& z ) throw()
//
//  asinh( z ) = i * asin( -i * z )
//
{
  l_cinterval res = asin( l_cinterval( Im(z), -Re(z) ) );
  return l_cinterval( -Im(res), Re(res) );
}
//
//-- end asinh ----------------------------------------------------------------


//-- acosh ---------------------------------------------------------- 040908 --
//
l_cinterval acosh( const l_cinterval& z ) throw()
//
//  acosh( z ) = i * acos( z ) = +/- i * ( pi / 2 - asin( z ) )
//
{
    l_interval
	rez = Re(z),
	imz = Im(z);

    l_real
	irez = Inf(rez),
	srez = Sup(rez),
	iimz = Inf(imz),
	simz = Sup(imz);

    l_interval
	hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

    l_real
	resxl, resxu, resyl, resyu;

  // cinterval res;

  //
  //  1st: check for singularities
  //
    if( ( iimz <= 0.0 && simz >= 0.0 ) && ( irez < 1.0 ) )
	cxscthrow(STD_FKT_OUT_OF_DEF("l_cinterval acosh( const l_cinterval& z ); z contains singularities."));
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
	l_cinterval res = acos(z);
	return l_cinterval( -Im(res),Re(res) );
    }
  //
  //  3rd: z in lower half plane
  //  acosh( z ) =  - i * ( pi / 2 - asin( z ) )
  //
    if( simz < 0.0 )
    {
//      cinterval res = HALFPI() - asin( z );
	l_cinterval res = acos(z);  // Blomquist, 14.06.2005
	return l_cinterval( Im(res), -Re(res) );
    }
  //
  //  z intersects [1,infinity)
  //
  //  real part
  //  minimum on the left on real axes, maximum in lower or upper right corner
  //
    resxl = Inf( acosh( hxl ) );
    l_interval ytilde( max( -iimz, simz ) );
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

    return l_cinterval(l_interval( resxl, resxu ),l_interval( resyl, resyu ));
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

void re_atan( const l_interval& y, l_interval& x, l_interval& res )
// Calculating an inclusion res of 1 - y^2 - x^2;
// If |1 - y^2 - x^2| is too small, then x = (2^d * x) and
// res = 2^d * (1 - y^2 - x^2); 
// x is always a point interval.
// Blomquist, 26.02.2007;
{
    if ( x==0.0 ) res = 1.0;
    else // x <> 0:
    {
	interval y1,x1;
	int ex,ex_x,d;
	l_interval y_(y);
	l_real lr;
	y1 = y;  x1 = x;
	y1 = 1 - sqr(y1) - sqr(x1);
	if ( expo(Sup(y1))<-2 )
	{
	    res = 1 - sqr(y) - sqr(x);
	    lr = Sup(abs(res));
	    ex = expo_gr(lr);
	    if (ex < -10000) res = 0;
	    else
	    {
		lr = Sup(x);
		ex_x = expo_gr(lr);
		d = 1 - ex;
		ex_x = 1022 - ex_x;
		if (ex_x<d) d = ex_x;  // d: Minimum{1-ex,1022-ex_x}
		if (d>1022) d = 1022;  // d: Minimum{1-ex,1022-ex_x,1022}
		if (d%2 != 0) d--;     // d divisible by 2;
		ex = d/2;
		times2pown(x,ex);
		times2pown(y_,ex);
		res = comp(0.5,d+1) - sqr(x) - sqr(y_);
		times2pown(x,ex);
	    }
	}
	else res = 1 - sqr(y) - sqr(x);
    }
} // re_atan

void re_vert( const l_real& x, const l_interval& hx,
	      const l_real& rew_inf, const l_real& rew_sup,
 	      l_real& resxl, l_real& resxu ) //---------------------- 040729 --
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
	l_interval hx2(hx);
	times2pown(hx2,1);
	if( x > 0.0 )
	//  w in quadrants I and/or II
	//  atan is the inverse function of tan(t), t in (-pi/2,pi/2).
	{
	    resxl = rew_sup > 0.0 ? Inf( Atan( hx2,rew_sup )/2.0 )
		    : ( rew_sup < 0.0 ? Inf( (Atan( hx2,rew_sup ) + Pi_l_interval() )/2.0 )
		                          : Inf( Pid4_l_interval() ) );

	    resxu = rew_inf > 0.0 ? Sup( Atan( hx2,rew_inf )/2.0 )
		        : ( rew_inf < 0.0 ? Sup( (Atan( hx2,rew_inf ) + Pi_l_interval())/2.0 )
			                  : Sup( Pid4_l_interval() ) );
	}
	else
	//  w in quadrants III and/or IV
	{
	    resxl = rew_inf < 0.0 ? Inf( (Atan( hx2,rew_inf ) - Pi_l_interval())/2.0 )
                            : ( rew_inf > 0.0 ? Inf( Atan( hx2,rew_inf )/2.0 )
			                      : -Sup( Pid2_l_interval()/2.0 ) );
	    resxu = rew_sup < 0.0 ? Sup( (Atan( hx2,rew_sup ) - Pi_l_interval())/2.0 )
                            : ( rew_sup > 0.0 ? Sup( Atan( hx2,rew_sup )/2.0 )
			                      : -Inf( Pid4_l_interval() ) );
	}
    }
} //  re_vert

l_interval Aux_1_atan(const l_real& x)
// x > minreal;  minreal = 2^(-1074);
// Calculating: ln[ 1+2/(sqrt(1+x^2)-1) ], [x] = x,
// [x] is a point interval !
// Tested in detail, Blomquist; 15.02.07;
{
    const int exOv = +510;
    const int exUn = 1;     

    l_interval res,
	ix(x),      // ix is point interval with x>0;
	r,t;
    int ex(expo_gr(x));

    if (ex>=exOv) 
    { // preventing overflow
	r = 1/ix;  
	res = lnp1( 2*r*( sqrt1px2(r) + r) );
    } else
	if (ex<=exUn) 
	{ // x < 2^(1)
	    res = sqr(ix);
            // Calculating  r := 2*ln(ix)  for  x --> 0;
	    if (ex<-1000) {
		Times2pown(ix,-ex);
		r = ln(ix) + ex*Ln2_l_interval();
	    } else r = ln(ix);
	    times2pown(r,1);
	    t = 1+sqrt(1+res);
	    times2pown(t,1);  // t = 2*{1+sqrt(1+x^2)}
	    res = ln(t + res) - r;
	} else 
	{ // normal calculation
	    res = sqrtp1m1( sqr(ix) ); // res = sqrt(1+x^2)-1
	    stagprec++;
	    res = 2/res;
	    stagprec--;
	    res = lnp1(res);  // res = ln[1 + 2/(sqrt(1+x^2)-1) ]
	}

    return res;
} // Aux_1_atan


l_interval Q_atan_UPSIGN(const l_interval& x, const l_interval& y)
{
// x: abs(Re(z));  x: a real interval x = [x1,x2], with 0<=x1<x2.
// y: Inf(Im(z));  y is point interval.
// Q_atan_UPSIGN returns an inclusion of ln[1 + 4y/(x^2+(1-y)^2)]
// Tested in detail; Blomquist, 14.02.2007;
    const int n = 511; 
    l_interval res(0.0),t,t1,t2;
    l_real lr;
    int ex_x,ex_y,ex,s;
    if (y==1.0) 
	if(Sup(x)<1) 
	{   // Now:   y = 1 and minreal <= x1 <= x2 < +1  :
	    lr = Inf(x);   
	    // Calculation of t = ln(x) with scaling:
            // ex = expo(x);
	    // ln(x) = ln(x*2^-ex) + ex * ln(2),
	    // if Inf(x) is too small, i.e. expo(lr[1])<-600
	    ex = expo_gr(lr); // ex < 0;
	    if (ex<-600) 
	    {
		t1 = x;
		Times2pown(t1,-ex);  // t1 = 2^(-ex) * x;
		t = ln(t1) + ex*Ln2_l_interval(); // t = ln(x)
	    } else t = ln(x); // t: ln(x)
	    times2pown(t,1);  // t: 2*ln(x);
	    res = ln(4+sqr(x)) - t; // Blomquist, 14.02.2007.
	}
	else  res = lnp1( sqr(2/x) );
    else
    {  // Now we have:   y>=0 and y<>1
	lr = Sup(x);
	ex_x = expo_gr(lr);
	lr = Sup(y);
	ex   = expo_gr(lr);
	if (ex<-100000) return res;  // y = 0  --->  res = 0.0;
        // Now we have:   ( y>0  and  y<>1 )
	ex_y = ex;
	if (ex_x>ex) ex = ex_x;
        // ex = Maximum{expo(x),expo(y)} 
	if (ex>n) 
	{   // Now scaling to prevent overflow by calculating
            // the denominator: x^2 + (1-y)^2 :
	    s = n-ex-1;
            // s is chosen in such a way that x^2 + (1-y)^2 will be maximal
            // after scaling, but an overflow must not be realized!
	    t = x;  times2pown(t,s);  // fast scaling with 2^(s)
	    t1 = y; times2pown(t1,s); // fast scaling with 2^(s)
	    t2 = sqr(t) + sqr(comp(0.5,s+1)-t1); // t2: denominator
	    stagprec++;
	    t = y / t2; // scaled quotient 
	    stagprec--;
	    times2pown(t,2*s+2); // back-scaling with 2^(s+2); 
                                 //  '+2': factor 4 !!
	    if (Inf(t)<0) SetInf(t,0.0);
	    res = lnp1(t);
	}
	else
	{   // Calculating x^2 + (1-y)^2 now without overflow!!!
            // Now we have:   ( y>0  and  y<>1 )
	    t = 1 - y;
	    lr = Inf(t);
	    ex = expo_gr(lr);
	    if (ex_x>ex) ex = ex_x; // It holds:  -1073 <= ex;
	    ex--;
            // Now x^2 + (1-y)^2 is scaled, so that:
            // [x*2^(-ex)]^2 + [(1-y)*2^(-ex)]^2 '=' 1
	    Times2pown(t,-ex); // t = (1-y)*2^(-ex)
	    t1 = x;
	    Times2pown(t1,-ex); // t = x * 2^(-ex)
	    res = sqr(t1) + sqr(t); // res: scaled denominator:
            // res = (x*2^-ex)^2 + [(1-y)*2^-ex]^2;
	    s = ex_y+2-2*ex;
            // For  x --> 0 and y --> 1  it holds:
            // 2^s is the order of 2^(-2ex) * [x^2+(1+y)^2].
	    if (s < 1020)
	    {   // No overflow by calculating y * 2^(-2*ex+2)
                // is not possible:
		t = y;
		times2pown(t,-2*ex+2); // t = y * 2^(-2*ex+2)
		stagprec++;
		t = t / res;
		stagprec--;
		res = lnp1(t);
	    }
	    else // s >= 1020;
	    {    // For  x --> 0 and y --> 1  it holds:
                 // 2^s is the order of 2^(-2ex) * [x^2+(1+y)^2]. 
                 // Thus, for s >= 1020 an overflow by calculating
                 // 2^(-2ex) * [x^2+(1+y)^2] is possible and an
                 // appropriate scaling is necessary:
		ex_y = s-1020;
		if (ex_y%2 != 0) ex_y++; // ex_y = d, now divisible by 2.
		s = ex_y/2;  // s = d/2;
		t1 = x;
		times2pown(t1,-ex-s);
		t2 = 1+y;
		times2pown(t2,-ex-s);
		t = sqr(t1) + sqr(t2);
		stagprec++;
		t = t / res;
		stagprec--;
		res = ex_y * Ln2_l_interval() + ln(t);
	    }
	} 
    }

    return res;
} // Q_atan_UPSIGN


l_cinterval atan( const l_cinterval& z ) throw() //----- 040912 --
{
    l_interval
	rez = Re(z),
	imz = Im(z);

    l_real
	irez = Inf(rez),
	srez = Sup(rez),
	iimz = Inf(imz),
	simz = Sup(imz);

    const int n = 511; // For possible scaling

    l_interval
	hxl(irez), hxu(srez), hyl(iimz), hyu(simz); // all theses intervals are point intervals!

    l_real
	resxl, resxu, resyl, resyu;
  //
  //  1st: check for singularities
  //
    if( ( irez <= 0.0 && srez >= 0.0 ) && ( iimz <= -1.0 || simz >= 1.0 ) )
	cxscthrow(STD_FKT_OUT_OF_DEF("l_cinterval atan( const l_cinterval& z ); z contains singularities."));
  //
  //  2nd: real part
  //  Re( atan( z ) ) = Arg( w ) / 2, where w = 1 - x^2 - y^2 + i * 2x )
  //
  //  evaluate atan on vertical boundaries
  //
    l_interval
//      y_sqr = sqr( imz ),
//      rew_l = (1 - y_sqr) - sqr( hxl ),  // Blomquist; before: rew_l = 1 - sqr(hxl) - y_sqr, 
//      rew_u = (1 - y_sqr) - sqr( hxu );  // Blomquist; before: rew_u = 1 - sqr(hxu) - y_sqr; 
	rew_l, rew_u;

//  ------------------------------ Blomquist ---------------------------------
//  ----------       Improvements for Im(z) = [1,1]  or  Im(z) = [-1,-1]  
    bool sqrImz_1 = (iimz==simz) && (iimz==1.0 || iimz==-1.0); 
                             // Test for Im(z) = [1,1] or [-1,-1]

    if (sqrImz_1) { 
	rew_l = -abs(hxl);  hxl = l_interval(sign(irez)); 
	rew_u = -abs(hxu);  hxu = l_interval(sign(srez));
    }
    else {
	int ex,s;
	l_interval imz_, scf;
	int ex1 = expo_gr(iimz);   
        int ex2 = expo_gr(simz); 
	if (ex2>ex1) ex1 = ex2;

	ex = expo_gr(irez); 
	if(ex1>ex) ex = ex1; // ex = Maximum
	if (ex>n) { // Scaling necessary
	    s = n - ex - 1;
	    scf = l_interval(comp(0.5,s+1)); // scf: scaling factor 2^s
	    times2pown(scf,s); //  scf = 2^(2*s)
	    times2pown(hxl,s); // hxl = hxl * 2^s
	    imz_ = imz;
	    times2pown(imz_,s); // imz_ = imz_ * 2^s
	    rew_l = (scf - sqr(imz_)) - sqr(hxl); // here now without overflow!! 
	    times2pown(hxl,s); // hxl = hxl * 2^s
	} else // rew_l = (1 - sqr( imz )) - sqr( hxl );
	    re_atan(imz,hxl,rew_l); 

	ex = expo_gr(srez);  // 
	if(ex1>ex) ex = ex1; // Maximum
	if (ex>n) { // Scaling necessary
	    s = n - ex - 1;
	    scf = l_interval(comp(0.5,s+1)); // scf: scaling factor 2^s
	    times2pown(scf,s); //  scf = 2^(2*s)
	    times2pown(hxu,s); // hxu = hxu * 2^s
	    imz_ = imz;
	    times2pown(imz_,s); // imz_ = imz_ * 2^s
	    rew_u = (scf - sqr(imz_)) - sqr(hxu); // here now without overflow!!
	    times2pown(hxu,s); // hxu = hxu * 2^s
	} else // rew_u = (1 - sqr( imz )) - sqr( hxu );
	    re_atan(imz,hxu,rew_u);
    }
//  ------------------------------ Blomquist; 22.02.05; --------------------  

  //
  //  left boundary
  //
    l_real rew_inf = Inf( rew_l );
    l_real rew_sup = Sup( rew_l );
    re_vert( irez, hxl, rew_inf, rew_sup, resxl, resxu );

  //
  //  right boundary
  //
    rew_inf = Inf( rew_u );
    rew_sup = Sup( rew_u );
    l_real res_l, res_u;
    re_vert( srez, hxu, rew_inf, rew_sup, res_l, res_u );

    if (res_l<resxl) resxl = res_l;
    if (res_u>resxu) resxu = res_u;

  //
  //  look for extremal values on horizontal boundaries
  //  since atan( x+iy ) = atan( x-iy ),
  //  intersections can be considered in the upper half plane
  //
    l_real abs_y_min = Inf( abs( imz ) );

    if( abs_y_min > 1.0 )
    {
	l_interval
	    abs_hyl = l_interval( abs_y_min ),
//      abs_hxl = sqrt( sqr( abs_hyl ) - 1.0 );
	    abs_hxl = sqrtx2m1(abs_hyl);  // Blomquist;

	if( Sup( abs_hxl ) > irez && Inf( abs_hxl ) < srez )
	//  extremal curve intersects lower boundary of x+i|y| in quadrant I
	//  intersection in Q I or Q IV: update minimum
	//	resxl = inf( atan( abs_y_min / abs_hxl ) ) / 2.0;
	    resxl = Inf( (Pi_l_interval() - atan( 1.0 / abs_hxl ))/2.0 );
	else if( -Inf( abs_hxl ) > irez && -Sup( abs_hxl ) < srez )
	//  extremal curve intersects lower boundary of x+i|y| in quadrant II
	//  intersection in Q II or Q III: update maximum
	    resxu = Sup( (atan( 1.0 / abs_hxl ) - Pi_l_interval())/2.0 );
    }

  //  3rd: imaginary part
  //  Im( atan( z ) ) = +/- Ln( 1 +/- 4y/( x^2 + (1 -/+ y)^2 ) ) / 4
  //
  //  evaluate atan on horizontal boundaries
    l_interval
	abs_rez = abs(rez), 
	im_atan_l, im_atan_u;

    if( iimz < 0.0 )
//    im_atan_l = -ln( 1 - 4 * hyl / ( x_sqr + sqr( 1 + hyl ) ) ) / 4.0;
//    im_atan_l = -lnp1(-4 * hyl / ( x_sqr + sqr( 1 + hyl ) )) / 4.0;  // Blomquist
	im_atan_l = -Q_atan_UPSIGN(abs_rez,-hyl);  // Blomquist 
    else
//    im_atan_l = ln( 1 + 4 * hyl / ( x_sqr + sqr( 1 - hyl ) ) ) / 4.0;
	im_atan_l = Q_atan_UPSIGN(abs_rez,hyl);  // Blomquist
    times2pown(im_atan_l,-2); // Division by 4

    if( simz < 0.0 )
//    im_atan_u = -ln( 1 - 4 * hyu / ( x_sqr + sqr( 1 + hyu ) ) ) / 4.0;
//    im_atan_u = -lnp1(-4 * hyu / ( x_sqr + sqr( 1 + hyu ) ) ) / 4.0; // Blomquist
	im_atan_u = -Q_atan_UPSIGN(abs_rez,-hyu);  // Blomquist
    else
//    im_atan_u = ln( 1 + 4 * hyu / ( x_sqr + sqr( 1 - hyu ) ) ) / 4.0;
	im_atan_u = Q_atan_UPSIGN(abs_rez,hyu);  // Blomquist
    times2pown(im_atan_u,-2); // Division by 4

    resyl = min( Inf( im_atan_l ), Inf( im_atan_u ) );
    resyu = max( Sup( im_atan_l ), Sup( im_atan_u ) );
  //
  //  look for extremal values on vertical boundaries,
  //  if vertical boundaries intersect extremal curves
  //
    l_real abs_x_min = Inf( abs( rez ) );
    l_interval
	x_extr = l_interval( abs_x_min ),
//    y_extr = sqrt( 1.0 + sqr( x_extr ) );
	y_extr = sqrt1px2(x_extr);                     // Blomquist;

    if( Inf( y_extr ) < simz && Sup( y_extr ) > iimz )
    //  extremal curve intersects left boundary of |x|+iy in quadrant I
    //  update maximum
    //  resyu = Sup( ln( 1 + 4 * y_extr / ( sqr( x_extr ) + sqr( 1 - y_extr ) ) ) ) / 4.0;
    //	resyu = Sup( Aux_1_atan(abs_x_min)/4.0 );    // Blomquist
    {
	rez = Aux_1_atan(abs_x_min);
	times2pown(rez,-2);
	resyu = Sup(rez);
    }

    if( -Sup( y_extr ) < simz && -Inf( y_extr ) > iimz )
    //  extremal curve intersects left boundary of |x|+iy in quadrant IV
    //  update minimum
    //  resyl = -Sup( ln( 1 + 4 * y_extr / ( sqr( x_extr ) + sqr( 1 - y_extr ) ) ) ) / 4.0;
    //	resyl = -Sup( Aux_1_atan(abs_x_min)/4.0 );  // Blomquist
    {
	rez = Aux_1_atan(abs_x_min);
	times2pown(rez,-2);
	resyl = -Sup(rez);
    }

    return l_cinterval( l_interval( resxl, resxu ), l_interval( resyl, resyu ) );

}
//
//-- end atan -----------------------------------------------------------------


//-- acot ----------------------------------------------------------- 040912 --
//
//  Analytic inverse cotangent function
//  acot( z ) = atan( 1/z )
//  The code of acot( z ) is almost identical to the code of atan( z )
//
l_cinterval acot( const l_cinterval& z ) throw()
{
 l_interval
    rez = Re(z),
    imz = Im(z);

  l_real
    irez = Inf(rez),
    srez = Sup(rez),
    iimz = Inf(imz),
    simz = Sup(imz);

  const int n = 511; // For possible scaling

 l_interval
    hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

  l_real
    resxl, resxu, resyl, resyu;
  //
  //  1st: check for singularities
  //
  if( ( irez <= 0.0 && srez >= 0.0 ) && ( iimz <= 1.0 && simz >= -1.0 ) )
    cxscthrow(STD_FKT_OUT_OF_DEF("l_cinterval acot( const l_cinterval& z ); z contains singularities."));
  //
  //  2nd: real part
  //  Re( atan(  z  ) )   = Arg( w ) / 2, where w = 1 - x^2 - y^2 + i * 2x )
  //  Re( atan( 1 / z ) ) = Arg( w ) / 2, where w = x^2 + y^2 - 1 + i * 2x )
  //
  //  evaluate acot on vertical boundaries
  //
 l_interval
//    y_sqr = sqr( imz ),
//    rew_l = (y_sqr - 1) + sqr(hxl),
//    rew_u = (y_sqr - 1) + sqr(hxu);
//    rew_l = (sqr( hxl )-1) + y_sqr, 
//    rew_u = (sqr( hxu )-1) + y_sqr; 
      rew_l, rew_u;
//  ------------------------------ Blomquist ------------------------------
//  ----------       Improvements for Im(z) = [1,1]  or  Im(z) = [-1,-1] 
  bool sqrImz_1 = (iimz==simz) && (iimz==1.0 || iimz==-1.0); 
                     // Test for Im(z) =  [1,1] or [-1,-1]

  if (sqrImz_1) { 
      rew_l = abs(hxl);  hxl = l_interval(sign(irez)); 
      rew_u = abs(hxu);  hxu = l_interval(sign(srez));
 }
  else {
      int ex,s;
//      l_real scf; // Scaling factor
      l_interval imz_, scf;
      int ex1 = expo_gr(iimz);  int ex2 = expo_gr(simz); 
      if (ex2>ex1) ex1 = ex2;

      ex = expo_gr(irez);
      if(ex1>ex) ex = ex1; // Maximum
      if (ex>n) { // Scaling necessary
	  s = n - ex - 1;
	  scf = l_interval(comp(0.5,s+1)); // scf: scaling factor 2^s
	  times2pown(scf,s); // scf = 2^(2*s);
	  times2pown(hxl,s); // hxl = hxl * 2^s
	  imz_ = imz;
	  times2pown(imz_,s); // imz_ = imz_ * 2^s
	  rew_l = (sqr(imz_) - scf) + sqr(hxl); // here now without overflow!!
	  times2pown(hxl,s); // hxl = hxl * 2^s
      } else // rew_l = (sqr( imz ) - 1.0) + sqr( hxl );
      {
	  re_atan(imz,hxl,rew_l);
	  rew_l = -rew_l;
      }

      ex = expo_gr(srez);
      if(ex1>ex) ex = ex1; // Maximum
      if (ex>n) { // Scaling necessary
	  s = n - ex - 1;
	  scf = l_interval(comp(0.5,s+1)); // scf: scaling factor 2^s
	  times2pown(scf,s); // scf = 2^(2*s);
	  times2pown(hxu,s); // hxu = hxu * 2^s
	  imz_ = imz;
	  times2pown(imz_,s); // imz_ = imz_ * 2^s
	  rew_u = (sqr(imz_) - scf) + sqr(hxu); // here now without overflow!!
	  times2pown(hxu,s); // hxu = hxu * 2^s
      } else // rew_u = (sqr( imz )-1.0) + sqr( hxu );
      {
	  re_atan(imz,hxu,rew_u);
	  rew_u = -rew_u;
      }
  }
//  ------------------------------ Blomquist; 22.02.05; ------------------  

  //
  //  left boundary
  //
  l_real rew_inf = Inf( rew_l );
  l_real rew_sup = Sup( rew_l );
  re_vert( irez, hxl, rew_inf, rew_sup, resxl, resxu );
  //
  //  right boundary
  //
  rew_inf = Inf( rew_u );
  rew_sup = Sup( rew_u );
  l_real res_l, res_u;
  re_vert( srez, hxu, rew_inf, rew_sup, res_l, res_u );

  if (res_l<resxl) resxl = res_l;
  if (res_u>resxu) resxu = res_u;

  //
  //  look for extremal values on horizontal boundaries
  //  since acot( x+iy ) = acot( x-iy ),
  //  intersections can be considered in the upper half plane
  //
  l_real abs_y_min = Inf( abs( imz ) );

  if( abs_y_min > 1.0 )
  {
      l_interval
	  abs_hyl = l_interval( abs_y_min ),
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
 l_interval
//  x_sqr = sqr( rez ), // overflow is avoided by calling Q_atan_UPSIGN(...)
    im_atan_l, im_atan_u,
    abs_rez = abs(rez);  // Blomquist;
  if( iimz < 0.0 )
//    im_atan_l = -ln( 1 - 4 * hyl / ( x_sqr + sqr( 1 + hyl ) ) ) / 4.0;
      im_atan_l = -Q_atan_UPSIGN(abs_rez,-hyl);  // Blomquist
  else
//    im_atan_l = ln( 1 + 4 * hyl / ( x_sqr + sqr( 1 - hyl ) ) ) / 4.0;
      im_atan_l = Q_atan_UPSIGN(abs_rez,hyl);  // Blomquist
  times2pown(im_atan_l,-2); // Division by 4;

  if( simz < 0.0 )
//    im_atan_u = -ln( 1 - 4 * hyu / ( x_sqr + sqr( 1 + hyu ) ) ) / 4.0;
      im_atan_u = -Q_atan_UPSIGN(abs_rez,-hyu);  // Blomquist
  else
//    im_atan_u = ln( 1 + 4 * hyu / ( x_sqr + sqr( 1 - hyu ) ) ) / 4.0;
      im_atan_u = Q_atan_UPSIGN(abs_rez,hyu);  // Blomquist
  times2pown(im_atan_u,-2); // Division by 4;

  resyl = min( Inf( im_atan_l ), Inf( im_atan_u ) );
  resyu = max( Sup( im_atan_l ), Sup( im_atan_u ) );
  //
  //  look for extremal values on vertical boundaries,
  //  if vertical boundaries intersect extremal curves
  //
  l_real abs_x_min = Inf( abs( rez ) );
  l_interval
      x_extr = l_interval( abs_x_min ),
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

  return l_cinterval(l_interval( resxl, resxu ), l_interval( -resyu, -resyl ) );

}
//
//-- end acot -----------------------------------------------------------------


//-- atanh ---------------------------------------------------------- 040912 --
//
l_cinterval atanh( const l_cinterval& z ) throw()
//
//  atanh( z ) = - i * atan( i * z )
//
{
  l_cinterval res = atan( l_cinterval( -Im(z), Re(z) ) );
  return l_cinterval( Im(res), -Re(res) );
}
//
//-- end atanh ----------------------------------------------------------------

//-- acoth ---------------------------------------------------------- 040912 --
//
l_cinterval acoth( const l_cinterval& z ) throw()
//
//  acoth( z ) = i * acot( i * z )
//
{
  l_cinterval res = acot( l_cinterval( -Im(z), Re(z) ) );
  return l_cinterval( -Im(res), Re(res) );
}
//
//-- end acoth ----------------------------------------------------------------

} // namespace cxsc

/*

  End of File: l_cimath.cpp

*/
