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

/* CVS $Id: lx_cinterval.cpp,v 1.12 2014/01/30 17:23:47 cxsc Exp $ */

#include "lx_cinterval.hpp"

namespace cxsc {
	
// -----------------------------------------------------------------------
// --------- Functions and operators related to type lx_cinterval --------
// -----------------------------------------------------------------------

l_cinterval & l_cinterval::operator = (const lx_cinterval &a) throw()
// This operator is declared in l_cinterval.hpp !!
{
    l_interval lre,lim;
    lre = Re(a);
    lim = Im(a);
    l_cinterval lc(lre,lim);

    return (*this) = lc;
}

cinterval & cinterval::operator = (const lx_cinterval &a) throw()
// This operator is declared in cinterval.hpp !!
{
    l_cinterval lc;
    cinterval z;

    lc = a;
    z = lc;

    return (*this) = z;
}

// ------------------------ Input --------------------------------------
/*
std::string & operator >> (std::string &s, lx_cinterval &a) throw()
// Funktionierende Version mit Zweier-Exponenten
{
    lx_interval Lar, Lai;

    s = skipwhitespacessinglechar (s, '(');
    s = s >> SaveOpt >> Lar;
    s = skipwhitespacessinglechar (s, ',');
    s = s >> Lai >> RestoreOpt;
    s = skipwhitespacessinglechar (s, ')');
    s = skipwhitespaces (s); // erasing whitespaces
    a = lx_cinterval(Lar,Lai);

    return s;
}
*/
std::string & operator >> (std::string &s, lx_cinterval &a) throw()
// Funktionierende Version mit Zehner-Exponenten, z.B.
// s = ( {+1, [4,4.01]} , {-2, [7,7]} ) 
{
    lx_interval Lar, Lai;
    string str(s);
    int i;

    str = skipwhitespacessinglechar (str, '(');
    i = str.find("}");
    str.erase(i+1);
    str >> SaveOpt >> Lar;

    i = s.find("}");
    s.erase(0,i+1);
    s = skipwhitespacessinglechar (s, ',');
    s >> Lai >> RestoreOpt;
    s = "";

    a = lx_cinterval(Lar,Lai);

    return s;
}

void operator >> (const std::string &s, lx_cinterval &a) throw()
{
// Writes strings s to variable a of type lx_cinterval;
    std::string r(s);
    r >> a;
}

void operator >> (const char *s, lx_cinterval &a) throw()
{
    std::string r(s);
    r >> a;
}

std::istream & operator >> (std::istream &s, lx_cinterval &a) throw()
{
    lx_interval Lar, Lai;
    char c;
	 
	 std::cerr << "Real part: {Exponent to base 10, [a,b]} = ?" 
			     << std::endl;
	 s >> Lar;
	 std::cerr << "Img. part: {Exponent to base 10, [a,b]} = ?" 
			     << std::endl;	
	 s >> Lai >> RestoreOpt;
    a = lx_cinterval(Lar,Lai);

    if (!waseolnflag)
    {
		skipeolnflag = false; inpdotflag = true;
		c = skipwhitespaces (s);
		if (inpdotflag && c != ')')
	  		s.putback(c);
    }
    return s;
}


// -----------------------------------------------------------------------------
// ---------- Elementary functions related to type lx_cinterval ----------------
// -----------------------------------------------------------------------------

// --------------------------------------------------------------------------
// ------------------------- Help functions ---------------------------------
// --------------------------------------------------------------------------

	int max(int a, int b)
	{
		int res(a);
		if (b>a) res = b;
		return res;
	}


// --------------------------------------------------------------------------
// --------------------------- sqr function ---------------------------------
// --------------------------------------------------------------------------

lx_cinterval sqr(const lx_cinterval &z) throw()
{
	lx_interval rez(Re(z)), reza(abs(rez)),
					imz(Im(z)), imza(abs(imz));
	lx_real
		irez = Inf(reza),
		srez = Sup(reza),
		iimz = Inf(imza),
		simz = Sup(imza);
	lx_interval
		hxl(irez), hxu(srez), hyl(iimz), hyu(simz);
	lx_real resxl, resxu;
	resxl = Inf( (hxl-hyu)*(hxl+hyu) );
	resxu = Sup( (hxu-hyl)*(hxu+hyl) );

	hxl = rez * imz;
	times2pown(hxl,1);

	return lx_cinterval( lx_interval(resxl,resxu), hxl );
} // sqr(...)

// -------------------------------------------------------------------------
// ---------------------- Begin square root --------------------------------
// -------------------------------------------------------------------------

	lx_interval Sqrt_zpx_m2( const lx_interval &x, const lx_interval &y )
    // z = x + i*y;
    // Calculating res = sqrt( 2*(|z| + |x|) );
	{
		lx_interval res;
		res = sqrtx2y2(x,y) + abs(x);  // New: 22.05.2008;
		times2pown(res,1);
		return sqrt(res);
	}

	lx_interval Sqrt_zpx_d2( const lx_interval &x, const lx_interval &y )
    // z = x + i*y;
    // Calculating res = sqrt( (|z| + |x|)/2 );
	{
		lx_interval res;
		res = sqrtx2y2(x,y) + abs(x);  // New: 22.05.2008;
		times2pown(res,-1);
		return sqrt(res);
	}

lx_interval Re_Sqrt_Point(const lx_interval &rez, const lx_interval &imz)
    // Real part of analytic square root of a POINT INTERVAL only.
    // Do not use this as a general function
    // - it is only a subroutine for sqrt(...)
    // Blomquist, 03.05.2007;
{
	lx_interval res;
	lx_real irez = Inf( rez ), iimz = Inf( imz );

	if (eq_zero(iimz))
		res = (ge_zero(irez))? sqrt(rez) : lx_interval(0);
	else // iimz <> 0
		res = (ge_zero(irez))? Sqrt_zpx_d2(rez,imz) : 
				 abs(iimz)/Sqrt_zpx_m2(rez,imz);
	return res;
}

lx_interval Im_Sqrt_Point(const lx_interval &rez, const lx_interval &imz)
    // Imaginary part of analytic square root of a POINT INTERVAL only.
    // Do not use this as a general function
    // - it is only a subroutine for sqrt(...)
    // Blomquist, 03.05.2007;
{
	lx_interval res;
	lx_real irez = Inf( rez ), iimz = Inf( imz );

	if (eq_zero(iimz))
		res = (ge_zero(irez))? lx_interval(0) : sqrt(-rez);
	else
		if (ge_zero(irez)) res = iimz/Sqrt_zpx_m2(rez,imz);
		else 
		{
			res = Sqrt_zpx_d2(rez,imz);
			if (se_zero(iimz)) res = -res;
		}
	return res;
}

lx_cinterval sqrt(const lx_cinterval& z) throw()
{
	lx_cinterval y;
	lx_real
		irez = Inf( Re(z) ),
		srez = Sup( Re(z) ),
		iimz = Inf( Im(z) ),
		simz = Sup( Im(z) );
	lx_interval 
		hxl( irez ), hxu( srez ), hyl( iimz ), hyu( simz );
	lx_real
		resxl, resxu, resyl, resyu;

	if ( sm_zero(irez) && sm_zero(iimz) && ge_zero(simz) )
		cxscthrow(STD_FKT_OUT_OF_DEF(
		"lx_cinterval sqrt(const lx_cinterval& z); z not in principal branch."));

	if (ge_zero(iimz))
	{
		resxl = Inf( Re_Sqrt_Point( hxl,hyl ) );
		resxu = Sup( Re_Sqrt_Point( hxu,hyu ) );

		resyl = Inf( Im_Sqrt_Point( hxu,hyl ) );
		resyu = Sup( Im_Sqrt_Point( hxl,hyu ) );
	}
	else
		if (se_zero(simz))
		{
			resxl = Inf( Re_Sqrt_Point( hxl,hyu ) );
			resxu = Sup( Re_Sqrt_Point( hxu,hyl ) );

			resyl = Inf( Im_Sqrt_Point( hxl,hyl ) );
			resyu = Sup( Im_Sqrt_Point( hxu,hyu ) );
		}
		else
		{
			resxl = Inf( sqrt(hxl) );
			resxu = ( -iimz>simz ? Sup( Re_Sqrt_Point( hxu,hyl ) )
						: Sup( Re_Sqrt_Point( hxu,hyu ) ) );

			resyl = Inf( Im_Sqrt_Point( hxl,hyl ) );
			resyu = Sup( Im_Sqrt_Point( hxl,hyu ) );
		}
		y = lx_cinterval( lx_interval(resxl,resxu), lx_interval(resyl,resyu) );

		return y;
} // sqrt(...)
	 
//  sqrt_all(z) computes a list of 2 intervals containing all square roots of z
//
std::list<lx_cinterval> sqrt_all( const lx_cinterval& z ) throw()
{
	lx_real
		irez = Inf( Re(z) ),
		srez = Sup( Re(z) ),
		iimz = Inf( Im(z) ),
		simz = Sup( Im(z) );
	lx_interval hxl( irez ), hxu( srez ), hyl( iimz ), hyu( simz );
	lx_real resxl,  resxu, resyl, resyu;
	lx_cinterval w;

	if( irez < 0.0 && iimz <= 0.0 && simz >= 0.0 )
   //  z contains negative real values
	{
		if( iimz == 0.0 )
	    //  z in upper half plane
	    //  principal values can be used
		{
	    //  min( Re ( sqrt( z ) ) ) in lower left  corner
	    //  max( Re ( sqrt( z ) ) ) in upper right corner
			resxl = Inf( Re_Sqrt_Point( hxl, hyl ) );
			resxu = Sup( Re_Sqrt_Point( hxu, hyu ) );
	   //  min( Im ( sqrt( z ) ) ) in lower right corner
	   //  max( Im ( sqrt( z ) ) ) in upper left  corner
			resyl = Inf( Im_Sqrt_Point( hxu, hyl ) );
			resyu = Sup( Im_Sqrt_Point( hxl, hyu ) );
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
				resxl = 0;
				resxu = Sup( Re_Sqrt_Point( hxu, hyl ) );
	         //  min( Im ( sqrt( z ) ) ) in lower left  corner
	         //  max( Im ( sqrt( z ) ) ) in upper right corner
				resyl = Inf( Im_Sqrt_Point( hxl, hyl ) );
				resyu = ( srez > 0.0 ? lx_real(0.0,l_real(0)) : -Inf(sqrt(-hxu) ) );
			}
			else
	      //  0 is interior point of Im( z )
			{
				if( srez <= 0.0 )
				{
		       //  0 is no interior point of Re (z )
		       //  in quadrant III,    Re( z ) = Im_Sqrt_Point(|x|,y),
		       //	                   Im( z ) = Re_Sqrt_Point(|x|,y)
		  //  min( Re ( sqrt( z ) ) ) in lower right corner = quadrant II/III
		  //  max( Re ( sqrt( z ) ) ) in upper right corner = quadrant II
					resxl = Inf( Im_Sqrt_Point(-hxu, hyl ) );
					resxu = Sup( Re_Sqrt_Point( hxu, hyu ) );
		  //  min( Im ( sqrt( z ) ) ) on real line
		  //  max( Im ( sqrt( z ) ) ) in lower or upper left corner
					resyl = Inf( sqrt( -hxu ) );
					resyu = ( - iimz > simz ? Sup( Re_Sqrt_Point( -hxl, hyl ) ) :
								Sup( Im_Sqrt_Point( hxl, hyu ) ) );
				}
				else
		      //  0 is interior point of z
		      //  here, the principal values apply in all corner points
				{
					resxl = 0;
		         //  max( Re ( sqrt( z ) ) ) in lower or upper right corner
					resxu = ( - iimz > simz ? Sup( Re_Sqrt_Point( hxu, hyl ) ) :
							Sup( Re_Sqrt_Point( hxu, hyu ) ) );
		         //  min( Im ( sqrt( z ) ) ) in lower left corner
		         //  max( Im ( sqrt( z ) ) ) in upper left corner
					resyl = Inf( Im_Sqrt_Point( hxl, hyl ) );
					resyu = Sup( Im_Sqrt_Point( hxl, hyu ) );
				}
			}
		}
		w = lx_cinterval( lx_interval(resxl,resxu), lx_interval(resyl,resyu ) );
	}
	else
   //  sqrt( z ) is well-defined
		w = sqrt( z );

	std::list<lx_cinterval> res;
	res.push_back(  w );
	res.push_back( -w );

	return res;
}
//-- end sqrt_all -------------------------------------------------------------
// ---------------------- End square root -----------------------------------

// ***************************************************************************
// ***************************************************************************
// ***                      Single-valued functions                       ****
// ***************************************************************************
// ***************************************************************************


// ***************************************************************************
// *** Power operator  pow  is not listed here, since it relies on the    ****
// *** (multi-valued) logarithm                                           ****
// ***************************************************************************


// ***************************************************************************
// *** The hyperbolic functions exp, sin, cos, sinh, cosh are separable:  ****
// *** Their real and imaginary parts are products of real functions      ****
// ***************************************************************************
// ***   With Re(z)=x, Im(z)=y :                                          ****
// ***                                                                    ****
// ***        exp   :   Re(exp(z)) = exp(x) * cos(y)                      ****
// ***                  Im(exp(z)) = exp(x) * sin(y)                      ****
// ***                                                                    ****
// ***        sin   :   Re(sin(z)) = sin(x) * cosh(y)                     ****
// ***                  Im(sin(x)) = cos(x) * sinh(y)                     ****
// ***                                                                    ****
// ***        cos   :   Re(cos(z)) = cos(x) * cosh(y)                     ****
// ***                  Im(sin(x)) = -sin(x) * sinh(y)                    ****
// ***                                                                    ****
// ***        sinh  :   Re(sinh(z)) = sinh(x) * cos(y)                    ****
// ***                  Im(sinh(z)) = cosh(x) * sin(y)                    ****
// ***                                                                    ****
// ***        cosh  :   Re(cosh(z)) = cosh(x) * cos(y)                    ****
// ***                  Im(cosh(z)) = sinh(x) * sin(y)                    ****
// ***                                                                    ****
// ***************************************************************************

// -------------------------- exp(...) --------------------------------------

	lx_cinterval exp(const lx_cinterval& z) throw()
	{
		int stagsave = stagprec,
  stagmax = 39;
  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval lreal(Re(z)), limg(Im(z));
  lx_cinterval y;
  lx_interval
		  A( exp(lreal) ),
			  B( limg );
			  y = lx_cinterval( A*cos( B ) , A*sin( B ) );

			  stagprec = stagsave;
			  y = adjust(y);

			  return y;
	}

	lx_cinterval exp2(const lx_cinterval& z) throw()
	{
		return exp(z*Ln2_lx_interval());
	}

	lx_cinterval exp10(const lx_cinterval& z) throw()
	{
		return exp(z*Ln10_lx_interval());
	}

// -------------------------- cos(...) --------------------------------------

	lx_cinterval cos(const lx_cinterval& z) throw()
	{
		int stagsave = stagprec,
  stagmax = 39;
  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval lreal(Re(z)),limg(Im(z));
  lx_cinterval y;
  lx_interval
		  A( lreal ),
			  B( limg );
			  y = lx_cinterval( cos( A )*cosh( B ) , -sin( A )*sinh( B ) );

			  stagprec = stagsave;
			  y = adjust(y);

			  return y;
	}

// -------------------------- sin(...) --------------------------------------

	lx_cinterval sin(const lx_cinterval& z) throw()
	{
		int stagsave = stagprec,
  stagmax = 39;
  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval lreal(Re(z)),limg(Im(z));
  lx_cinterval y;
  lx_interval
		  A( lreal ),
			  B( limg );
			  y = lx_cinterval( sin( A )*cosh( B ) , cos( A )*sinh( B ) );

			  stagprec = stagsave;
			  y = adjust(y);

			  return y;
	}

// -------------------------- cosh(...) -------------------------------------

	lx_cinterval cosh(const lx_cinterval& z) throw()
	{
		int stagsave = stagprec,
  stagmax = 39;
  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval lreal(Re(z)),limg(Im(z));
  lx_cinterval y;
  lx_interval
		  A( lreal ),
			  B( limg );
			  y = lx_cinterval( cos( B )*cosh( A ) , sin( B )*sinh( A ) );

			  stagprec = stagsave;
			  y = adjust(y);

			  return y;
	}

// -------------------------- sinh(...) -------------------------------------

	lx_cinterval sinh(const lx_cinterval& z) throw()
	{
		int stagsave = stagprec,
  stagmax = 39;
  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval lreal(Re(z)),limg(Im(z));
  lx_cinterval y;
  lx_interval
		  A( lreal ),
			  B( limg );
			  y = lx_cinterval( cos( B )*sinh( A ) , sin( B )*cosh( A ) );

			  stagprec = stagsave;
			  y = adjust(y);

			  return y;
	}

//-----------------------------------------------------------------------------
	//
//  Part I: Multi-valued functions
	//
//-----------------------------------------------------------------------------
	//
	//

lx_interval Atan(const lx_interval& y,const lx_interval& x)
// Calculating an inclusion of atan(y/x) with x<>[0,0];
// This help function must only be used for POINT intervals y,x !!
// This function avoids internal overflow by calculating y/x.
// Blomquist, 25.05.2008;
{ 
	lx_interval res(0),u;
	l_interval yl(li_part(y));
	int ex_yl(expo_gr(yl)),
		signy(sign(Inf(y))), signx(sign(Inf(x)));
	real ex_y(expo(Inf(y))), ex_x(expo(Inf(x))); 
	bool neg;

	if (ex_yl > -1000000)
	{ // y != 0;
		neg = signy * signx == -1;
		if (ex_y > ex_x+4197)
		{   // We now assume y>0 and x>0:
			res = Pi_lx_interval();
			times2pown(res,-1);
			if (neg) res = -res;
		}
		else
		{
			if (ex_x >= 9007199254738946.0)
			// Now y/x can generate an error!
			{
				if (ex_y<-5217)
				{
						res = lx_interval( lx_real(0,l_real(0)),
												lx_real(-Max_Int_R,l_real(minreal)) );
						if (neg) res = -res;
					}
					else // ex_y >= -5217
					{
						res = x;
						times2pown(res,-2045);
						u = y;
						times2pown(u,-2045);
						res = atan(u/res); // Division now without an error!
					}
				}
				else 
					if (ex_x <= -9007199254738891.0)
					{
						res = x;
						times2pown(res,2101);
						u = y;
						times2pown(u,2101);
						res = atan(u/res);
					}
					else
						res = atan(y/x);
		}
	}

	return res;
}

lx_interval Atan(const lx_interval& y, const lx_real& x)
{
	return Atan(y,lx_interval(x));
}

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
lx_interval Arg(const lx_cinterval& z) throw()
{
	lx_real 
		srez = Sup( Re(z) ),
		irez = Inf( Re(z) ), 
		simz = Sup( Im(z) ),
		iimz = Inf( Im(z) );

	lx_interval
		hxl(irez), hxu(srez), hyl(iimz), hyu(simz), Pid2; // Point intervals

	lx_real resl, resu;

	Pid2 = Pid2_lx_interval();

	if( iimz > 0.0 )
    //  case I: Im(z) > 0
	{
		resl = ( srez > 0.0 ? Inf( Atan( hyl,hxu ) ) : 
	       	 	( srez < 0.0 ? Inf( Atan( hyu,hxu ) + Pi_lx_interval() ) : 
					Inf( Pid2 ) ) );
		resu = ( irez > 0.0 ? Sup( Atan( hyu,hxl ) ) : 
		 		   ( irez < 0.0 ? Sup( Atan( hyl,hxl ) + Pi_lx_interval() ) : 
					Sup( Pid2 ) ) );
		return lx_interval( resl, resu );
	}
	else
	{
		if( simz < 0.0 )
		//  case II: Im(z) < 0
		{
	    	resl = ( irez < 0.0 ? Inf( Atan( hyu,hxl ) - Pi_lx_interval() ) :
				 	   ( irez > 0.0 ? Inf( Atan( hyl,hxl ) ) : -Sup( Pid2 ) ) );
	  		resu = ( srez < 0.0 ? Sup( Atan( hyl,hxu ) - Pi_lx_interval() ) :
			  			( srez > 0.0 ? Sup( Atan( hyu,hxu ) ) : -Inf(Pid2) ) );
	  		return lx_interval( resl, resu );
		}
		else
		// 0 in Im(z)
		{
			if( irez > 0.0 )
	    	//  case III: Re(z) > 0
	   	 //  z contains positive real values
			{
				resl = iimz < 0.0 ? Inf( Atan( hyl,hxl ) ) : lx_real(0.0);
				return lx_interval( resl, Sup( Atan( hyu,hxl ) ) );
			}
			else
	    	//  z contains nonpositive real numbers
			{
				if( irez < 0.0 )
				{
		  		//  case IV: z contains negative real numbers
					cxscthrow (STD_FKT_OUT_OF_DEF("lx_interval Arg(const lx_cinterval& z); z contains negative real numbers"));
					return lx_interval(0.0);
				}
				else
				//  case V: 0 in z, but z doesn't contain negative real numbers
				{
					if( srez > 0.0 )
		    		//  diam( Re(z) > 0.0 )
					{
						resl = iimz < 0.0 ? -Sup(Pid2) : lx_real(0.0);
						resu = simz > 0.0 ?  Sup(Pid2) : lx_real(0.0);
						return lx_interval( resl, resu );
					}
					else
		    		//  Re(z) == 0.0
					{
						if( eq_zero(iimz) && eq_zero(simz) )
						//  Z == 0
							return lx_interval(0.0);
						else
						{
							resl = ( iimz < 0.0 ? - Sup(Pid2) : Inf(Pid2) );
							resu = ( simz > 0.0 ? Sup(Pid2) : -Inf(Pid2) );
							return lx_interval( resl, resu );
						}
					}
				}
			}
		}
	}
}
	//
//-- end Arg ------------------------------------------------------------------

//-- arg: non-analytic argument function --------------------------------------
	//
//  (i)   arg(Z) is defined for all Z in IC.
//  (ii)  arg([0,0]) = 0.
//  (iii) arg(Z) subset [-pi,3*pi/2).
//  (iv)  arg(Z) == Arg(Z) if Arg(Z) is well-defined.
	//

lx_interval arg( const lx_cinterval& z ) throw()
{
	lx_real
		srez = Sup( Re(z) ),
		irez = Inf( Re(z) ),
		simz = Sup( Im(z) ),
		iimz = Inf( Im(z) );

	lx_real resl, resu;
	lx_interval Pid2;
	Pid2 = Pid2_lx_interval();

	if( sm_zero(irez) && se_zero(iimz) && ge_zero(simz) )
   //  z contains negative real values
	{
		if( gr_zero(srez) )
	   //  0 in z and 0 interior point of Re(z)
		{
			resl = ( sm_zero(iimz)? - Sup( Pi_lx_interval() ) : lx_real(0.0) );
	    	resu = ( ( sm_zero(iimz) && eq_zero(simz) ) ? lx_real(0.0) :
				 		Sup( Pi_lx_interval() ) );
		 	return lx_interval( resl, resu );
		}
		else
		{ // srez <= 0.0
			if( iimz == simz )
	    	//  z is real interval containing no positive values
				return Pi_lx_interval();
			else
	    	// sup( Re(z) ) <= 0, diam( Im(z) ) > 0
			{
				if( eq_zero(srez) )
				{
		    		resl = ( gr_zero(simz)? Inf(Pid2) : 
					 			-Sup( Pi_lx_interval() ) );
			 		resu = ( sm_zero(iimz)? 
			     		( gr_zero(simz)? Sup(3*Pid2) : -Inf(Pid2) ) : 
					 		Sup( Pi_lx_interval() ) );
			 		return lx_interval( resl, resu );
				}
				else
				//   sup( Re(z) ) < 0, diam( Im(z) ) > 0
				{
					lx_interval hyl(iimz), hyu(simz);
					resl = ( simz > 0.0 ? Inf( Atan( hyu,srez ) + Pi_lx_interval() ) : 
												 -Sup( Pi_lx_interval() ) );
					resu = ( iimz < 0.0 ? ( simz > 0.0 ? Sup( Atan( hyl,srez ) +
						 Pi_lx_interval() ) : Sup( Atan( hyl,srez ) -
						 Pi_lx_interval() ) ) : Sup( Pi_lx_interval() ) );
					return lx_interval( resl, resu );
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

//-------------------- sqrt(z,n): analytic n-th root --------------------------
	//
	lx_interval Re_Sqrt_point(const lx_interval& rez, const lx_interval& imz,
									  int n ) 
			// before: unsigned int n  ---------- 040624 --
			//
//  Real part of analytic n-th root.
//  Do not use this as a general function
//  - it's only a subroutine for sqrt(z,n) and sqrt_all(z,n).
//  The calculation is validated but largely overestimating
//  if (rez+I*imz) is not a complex number.
//  Blomquist, 25.05.2008;
	{
		lx_interval a = sqrtx2y2(rez,imz);
		if( eq_zero(Sup(a)) )
    //  a == 0
			return lx_interval(0);
		else
			return sqrt(a, n ) * 
					cos( Arg( lx_cinterval( rez, imz ) ) / n );
	}

	lx_interval Im_Sqrt_point(const lx_interval& rez, const lx_interval& imz,
									  int n ) // before: unsigned int n  --- 040624 --
			//
//  Imaginary part of analytic n-th root.
//  Do not use this as a general function
//  - it's only a subroutine for sqrt(z,n) and sqrt_all(z,n).
//  The calculation is validated but largely overestimating
//  if (rez+I*imz) is not a complex number.
//  Blomquist, 25.05.2008;
	{
		lx_interval a = sqrtx2y2(rez,imz);
		if( eq_zero(Sup(a)) )
    //  a == 0
			return lx_interval(0);
		else
			return sqrt(a, n ) * 
					sin( Arg( lx_cinterval( rez, imz ) ) / n );
	}

lx_cinterval sqrt(const lx_cinterval& z, int n ) throw() 
			//
//  Analytic n-th root function
//  sqrt(z,n) is undefined if z contains negative real numbers.
			//
{
	if( n == 0 ) return lx_cinterval(lx_interval(1));
	if( n == 1 ) return z;
	if( n == 2 ) return sqrt(z);
	else
	{
		lx_real
			irez = Inf( Re(z) ),
			srez = Sup( Re(z) ),
			iimz = Inf( Im(z) ),
			simz = Sup( Im(z) );
		lx_interval hxl( irez ), hxu( srez ), hyl( iimz ), hyu( simz );
		lx_real resxl,  resxu, resyl, resyu;

		if( sm_zero(irez) && se_zero(iimz) && ge_zero(simz) )
		{
	   //  z contains negative real values
			cxscthrow(STD_FKT_OUT_OF_DEF("lx_cinterval sqrt(const lx_cinterval& z, int n ); z contains negative real values."));
			return z;
		}
		else 
		{
			if( sm_zero(simz) )
			{
	      //  z in lower half plane
				lx_cinterval hres = sqrt( lx_cinterval( Re(z), -Im(z) ), n );
				return lx_cinterval( Re(hres), -Im(hres) );
			}
			else
			{
				if( iimz > 0.0 ) 
				{
		  		//  z in upper half plane
					lx_interval tangle = tan( ( Pi_lx_interval() * n )
					     / ( 2 * ( n-1 ) ) );
					lx_real tanglel = Inf( tangle ), tangleu = Sup( tangle );
				//
		      //  min( Re( Root( z ) ) )
				//
					if ( ge_zero(irez) ||  Sup( hyl / irez ) <= tanglel )
		         //  lower boundary right of phi = n*Pi/(2n-2)
		         //  min( Re( Root( z ) ) ) in lower left corner
						resxl = Inf( Re_Sqrt_point( hxl, hyl, n ) );
					else
					{
						if( sm_zero(srez) && Inf( hyl / srez ) >= tangleu )
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
					if ( ge_zero(irez) || Sup( hyu / irez ) <= tanglel )
		         //  upper boundary right of phi = n*Pi/(2n-2)
		         //  max( Re( Root( z ) ) ) in upper right corner
						resxu = Sup( Re_Sqrt_point( lx_interval(srez),
										 			lx_interval(simz), n ) );
					else
					{
						if ( sm_zero(srez) && Inf( hyu / srez ) >= tangleu )
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
					if ( ge_zero(srez) || Sup( hyl / srez ) <= tanglel )
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
					if( ge_zero(irez) || Sup( hyl / irez ) <= tanglel )
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
		    		resxl = ( eq_zero(irez)? lx_real(0.0) : Inf( sqrt( hxl, n ) ) );
			 		resxu = ( - iimz > simz ? Sup( Re_Sqrt_point( hxu, hyl, n ) ) :
					 								  Sup( Re_Sqrt_point( hxu, hyu, n ) ) );
		         //  min( Im ( sqrt( z ) ) ) in lower left corner
		         //  max( Im ( sqrt( z ) ) ) in upper left corner
			 		resyl = Inf( Im_Sqrt_point( hxl, hyl, n ) );
			 		resyu = Sup( Im_Sqrt_point( hxl, hyu, n ) );
				}
				return lx_cinterval( lx_interval( resxl, resxu ),
											lx_interval( resyl, resyu ) );
			}
		}
	}
}
//
//-- end sqrt(z,n) -----------------------------------------------------------

//-- sqrt_all ------------------------------------------------------- 040628 --
	//
	std::list<lx_cinterval> sqrt_all( const lx_cinterval& z, int n ) throw()
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
		std::list<lx_cinterval> res;

		if( n == 0 )
		{
			res.push_back( lx_cinterval( lx_interval(0,l_interval(1)),
								lx_interval(0,l_interval(0)) ) );
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
			lx_interval
					arg_z = arg( z ), root_abs_z = sqrt( abs(z), n );

			for(int k = 0; k < n; k++)
			{
				lx_interval arg_k = ( arg_z + 2 * k * Pi_l_interval() ) / n;

				res.push_back( lx_cinterval( root_abs_z * cos( arg_k ),
									root_abs_z * sin( arg_k ) ) );
			}
			return res;
		}
	}
	//
//-- end sqrt_all -------------------------------------------------------------


//-----------------------------------------------------------------------------
	//
//  Section 2: Logarithms
	//
//-----------------------------------------------------------------------------

// ------------------- Ln: analytic natural logarithm -------------------------
//
//  Ln(z) is undefined if z contains zero; z must not touch the negative real
//  axis from below;
//
lx_cinterval Ln(const lx_cinterval& z) throw()
{  // Blomquist, 24.09.2007;
	int stagsave = stagprec,
  	stagmax = 30; 
  	if (stagprec > stagmax) stagprec = stagmax;

  	lx_cinterval y;
  	lx_real
		srez = Sup( Re(z) ),
		simz = Sup( Im(z) ),
		iimz = Inf( Im(z) );
	lx_interval a1( abs(Re(z)) ), a2( abs(Im(z)) );
	if ( eq_zero(Inf(a1)) && eq_zero(Inf(a2)) )
	//  z contains 0
		cxscthrow(STD_FKT_OUT_OF_DEF(
			"lx_cinterval Ln(const lx_cinterval& z); z contains 0"));
	if ( sm_zero(srez) && sm_zero(iimz) && ge_zero(simz) ) 
		cxscthrow(STD_FKT_OUT_OF_DEF(
			"lx_cinterval Ln(const lx_cinterval& z); z not allowed"));
	
	y = lx_cinterval( ln_sqrtx2y2(Re(z),Im(z)), arg(z) );

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
lx_cinterval ln(const lx_cinterval& z) throw()
{
	int stagsave = stagprec,
  		stagmax = 30; 
  	if (stagprec > stagmax) stagprec = stagmax;

  	lx_cinterval y;
  	lx_interval a1( abs(Re(z)) ), a2( abs(Im(z)) );
	if ( eq_zero(Inf(a1)) && eq_zero(Inf(a2)) ) 
	//  z contains 0
		cxscthrow(STD_FKT_OUT_OF_DEF
			("lx_cinterval ln(const lx_cinterval& z); z contains 0"));
	y = lx_cinterval( ln_sqrtx2y2(Re(z),Im(z)), arg(z) );

	stagprec = stagsave;
	y = adjust(y);

	return y; 
}
//
//-- end ln -------------------------------------------------------------------

// ---------------------- log2, log10: analytic functions ---------------------
//
//  log2(z),log10 are undefined if z contains zero; z must not touch the 
//  negative real axis from below;
//
lx_cinterval log2(const lx_cinterval& z) throw()
{  // Blomquist, 30.11.2008;
	return Ln(z) / Ln2_lx_interval();
}

lx_cinterval log10(const lx_cinterval& z) throw()
{  // Blomquist, 30.11.2008;
	return Ln(z) / Ln10_lx_interval();
}

//-----------------------------------------------------------------------------
//
//  Section 3: Power functions
//
//-----------------------------------------------------------------------------

//-- power_fast ---------------------------------------------------------------
//
//  Fast, validated power function for integer powers, based on Moivre's formula.
//  If n is not an integer value, an error message is generated.
//  Medium amount of overestimation.
//
lx_cinterval power_fast(const lx_cinterval& z, const real& n) throw()
{
	if( n == 0 ) return lx_cinterval(lx_interval(1));
	else if( n == 1 ) return z;
		  else if( n == -1 ) return 1 / z;
		       else if( n == 2 ) return sqr(z); 
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
    //  In power_fast(z,n), z is enclosed into a sector S of
    //  an annulus. S^n is again some sector S' of a (different) annulus.
    //  S' is calculated exactly (apart from roundoff), and then enclosed
    //  into a rectangle. There is a certain amount of overestimation
    //  but the calculations are quit cheap.
    //  The implementation is base on Moivre's formula.
	{
		lx_interval abs_z = abs(z);

		if( ((n < 0) && eq_zero(Inf(abs_z))) || !(Is_Integer(n)))
		//  z contains 0
			cxscthrow (STD_FKT_OUT_OF_DEF("lx_cinterval power_fast(const lx_cinterval& z, const real& n); z contains 0 or n is not integer."));

		if( eq_zero(Sup(abs_z)) ) return lx_cinterval( lx_interval(0));
		else
		{	
			lx_interval arg_z = arg(z);
			lx_interval abs_z_n = power(abs_z,n); // Legendre algorithm

			return lx_cinterval( abs_z_n * cos( n * arg_z ),
										abs_z_n * sin( n * arg_z ) );
		}
	}
}
//
//-- end power_fast -----------------------------------------------------------

//------------------------------------ power ----------------------------------
//
//  power function based on the Legendre algotizkm with an integer valued 
//  exponent n of type real.
//

	lx_cinterval power( const lx_cinterval& x, const real& n ) throw()
	{
		if( !(Is_Integer(n)) )
	//  n is not an integer
			cxscthrow(STD_FKT_OUT_OF_DEF("lx_cinterval power(const lx_cinterval& z, const real& n); n is not integer."));
	
		real one(1.0), zhi(2.0), N(n), r;
		double dbl;
		lx_cinterval y,neu,X(x);
	
		if (x == one) y = x;
		else 
			if (N == 0.0) y = one; 
		else 
		{
			if (N == 1) y = x;
			else 
				if (N == 2) y = sqr(x);
			else 
			{
				if (N < 0) 
				{
					X = 1.0/X;
					N = -N;
				}
		    // Initialisierung
				if ( !Is_Integer(N/2) )
					y = X;
				else y = one; 
				neu = sqr(X);   // neu = X*X;
				do {
					dbl = _double(N/zhi);
					dbl = floor(dbl);
					r = (real) dbl;
					if ( !Is_Integer( r/2 ) )
						y *= neu;
					zhi += zhi;
					if (zhi <= N) 
						neu = sqr(neu); // neu = neu * neu;
				} while (zhi <= N);
			}
		}
		return y;
	}
	//
//-- end power ----------------------------------------------------------------

//------------------------------------ pow ------------------------------------
	//
//  Analytic power function for a real interval exponent, based on Ln.
	//

	lx_cinterval pow( const lx_cinterval& z, const lx_interval& p ) throw()
	{
		return exp( p*Ln(z) );
	}

	//
//-- end pow ------------------------------------------------------------------

//------------------------------------ pow ------------------------------------
	//
//  Analytic power function for a complex interval exponent, based on Ln.
	//

	lx_cinterval pow(const lx_cinterval& z, const lx_cinterval& p) throw()
	{
		return exp( p*Ln(z) );
	}

//
//-- end pow ------------------------------------------------------------------

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

//---------------------------------- tan --------------------------------------
//
//  Complex tangent function
//

void horizontal_check( const lx_interval& hy, const lx_interval& cosh_2y,
							  lx_real irez, lx_real srez, const lx_interval& hxl,
								const lx_interval& hxu, lx_real& resxl, lx_real& resxu )
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
	lx_interval hlp, hlp1;

	hlp = Pid2_lx_interval();
	hlp1 = lx_interval(srez) - lx_interval(irez);
	if (Inf(hlp1) > Inf(hlp))
   //  2 intersections
		both = true;
	else
	{
		lx_interval
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
				lx_interval
					sin_2xl = sin( 2 * hxl ),
					sin_2xu = sin( 2 * hxu );

				if( !Disjoint( lx_interval(0), res_l ) )
		      //  intersection on the left boundary
				{
					if( ge_zero(Inf(sin_2xl)) )
		         // "left" intersection
					{
						left = true;
			         //  remove the intersection by changing res_l!
						res_l = -lx_interval(1);
					}
					else
					{
						if( se_zero(Sup(sin_2xl)) )
			         // "right" intersection
						{
							right = true;
			            //  remove the intersection by changing res_l!
							res_l =  lx_interval(1);
						}
						else
			         //  zero is interior point of sin_2xl
			         //  if the real sine function has optimal precision,
			         //  this case should never happen
							both = true;
					}
				}

				if( !Disjoint( lx_interval(0), res_u ) )
		      //  intersection on the right boundary
				{
					if( ge_zero(Inf(sin_2xu)) )
		         // "left" intersection
					{
						left = true;
		            //  remove the intersection by changing res_u!
						res_u = lx_interval(1);
					}
					else
					{
						if( se_zero(Sup(sin_2xu)) )
			         // "right" intersection
						{
							right = true;
			            //  remove the intersection by changing res_u!
							res_u = -lx_interval(1);
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
	lx_interval re_tan = 1 / sinh( 2 * abs( hy ) );

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


lx_cinterval Tan(const lx_cinterval& z) throw() 
{
	lx_cinterval y;
	lx_interval
		rez = Re(z),   // rez = z.re(),
		imz = Im(z);   // imz = z.im();

	lx_real
		irez = Inf(rez),
		srez = Sup(rez),
		iimz = Inf(imz),
		simz = Sup(imz);

	lx_interval hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

	lx_real resxl, resxu, resyl, resyu;
	//
   //  1st: check for poles
	//
	if( ( !Disjoint( lx_interval(0), imz ) ) && 
		 ( !Disjoint( lx_interval(0), cos( rez ) ) ) )
		cxscthrow (STD_FKT_OUT_OF_DEF(
			"lx_cinterval tan(const lx_cinterval& z); Pole(s) in z"));
	//
   //  2nd: real part
	//
   //  evaluate tan on vertical boundaries
	//
	lx_interval
		cos_2rez   = cos( 2 * rez ), sinh_imz_2 = sqr( sinh( imz ) );

	lx_interval
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
		lx_interval
			cosh_2yl = - 1 / cosh( 2 * hyl ),
			cosh_2yu = - 1 / cosh( 2 * hyu );

		if( !Disjoint( cos_2rez, cosh_2yl ) && iimz != 0.0 )
	   //extremal curve intersects lower boundary
			horizontal_check(hyl,cosh_2yl,irez,srez,hxl,hxu,resxl,resxu);

		if( !Disjoint( cos_2rez, cosh_2yu ) && simz != 0.0 )
	   //extremal curve intersects upper boundary
			horizontal_check(hyu,cosh_2yu,irez,srez,hxl,hxu,resxl,resxu);
	}
	//
   //  3rd: imaginary part
	//
   //  evaluate tan on horizontal boundaries
	//
	lx_interval cos_rez_2 = sqr( cos( rez ) );

	lx_interval
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
	lx_interval
		cos_2xl = cos( 2 * hxl ), cos_2xu = cos( 2 * hxu );
	lx_interval im_tan;

	if( iimz < 0.0 )
   //  intersection in lower half plane?
	{
		lx_interval
			imz_h = lx_interval( iimz, min( simz, lx_real(0.0) ) ),
			cosh_2imz = - 1 / cosh( 2 * imz_h );

		if( ( !Disjoint( cosh_2imz, cos_2xl ) ) )
	   //extremal curve intersects left boundary
	   //in this case, sin( 2 * xl ) <> 0.0 (no poles here!)
		{
			im_tan = - 1 / abs( sin( 2 * hxl ) );
			resyl = min( resyl, Inf( im_tan ) );
			resyu = max( resyu, Sup( im_tan ) );
		}
		if( ( !Disjoint( cosh_2imz, cos_2xu ) ) )
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
		lx_interval
			imz_h = lx_interval( max( iimz, lx_real(0.0) ), simz ),
			cosh_2imz = - 1 / cosh( 2 * imz_h );

		if( ( !Disjoint( cosh_2imz, cos_2xl ) ) )
	   //extremal curve intersects left boundary
	   //in this case, sin( 2 * xl ) <> 0.0 (no poles here!)
		{
			im_tan = + 1 / abs( sin( 2 * hxl ) );
			resyl = min( resyl, Inf( im_tan ) );
			resyu = max( resyu, Sup( im_tan ) );
		}
		if( ( !Disjoint( cosh_2imz, cos_2xu ) ) )
	   //extremal curve intersects right boundary
	   //in this case, sin( 2 * xu ) <> 0.0 (no poles here!)
		{
			im_tan = + 1 / abs( sin( 2 * hxu ) );
			resyl = min( resyl, Inf( im_tan ) );
			resyu = max( resyu, Sup( im_tan ) );
		}
	}

	y = lx_cinterval( lx_interval(resxl,resxu ),lx_interval(resyl,resyu ) );

	return y;
} // Tan

lx_cinterval tan(const lx_cinterval& z) throw() 
{
// tan(z) has the poles z_k = pi*(1+2k)/2;   |k| in {0,1,2,3,...}.
// z = z_k + eps = pi*(1+2k)/2 + eps;  With  |eps|<<1  k can be calculated
// by:  k = Re(z)/pi - 0.5; With this k the complex value eps is given
// by:  eps = z - pi*(1+2k)/2;  pi = 3.1415926... ;
// It holds: 
// tan(z) = tan(z_k+eps) = tan[pi*(1+2k)/2 + eps] 
//        = tan[pi/2 + pi*k + eps] = tan[pi/2 + eps] = -1 / tan(eps);
// Definitions:  u = Re(eps);  u = abs(u);  v = Im(eps);  v = abs(v);
// if (Sup(u)<S && Sup(v)<S) tan(z) = -1 / Tan(eps);
// else tan(z) = Tan(z);     S = 1e-15;
// Thus, near the poles tan(z) is calculated in higher accuracy with
// -1 / Tan(eps);
// Blomquist, 27.09.2007;
		const real S = 1e-15;
		int stagsave = stagprec,
  stagmax = 39;
  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_cinterval y,eps;
  l_interval rezl,ul,vl;

  lx_interval
		  rez = Re(z), Pih;   // rez = z.re(),
  rezl = rez;
  interval re,u_,v_; 
  re = rezl;
  real x(mid(re)),k, pi(3.14159265358979323);
  lx_interval u,v; 

  int s;
  x = x/pi - 0.5;
  s = sign(x);
  k = (s>=0)? cutint(x+0.5) : cutint(x-0.5);
  if (k == 9007199254740992.0)
	  cxscthrow (STD_FKT_OUT_OF_DEF(
					 "lx_cinterval tan(const lx_cinterval& z); z out of range"));
  Pih = Pi_lx_interval();
  rez = k*Pih;
  times2pown(Pih,-1); // Pih = pi/2
  rez = rez + Pih;    // rez: inclusion of   k*pi + pi/2
  eps = z - rez;

  u = Re(eps);  u = abs(u);
  v = Im(eps);  v = abs(v);
  ul = u;       vl = v;
  u_ = ul;      v_ = vl;
  y = (Sup(u_)<S && Sup(v_)<S)? -lx_cinterval(1) / Tan(eps) : Tan(z);

  stagprec = stagsave;
  y = adjust(y);

  return y;
	} // tan()

	lx_cinterval cot(const lx_cinterval& z) throw() 
	{
// cot(z) has the poles z_k = k*pi;   |k| in {0,1,2,3,...}.
// z = z_k + eps = k*pi + eps;  With  |eps|<<1  k can be calculated
// by:  k = Re(z)/pi;   With this k the komplex value eps is given
// by:  eps = z - k*pi;    pi = 3.1415926... ;
// It holds: 
// cot(z) = cot(z_k+eps) = cot(k*pi + eps) 
//        = cot(eps) = 1 / tan(eps);
// Definitions:  u = Re(eps);  u = abs(u);  v = Im(eps);  v = abs(v);  
// if (Sup(u)<S && Sup(v)<S) cot(z) = 1 / tan(eps);
// else cot(z) = tan(pi/2 - z);     S = 1e-15;
// Thus, near the poles cot(z) is calculated in higher accuracy with
// 1 / tan(eps);
// Blomquist, 28.09.2007;
		const real S = 1e-15;
		int stagsave = stagprec,
  stagmax = 39; 
  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_cinterval y,eps;
  l_interval rezl,ul,vl;

  lx_interval 	rez = Re(z);
  rezl = rez;
  interval re,u_,v_;
  re = rezl;
  real x(mid(re)), k, pi(3.14159265358979323);
  lx_interval u,v;

  int s;
  x = x/pi;
  s = sign(x);
  k = (s>=0)? cutint(x+0.5) : cutint(x-0.5);
  if (k == 9007199254740992.0)
	  cxscthrow (STD_FKT_OUT_OF_DEF(
					 "lx_cinterval cot(const lx_cinterval& z); z out of range"));
  eps = z - Pi_lx_interval()*k;
  u = Re(eps);  u = abs(u);
  v = Im(eps);  v = abs(v);
  ul = u;       vl = v;
  u_ = ul;      v_ = vl;
  if (Sup(u)<S && Sup(v)<S)
	  y = 1 / Tan(eps);
  else 
  {
	  u = Pi_lx_interval();
	  times2pown(u,-1); // u = pi/2
	  y = Tan(u - z);
  }

  stagprec = stagsave;
  y = adjust(y);

  return y;
	} // cot

//---------------------------------- tanh -----------------------------------
//
//  tanh( z ) = transp( i * tan( transp( i * z ) )
//
lx_cinterval tanh(const lx_cinterval& z) throw()
{
	int stagsave = stagprec,
  		stagmax = 39; 
  	if (stagprec > stagmax) stagprec = stagmax;

  	lx_cinterval res = tan( lx_cinterval( Im(z), Re(z) ) ),y;
	y = lx_cinterval( Im(res), Re(res) );

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

lx_cinterval coth(const lx_cinterval& z) throw()
{ // coth( z ) = i * cot( i * z );
	int stagsave = stagprec,
  		stagmax = 39; 
  	if (stagprec > stagmax) stagprec = stagmax;

  	lx_cinterval zh = lx_cinterval( -Im(z), Re(z) ); //  zh = i*z;
  	lx_cinterval res = cot(zh);
  	zh = lx_cinterval( -Im(res), Re(res) );

  	stagprec = stagsave;
  	zh = adjust(zh);

  	return zh;
}
//
//-- end coth -----------------------------------------------------------------

//-----------------------------------------------------------------------------
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

//------------------------------------ asin -----------------------------------
//
//  Analytic inverse sine function
//
lx_interval f_aux_asin(const lx_interval& x, const lx_interval& y) 
//
//  auxiliary function for evaluating the imaginary part of asin(z);
//  f_aux_asin(z) = ( |z+1| + |z-1| ) / 2;  z = x + i*y;
//                = ( sqrt[(x+1)^2 + y^2] + sqrt[(x-1)^2 + y^2] )/2; 
//  Blomquist, 01.10.2007;        
	{
		lx_interval res;
//  interval sqr_y = sqr( y );
//  interval res = ( sqrt( sqr( x + 1.0 ) + sqr_y ) + 
//                   sqrt( sqr( x - 1.0 ) + sqr_y ) ) / 2.0;
		res = abs(x);
		if (y != 0.0 || Inf(res) < 1.0) 
		{ 
			res = sqrtx2y2(x+1.0,y) + sqrtx2y2(x-1.0,y);  // Blomquist;
			times2pown(res,-1);
		}

  //  Now we correct a possible overestimation of the lower bound
  //  of res.
  //  It holds: 
  //  (sqrt[(x+1)^2 + y^2] + sqrt[(x-1)^2 + y^2])/2 >= Max(1,|x|);
		lx_real hlb = max( lx_real(1.0), abs( Sup(x) ) );
		if( Inf(res) < hlb ) //  this is an invalid overestimation!
			res = lx_interval( hlb, Sup(res) );
		return res;
	}

	lx_interval ACOSH_f_aux(const lx_interval& x, const lx_interval& y)
// Function for calculating the imaginary part of asin(z),
// z = x + i*y;
// Notations:
// f_aux_asin(x,y) = alpha(x,y);
// alpha(x,y) := ( sqrt[(x+1)^2 + y^2] + sqrt[(x-1)^2 + y^2] )/2;
// return values:
// 1.    res = acosh( f_aux_asin(x,y) ), if  |x|>2 or |y|>2;
// 2.    res = ln[ alpha + sqrt(alpha^2 - 1) ]; d = alpha(x,y) -1
//           = ln[1 + sqrt(d)*(sqrt(d) + sqrt(2+d))], if it holds:
//       |x|<=2 und |y|<=2;   x,y are point intervals only!       
// Blomquist, 02.06.2008;

	{
		lx_interval res, xa(abs(x)), ya(abs(y)), S1, S2, S3;
		lx_real rx((Inf(xa))), ry(Inf(ya));

		if (rx>2.0 || ry>2.0) res = acosh( f_aux_asin(x,y) );
		else 
		{ //  Now with improvements:
			if (rx == 1.0) 
			{
				res = ya/2;   S2 = res/2;
				S1 = ya*(0.5 + S2/(sqrt1px2(res)+1));  // S1 = delta
				res = lnp1(S1 + sqrt(S1*(2+S1)));
			}
			else 
				if (rx<1.0) // 0 <= x < +1;
			{
				S1 = 1+xa;
				S1 = 1/(sqrtx2y2(S1,ya) + S1);
				S2 = 1-xa;
				S2 = 1/(sqrtx2y2(S2,ya) + S2);
				S1 = S1 + S2;
				times2pown(S1,-1);   // S1 = {...}/2
				S2 = ya * sqrt(S1);  // S2 = sqrt(delta)
				res = lnp1( S2*(S2 + sqrt(2+sqr(ya)*S1)) );
			}
			else  // 1 < x <= 2;
			{
				if (ya == 0)
				{
					S1 = xa-1;
					if (eq_zero(Inf(S1)))
						S1 = lx_interval(Inf(lx_interval(-2097,l_interval(1))),
											  Sup(S1));
				}
				else // 1 < x <= 2  and  0 < |y| <= 2
				{
					S1 = xa-1;
					if (eq_zero(Inf(S1)))
						S1 = lx_interval(Inf(lx_interval(-2097,l_interval(1))),
											  Sup(S1));
					S2 = ya;
					times2pown(S2,1048);
					if (Inf(S1)>Inf(S2))
						S1 *= lx_interval( lx_real(1),Sup(One_p_lx_interval()) );
					else
					{
						res = sqr(ya);
						S3 = res / (sqrtx2y2(S1,ya) + S1);  // S3 = v;
						S2 = 1+xa;
						S2 = res / (sqrtx2y2(S2,ya) + S2);  // S2 = u;
						times2pown(S1,1);                   // S1 = w=2*(x-1)
						S1 = (S3 + S2) + S1;            // S1 = u+v+w
						times2pown(S1,-1);              // S1 = delta
					}
				}
				res = lnp1( S1+sqrt(S1*(2+S1)) );
			}
		}
		return res;
	} // ACOSH_f_aux

	lx_interval Asin_beta(const lx_interval& x, const lx_interval& y )
// Calculating the improved real part of asin(z); Blomquist 22.01.2007;
// Re(asin(z)) = asin[ 2x/(sqrt[(x+1)^2+y^2] + sqrt[(x-1)^2+y^2]) ]=asin[beta];
// Improvements for beta --> +1  and  beta --> -1 are necessary because of 
// nearly vertical tangents of the real asin(t)-function for |t|-->+1;
// This function should only be used for POINT intervals x,y!  z = x + i*y;
// Blomquist, 30.05.2008;
	{
		const real c1 = 0.75;
		bool neg_x, fertig(false); 
		lx_real Infxa(Inf(x)), L_y(Inf(y));
		int ex_xl(expo_gr(lr_part(Infxa))),
					 ex_yl(expo_gr(lr_part(L_y)));
					 real ex_x(expo(Infxa)), ex_y(expo(L_y));
					 lx_interval res,beta,abs_beta,delta,w_delta,xa,Ne,Kla;

					 if (ex_xl<-1000000)
						 res = 0;
					 else
					 {  // x != 0
						 if ( ex_yl>-1000000 && ex_y>=9007199254739972.0-ex_yl && 
												ex_x <= 2097-ex_xl )
						 {   // It holds:  |x/y| <= 2^(-9007199254737871)  and  |y| != 0
            // res in [0,eps]  or  res in [-eps,0], with
            // eps = 2^(-9007199254737870); 
							 xa = x;
							 neg_x = Inf(x)<0;
	    if (neg_x) xa = -xa; // Inf(xa) >0 :
		 res = lx_interval( lx_real(0), lx_real(-9007199254737870.0,l_real(1)) );
		 if (neg_x)
			 res = -res;
						 }
						 else
						 {
							 Kla = sqrtx2y2(1-x,y);
							 Ne = sqrtx2y2(1+x,y) + Kla;
							 beta = x / ( Ne/2 );
							 if (Inf(beta)<-1) Inf(beta) = -1;
							 if (Sup(beta)> 1) Sup(beta) = 1; 
							 abs_beta = abs(beta);
							 if (Inf(abs_beta) < c1) 
								 res = asin(beta); // Normal calculation 
							 else 
							 { // Inf(abs_beta)>=c1; Calculation now with improvements:
								 Ne = Ne + 2.0; 
                // Ne = 2 + sqrt[(1+x)^2 + y^2] + sqrt[(1-x)^2 + y^2];
								 xa = x;
								 neg_x = Inf(x)<0;
		if (neg_x) xa = -xa; // Inf(xa) >0 :
		Infxa = Inf(xa);
		if (Infxa > 1.0) 
		{ // x > +1;
			if (y == 0.0) 
				fertig = true;  // due to delta == 0;
			else
			{
				beta = xa - 1;
				delta = sqrt(beta);
				delta = lx_interval(lx_real(-2100,7.4699)) * delta;
                        // delta = 7.4699*(2^-2100)*sqrt(x-1);
				if (Sup(abs(y)) < Inf(delta)) 
					fertig = true;
				else
				{
					delta = sqr(y/beta); // delta = (y/(x-1))^2
					delta = beta * sqrtp1m1(delta);
					times2pown(Ne,-1);
					delta /= Ne;
				}
			}
		} 
		else
			if (Infxa == 1.0) 
		{ // x = 1;
			delta = abs(y);
			times2pown(delta,1);
			delta /= Ne;
		} 
		else 
		{ // 0.75 <= x < 1;
			if (y == 0.0)
				delta = 1 - xa; 
			else // 0.75 <= x < 1  and  y != 0:
			{
				beta = 1 - xa; // beta > 0;
				delta = Kla + beta;
				times2pown(delta,1);
				delta /= (2+Ne);
			}
		}
		res = Pi_lx_interval();
		times2pown(res,-1); // res = pi/2;
		if (!fertig)
			res -= asin( sqrt(delta*(2-delta)) ); 
		if (neg_x) res = -res;
							 }
						 }
					 } // x != 0

					 return res;
	} // Asin_beta(...)

lx_cinterval asin(const lx_cinterval& z) throw() 
{
	int stagsave = stagprec,
  		stagmax = 30;
  	if (stagprec>stagmax) stagprec = stagmax;

  	lx_cinterval res;
  	lx_interval rez = Re(z), imz = Im(z);

	lx_real 
		irez = Inf(rez),
		srez = Sup(rez),
		iimz = Inf(imz),
		simz = Sup(imz);

	lx_interval hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

	lx_real resxl, resxu, resyl, resyu;

	bool bl = (sm_zero(iimz)) && (gr_zero(simz)),
		  raxis = (eq_zero(iimz)) && (eq_zero(simz));

	//
   //  1st: check for singularities
	//
	if ( (irez<-1 && (bl || (sm_zero(iimz) && eq_zero(simz)))) || 
		(srez >1 && (bl || (eq_zero(iimz) && gr_zero(simz)))) )
		cxscthrow(STD_FKT_OUT_OF_DEF(
		"lx_cinterval asin(const lx_cinterval& z); z contains singularities."));

	//
   //  2nd: real part
	//
	if (sm_zero(iimz) && gr_zero(simz))
   //  z intersects [-1,1]
	{
		if (se_zero(irez))
			resxl = Inf( asin(hxl) );
		else
			resxl = Inf(Asin_beta(hxl,lx_interval( max(-iimz,simz) )) );
		if (sm_zero(srez))
			resxu = Sup(Asin_beta(hxu,lx_interval( max(-iimz,simz) )) );
		else
			resxu = Sup( asin(hxu) );
	}
	else
	{
		if ( (ge_zero(iimz) && ge_zero(irez)) || 
			  (se_zero(simz) && se_zero(irez)) )
	   //  left boundary in quadrants I or III
	   //  min( Re( z ) ) in upper left corner
			resxl = Inf( Asin_beta(hxl,hyu) ); // Blomquist, 19.06.2005;
		else
	   //  left boundary in quadrants II or IV
	   //  min( Re( z ) ) in lower left corner
			resxl = Inf( Asin_beta(hxl,hyl) ); // Blomquist, 19.06.2005;
		if ( (ge_zero(iimz) && ge_zero(srez)) || (se_zero(simz) && se_zero(srez) ) )
	   //  right boundary in quadrants I or III
	   //  max( Re( z ) ) in lower right corner
			resxu = Sup( Asin_beta(hxu,hyl) ); // Blomquist, 19.06.2005;
		else
	//  right boundary in quadrants II or IV
	//  max( Re( z ) ) in upper right corner
			resxu = Sup( Asin_beta(hxu,hyu) ); // Blomquist, 19.06.2005;
	}

	//
   //  3rd: imaginary part
	//
	if (raxis) 
	{ // Interval argument is now a subset of the real axis.
     // Blomquist, 16.06.2005;
		if (sm_zero(srez)) resyl =  Inf( ACOSH_f_aux( hxu, hyu ));
		else resyl = -Sup( ACOSH_f_aux( hxu, hyu ));
		if (gr_zero(irez)) resyu = -Inf( ACOSH_f_aux( hxl, hyu ));
		else resyu =  Sup( ACOSH_f_aux( hxl, hyu ));
	} 
	else 
		if (se_zero(simz))
      //  z in lower half plane
      //  min( Im( z ) ) in point with max |z|
      //  max( Im( z ) ) in point with min |z|
		{
			if (irez < -srez)
	      //  most of z in quadrant III
			{
				resyl = -Sup( ACOSH_f_aux( hxl, hyl ) );
				if (sm_zero(srez))
					resyu = -Inf( ACOSH_f_aux( hxu, hyu ) );
				else
					resyu = -Inf( ACOSH_f_aux( lx_interval(0),hyu ) );
			}
			else
	      //  most of z in quadrant IV
			{
				resyl = -Sup( ACOSH_f_aux( hxu, hyl ) );
				if (gr_zero(irez))
					resyu = -Inf( ACOSH_f_aux( hxl, hyu ) );
				else
					resyu = -Inf( ACOSH_f_aux( lx_interval(0),hyu ) );
			}
		}
		else 
			if (ge_zero(iimz))
         //  z in upper half plane
         //  min( Im( z ) ) in point with min |z|
         //  max( Im( z ) ) in point with max |z|
			{   
				if( irez < -srez )  // if( irez + srez < 0.0 )
	         //  most of z in quadrant II
				{
					resyu = Sup( ACOSH_f_aux( hxl, hyu ) );
					if (sm_zero(srez))
						resyl = Inf( ACOSH_f_aux( hxu, hyl ) );
					else
						resyl = Inf( ACOSH_f_aux( lx_interval(0), hyl ) );
				}
				else
	         //  most of z in quadrant I
				{
					resyu = Sup( ACOSH_f_aux( hxu, hyu ) );
					if (gr_zero(irez))
						resyl = Inf( ACOSH_f_aux( hxl, hyl ) );
					else
						resyl = Inf( ACOSH_f_aux( lx_interval(0), hyl ) );
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

		res = lx_cinterval( lx_interval(resxl,resxu),lx_interval(resyl,resyu) );
		stagprec = stagsave;
		res = adjust(res);
	return res;
}
//
//-- end asin -----------------------------------------------------------------

//-----------------------------------------------------------------------------
//----------------------------- acos ------------------------------------------

lx_interval BETA_xy(const lx_interval& x, const lx_interval& y)
// Calculation of beta = x / ( sqrt[(x+1)^2 + y^2] + sqrt[(x-1)^2 + y^2]/2 );
// Only for the internal use for calculating acos(beta);
	{
		lx_interval res,xa;
		l_interval xl,yl;
		real exx,exy;
		int exxl,exyl;
		interval z;
		bool neg;

		xa = x;
		xl = li_part(xa);  yl = li_part(y);
		neg = Inf(xl)<0;
		if (neg) xa = -xa;
		exx = expo(xa);  exxl = expo_gr(Inf(xl));
		exy = expo(y);   exyl = expo_gr(Inf(yl));

		if (exyl<-100000) // y = 0;
			if (Inf(xa)<1) res = xa;
		else res = 1.0;
		else // y != 0;
			if (exxl<-100000) // x = 0
				res = 0;
		else 
		{
			z = interval(exy) + (exyl-exxl-2103);
			if (exx<Inf(z)) res = 0;
			else 
				res = xa / f_aux_asin(xa,y);
		}
		if (Sup(res)>1) res = 1.0;

		if (neg) res = -res;
		return res;
	} // BETA_xy(...)

lx_interval Asin_arg(const lx_interval& x, const lx_interval& y )
// Asin_arg calculats for point intervals x,y with 
// beta := 2*|x| / ( sqrt((x+1)^2+y^2) + sqrt((x-1)^2+y^2) )
// and delta := 1 - |beta|  an inclusion of:  
//                arcsin[ sqrt(delta)*sqrt(2-delta) ].
// The point interval x may be negative!
// Blomquist, 06.06.2008;
	{
		const real c2 = -9007199254739968.0;
		lx_interval res,xa,F_xa,S1,S2,xm1,xp1,w_delta,T;
		lx_real Infxa;
		bool neg_x;
		xa = x;
		neg_x = Inf(x)<0;
    if (neg_x) xa = -xa; // Inf(xa) > 0 :
	 Infxa = Inf(xa);

	 if (Infxa > 1.0) 
	 { // x > +1;
		 if (y == 0.0) 
			 res = 0.0; 
		 else // y != 0;
		 {
			 F_xa = xa;
			 times2pown(F_xa,c2);

			 if (Inf(abs(y)) > Inf(F_xa))
            // y > x * 2^(-9007199254739968):
            // Now the calculation of w_delta according to (1) is 
            // possible and additional the calculation of the total
            // argument of the asin-function:
			 {
				 xm1 = xa-1.0;
				 if (Inf(xm1)<=0)
					 xm1 = lx_interval(Inf( lx_interval(-2097,l_interval(1)) ),
											 Sup(xm1));
				 xp1 = xa+1.0;
				 S1 = sqrtx2y2(xm1,y);
				 S2 = sqrtx2y2(xp1,y);
				 w_delta = sqrt(2 + S1 + S2);
				 w_delta = w_delta * sqrt(S1+xm1);
				 w_delta = Sqrt2_lx_interval()*abs(y) / w_delta;
				 T = (Sqrt2_lx_interval()-w_delta)*(Sqrt2_lx_interval()+w_delta);
				 T = w_delta * sqrt(T);
				 res = asin( T );
			 }
			 else 
			 {
                // y <= x * 2^(-9007199254739968)
				 if (Inf(xa) >= 2)
				 {     
					 real exx,exxl,exy,exyl;
					 interval z;
					 double dbl;

					 exx  = expo(xa);  exy = expo(y);
					 exxl = expo_gr(li_part(xa));
					 exyl = expo_gr(li_part(y));

					 z = interval(exy)-interval(exx) + (exyl-exxl+1079);
                    // 1079 = 1074 + 5; see documentation
					 exx = Sup(z);
					 if (exx<-Max_Int_R)
						 exyl = -Max_Int_R;
					 else // Sup(z) >= -9007199254740991
					 {
						 if (diam(z)==0) exyl = Sup(z);
						 else
						 {
							 dbl = floor(_double(exx));
							 exyl = real(dbl) + 1; // exyl is integer!
						 }
					 }

					 res = lx_interval(lx_real(0.0),
											 lx_real(exyl,l_real(minreal)));
				 }
				 else // 1 < Inf(xa) < 2; 
				 {
					 T = sqrt(xa*(xa-1.0));
					 if (Inf(T) <= 0)
					 {
						 times2pown(xa,-2097);
						 xa = sqrt(xa);
						 times2pown(xa,3000);
						 res = abs(y);
						 times2pown(res,3001);
						 res = res / xa;
					 }
					 else
					 {
						 res = abs(y);
						 times2pown(res,1);
						 res = res/T;
					 }
					 res = lx_interval(lx_real(0.0),Sup(res));
				 }
			 }
		 }
	 } 
	 else
		 if (Infxa == 1.0) 
	 { // x = 1;
		 T = abs(y);
		 res = 2 + T + sqrtx2y2(lx_interval(0,l_interval(2)),T);
		 times2pown(T,1);
		 T = T / res;
		 res = asin(sqrt(T*(2-T)));
	 } 
	 else 
	 { // 0.75 <= x < 1;
		 if (y == 0.0)
		 {
			 T = sqrt1mx2(xa);
			 res = asin(T);
		 }
		 else
		 {
			 T = 1 - xa;
			 S1 = sqrtx2y2(T,y);
			 S2 = sqrtx2y2(x+1,y);
			 res = S1 + T;
			 times2pown(res,1);
			 T = 2 + S1 + S2;
			 T = res / T;
			 T = sqrt( T*(2-T) );
			 res = asin(T);
		 }
	 }
	 return res;

	} // Asin_arg

lx_interval Acos_beta(const lx_interval& x, const lx_interval& y)
// Calculating the improved real part of acos(z); Blomquist 05.06.2005;
// Re(acos(z)) = acos[ 2x / (sqrt[(x+1)^2+y^2] + sqrt[(x-1)^2+y^2]) ]
	{
		const real c1 = 0.75;
		lx_interval res(0),beta;
		beta = BETA_xy(x,y);
		if (Sup(beta)<c1) 
			if (Sup(beta)<-c1) 
		{ // Improvement for beta --> -1
			res = Pi_lx_interval() - Asin_arg(x,y);
		} 
		else res = acos(beta); // Normal calculation
		else  // Sup(beta)>=c1
			res = Asin_arg(x,y);
		return res;
} // Acos_beta(...)

//
lx_cinterval acos(const lx_cinterval& z) throw()   
{
	lx_interval 
		rez = Re(z),
		imz = Im(z);

	lx_real
		irez = Inf(rez),
		srez = Sup(rez),
		iimz = Inf(imz),
		simz = Sup(imz);

	lx_interval
		hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

	bool bl = iimz< 0.0 && simz>0.0,
		raxis = iimz==0.0 && simz==0;
	lx_real resxl, resxu, resyl, resyu;
	//
   //  1st: check for singularities
	//
	if( (irez<-1 && (bl || (iimz<0  && simz==0))) || 
		 (srez >1 && (bl || (iimz==0 && simz>0))) )
		cxscthrow(STD_FKT_OUT_OF_DEF("lx_cinterval acos(const lx_cinterval& z); z contains singularities."));
						//
   //  2nd: real part
						//
   //  Blomquist, 05.02.2007;

	if( iimz < 0.0 && simz > 0.0 )
   //  z intersects [-1,1] on the x-axis
	{
		if( irez <= 0.0 ) resxu = Sup( acos( hxl ) );
		else resxu = Sup( Acos_beta(hxl,lx_interval( max(-iimz,simz) )) );

		if( srez < 0.0 ) 
			resxl = Inf( Acos_beta(hxu,lx_interval(max(-iimz,simz))) );
		else resxl = Inf( acos( hxu ) ); 
	}
	else
	{
		if (irez<0 && srez>0) 
	   // z intersects the posizive or negative y-axis
			if (ge_zero(iimz)) 
		   {
				resxl = Inf( Acos_beta(hxu,hyl) );
				resxu = Sup( Acos_beta(hxl,hyl) );
			} 
			else 
			{
				resxl = Inf( Acos_beta(hxu,hyu) );
				resxu = Sup( Acos_beta(hxl,hyu) );
			}
			else 
			{
				if ( (ge_zero(iimz) && ge_zero(irez)) || 
					  (se_zero(simz) && sm_zero(irez)) )
	         //  left boundary in quadrants I or III
	         //  min( Re( z ) ) in lower right corner
					resxl = Inf( Acos_beta(hxu,hyl) );
				else
	         //  left boundary in quadrants II or IV
	         //  min( Re( z ) ) in upper right corner
					resxl = Inf( Acos_beta(hxu,hyu) );

				if ( (ge_zero(iimz) && gr_zero(srez)) || 
					  (se_zero(simz) && se_zero(srez)) )
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
	if (raxis) 
	{ // Interval argument is now a subset of the real axis.
     // Blomquist, 16.06.2005;
		if (srez<0.0) resyl =  Inf( ACOSH_f_aux( hxu, hyu )); 
		else          resyl = -Sup( ACOSH_f_aux( hxu, hyu )); 
		if (irez>0.0) resyu = -Inf( ACOSH_f_aux( hxl, hyu ));
		else          resyu =  Sup( ACOSH_f_aux( hxl, hyu ));
	} 
	else 
		if( simz <= 0.0 )
      //  z in lower half plane
      //  min( Im( z ) ) in point with max |z|
      //  max( Im( z ) ) in point with min |z|
		{
			if (irez < -srez)
	      //  most of z in quadrant III
			{
				resyl = -Sup( ACOSH_f_aux( hxl, hyl ) );
				if( srez < 0.0 )
					resyu = -Inf( ACOSH_f_aux( hxu, hyu ) );
				else
					resyu = -Inf( ACOSH_f_aux( lx_interval(0), hyu ) );
			}
			else
         //  most of z in quadrant IV
			{
				resyl = -Sup( ACOSH_f_aux( hxu, hyl ) );
				if( irez > 0.0 )
					resyu = -Inf( ACOSH_f_aux( hxl, hyu ) );
				else
					resyu = -Inf( ACOSH_f_aux( lx_interval(0), hyu ) );
			}
		}
		else 
			if( ge_zero(iimz) )
         //  z in upper half plane
         //  min( Im( z ) ) in point with min |z|
         //  max( Im( z ) ) in point with max |z|
			{   
				if( irez < -srez )  // if( irez + srez < 0.0 )
	         //  most of z in quadrant II
				{
					resyu = Sup( ACOSH_f_aux( hxl, hyu ) );
					if( sm_zero(srez) )
						resyl = Inf( ACOSH_f_aux( hxu, hyl ) );
					else
						resyl = Inf( ACOSH_f_aux( lx_interval(0), hyl ) );
				}
				else
	         //  most of z in quadrant I
				{
					resyu = Sup( ACOSH_f_aux( hxu, hyu ) );
					if( irez > 0.0 )
						resyl = Inf( ACOSH_f_aux( hxl, hyl ) );
					else
						resyl = Inf( ACOSH_f_aux( lx_interval(0), hyl ) );
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

	return lx_cinterval( lx_interval( resxl, resxu ),-lx_interval( resyl, resyu ) );

} // acos(...)

//-- end acos -----------------------------------------------------------------

//-- asinh --------------------------------------------------------------------
//
lx_cinterval asinh( const lx_cinterval& z ) throw()
//
//  asinh( z ) = i * asin( -i * z )
//
{
	lx_cinterval res = asin( lx_cinterval( Im(z), -Re(z) ) );
	return lx_cinterval( -Im(res), Re(res) );
}
//
//-- end asinh ----------------------------------------------------------------

//-- acosh --------------------------------------------------------------------
//
lx_cinterval acosh( const lx_cinterval& z ) throw()
//
//  acosh( z ) = i * acos( z ) = +/- i * ( pi / 2 - asin( z ) )
//
{
	lx_interval
		rez = Re(z),
		imz = Im(z);

	lx_real
		irez = Inf(rez),
		srez = Sup(rez),
		iimz = Inf(imz),
		simz = Sup(imz);

	lx_interval hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

	lx_real resxl, resxu, resyl, resyu;

	//
   //  1st: check for singularities
	//
	if( ( se_zero(iimz) && ge_zero(simz) ) && ( irez < 1.0 ) )
		cxscthrow(STD_FKT_OUT_OF_DEF(
		"lx_cinterval acosh( const lx_cinterval& z ); z contains singularities."));
   //  With this restriction the complex interval argument and the real axis 
	//  must not have any common point, if irez < +1;
   //  So for example the negative real axis must not be touched from above if
	//  irez<1, although this should be possible if the principal branch is
	//  considered! So the above restriction is too widely in
   //  some cases;  Blomquist, 21.06.2005;
	//
   //  2nd: z in upper half plane (or on the real axis)
   //  acosh( z ) =  + i * ( pi / 2 - asin( z ) )
	//
	if( gr_zero(iimz) )
	{
		lx_cinterval res = acos(z);
		return lx_cinterval( -Im(res),Re(res) );
	}
	//
   //  3rd: z in lower half plane
   //  acosh( z ) =  - i * ( pi / 2 - asin( z ) )
	//
	if( sm_zero(simz) )
	{
   // cinterval res = HALFPI() - asin( z );
		lx_cinterval res = acos(z);  // Blomquist, 14.06.2005
		return lx_cinterval( Im(res), -Re(res) );
	}
	//
   //  z intersects [1,infinity)
	//
   //  real part
   //  minimum on the left on real axes, maximum in lower or upper right corner
	//
	resxl = Inf( acosh( hxl ) );
	lx_interval ytilde( max( -iimz, simz ) );
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

	return lx_cinterval(lx_interval( resxl, resxu ),lx_interval( resyl, resyu ));
}
//
//-- end acosh ----------------------------------------------------------------


//-- atan ---------------------------------------------------------------------
//

bool sign_test(const lx_interval& x, int s_org)
// Only for internal use by function re_atan(...)
{
	bool bl1,bl2,alter;
	if (diam(li_part(x))>0)
	{
		bl1 = eq_zero(Sup(x)) && (s_org==1);
		bl2 = eq_zero(Inf(x)) && (s_org==-1);
		alter = bl1 || bl2;
	}
	else
		alter = sign(Sup(li_part(x))) != s_org;
	return alter;
} // sign_test(...)

void re_atan( const lx_interval& y, lx_interval& x, lx_interval& res )
// Calculating an inclusion res of 1 - y^2 - x^2;
// x is always a point interval and y =[y1,y2], with y1 <= y2;
// Blomquist, 10.06.2008;
	{
		const real c1 = 4503599627369982.0; 
		lx_interval ya(abs(y)), One(lx_interval(0,l_interval(1))), xa;
		lx_real R,U;
		real ex_x,ex_xl,ex_y,ex_yl,s;
		int signx;
		bool alter;

		ex_x = expo(x);
		ex_xl = expo_gr(Inf(li_part(x)));
		ex_y = expo(y);
		ex_yl = expo_gr(Sup(li_part(ya)));
		ya = y;
		signx = sign(Inf(li_part(x)));

		if (ex_xl < -1000000) 
			res = 1.0;
    else // x <> 0:
		 if (ex_yl < -1000000) // y == [0,0] and x <> [0,0];
			 if(ex_x > c1)
	 {
		 res = 1/x - x;
		 x = -1;  // Eigentlich sollte hier x = +1 stehen ???
	 }
	 else  // y == [0,0] and x^2 without overflow
	 {
		 res = (1-x)*(1+x);
		 if (Sup(res)>1) // to avoid an overestimation
			 res = lx_interval(Inf(res),lx_real(1.0));
	 }
	 else // y <> [0,0] and x <> [0,0];
		 if (ex_x>c1) // |x| --> infty;
			 if (ex_y>c1) // |x| --> infty and |y| --> infty;
	 {
		 s = c1 - max(ex_x,ex_y); // s < 0; and s is an integer value!
		 times2pown_neg(One,s); // 2^s
		 times2pown_neg(x,s);   // x*2^s
		 times2pown(ya,s);      // y*2^s
		 res = (One-x)*(One+x) - sqr(ya);
		 times2pown_neg(x,s);   // x*2^2s
	 } 
	 else // |x| --> infty  and  |y| not--> infty;
	 {
		 res = sqr(y);  xa = abs(x);
		 if (res == 1.0)
			 res = xa;
		 else
		 {
			 res = res - 1.0; // res = y^2 - 1;
			 R = Sup(abs(res));
			 ex_y = expo(R);
			 ex_yl = expo_gr(lr_part(R));
			 ex_xl = ex_xl - 1051;
			 if (ex_y+ex_yl < 2*(ex_x+ex_xl))
				 res = xa * 
						 lx_interval(Inf(One_m_lx_interval()),
										 Sup(One_p_lx_interval()));
			 else 
				 res = res/xa + xa;
		 }
		 x = (Sup(x)>0)? -1.0 : 1.0;
	 }
	 else // |x| not to infty
		 if (ex_y>c1) // |y| --> infty;
	 {
		 s = c1 - ex_y; // 
		 res = (x-1)*(x+1);
		 R = Sup( abs(res) );
		 times2pown_neg(x,s);   // x*2^s
		 times2pown(ya,s);      // y*2^s
		 if (eq_zero(R))
			 res = sqr(ya);
		 else
		 {   // R > 0;
			 ya = sqr(ya); // ya = (y*2^s)^2;
			 times2pown_neg(res,s);
			 times2pown_neg(res,s);
                        // res = {(x-1)*(x+1)}*2^2s;
			 res = res + ya;
		 }
		 times2pown_neg(x,s); // (x*2^s)*2^s = x*2^(2s)
		 x = -x;
	 }
	 else // |x|,|y| not to +infty
		 if (ex_y<-c1) // |y| ---> 0, y<>[0,0];
			 if (abs(x)==1)
	 {
		 times2pown(x,9007199254738894.0);
		 x = -x;
		 times2pown(ya,4503599627369447.0);
		 res = sqr(ya);
	 }
	 else
	 {
		 res = (x-1)*(x+1) + sqr(ya);
		 x = -x;
	 }
	 else // |x|,|y| not to +infty and |y| not to 0 
		 if (ex_x < -c1) // now: |x| --> 0
	 {
		 res = (y-1)*(y+1);
		 if (res == 0.0) // alpha = -1/x;
		 {
			 x = -1;
			 res = x;
		 }
		 else
		 {
			 res += sqr(x);
			 x = -x;
		 }
	 }
	 else // x^2 and y2 can now be evaluated without any
                             // integer overflow!
		 res = (1-x)*(1+x) - sqr(y);
	 alter = sign_test(x,signx);
	 if (alter)
	 {
		 x = -x;     
		 res = -res; 
	 }
} // re_atan

void re_vert( const lx_real& x, const lx_interval& hx,
				  const lx_real& rew_inf, const lx_real& rew_sup,
		        lx_real& resxl, lx_real& resxu )
//
//  Subroutine of analytic inverse tangent function.
//  Evaluate real part on a vertical boundary.
//
{
	if( eq_zero(x) )
   //  singularities have been handled before, hence Re( w ) > 0
	{
		resxl = 0.0;
		resxu = 0.0;
	}
	else
	{
		lx_interval hx2(hx),Pid2,Pid4;
		Pid4 = Pid4_lx_interval();
		times2pown(hx2,1);
		if( gr_zero(x) )
	   //  w in quadrants I and/or II
	   //  atan is the inverse function of tan(t), t in (-pi/2,pi/2).
		{
			resxl = gr_zero(rew_sup)? Inf( Atan( hx2,rew_sup )/2.0 )
						: ( sm_zero(rew_sup)? Inf( (Atan( hx2,rew_sup ) + Pi_lx_interval() )/2.0 )
					: Inf( Pid4 ) );

			resxu = gr_zero(rew_inf)? Sup( Atan( hx2,rew_inf )/2.0 )
					: ( sm_zero(rew_inf)? Sup( (Atan( hx2,rew_inf ) + Pi_lx_interval())/2.0 )
					: Sup( Pid4 ) );
		}
		else
	   //  w in quadrants III and/or IV
		{
			resxl = sm_zero(rew_inf)? Inf( (Atan( hx2,rew_inf ) - Pi_lx_interval())/2.0 )
					: ( gr_zero(rew_inf)? Inf( Atan( hx2,rew_inf )/2.0 )
					: -Sup( Pid4 ) );
			resxu = sm_zero(rew_sup)? Sup( (Atan( hx2,rew_sup ) - Pi_lx_interval())/2.0 )
					: ( gr_zero(rew_sup)? Sup( Atan( hx2,rew_sup )/2.0 )
					: -Inf( Pid4 ) );
		}
	}
} //  re_vert

lx_interval T_atan(const lx_real& x)
// Calculating an inclusion of 
//              ln[ 1+2/(sqrt(1+x^2)-1) ]. 
// x will be handeld as a point interval [x]=[x,x], with x>0.
// Blomquist, 10.06.2008;
{
	const real c1 = 4503599627367859.0;
	lx_interval res,
  	ix(x);      // ix is point interval with x>0;
  	real ex_x(expo(ix));

  	if (ex_x<-c1) 
  	{
	  	res = sqrt1px2(ix) + 1;
	  	times2pown(res,1);
	  	res += sqr(ix); // res = 2(1+sqrt(1+x^2)) + x^2;
	  	ix = ln(ix);
	  	times2pown(ix,1);       // ix = 2*ln(x);
	  	res = ln(res) - ix;
  	}
  	else  // -c1 <= ex_x 
	  	if (ex_x<1150) // -c1 <= ex_x < 1150
		  	res = lnp1(2/sqrtp1m1(sqr(ix)));
  		else res = lnp1(2/(sqrt1px2(ix)-1));

  	return res;
} // T_atan

lx_interval Q_atan(const lx_interval& x, const lx_interval& y)
{
// x: abs(Re(z)); So x is a real interval x = [x1,x2], with 0<=x1<x2.
// y: Inf(Im(z)); So y is a point interval, with y >= 0.
// Q_atan returns an inclusion of ln[1 + 4y/(x^2+(1-y)^2)]
// Tested in detail; Blomquist, 13.06.2008;

	const real c1 = 4503599627367859.0,
  	c2 = 9007199254740990.0,
  	c3 = 4503599627370495.0;
  const lx_real S = lx_real(-c2,MinReal);

  double Dbl;
  int r;
  lx_real R;
  lx_interval res(0.0);
  real ex_y(expo(y)), ex_yl(expo_gr(li_part(y))),
				ex_x(expo(Sup(x))),ex_xl,exx,exxl,dbl,s,up,exy,exyl;
				interval dbli,z;

  ex_xl = expo_gr(lr_part(Sup(x)));
  if (ex_xl<-100000) ex_x = 0;
  if (ex_yl>-100000) // y = [y,y], y>0;
  {
    if (y==1.0) 
    {  // y = 1;
      if (ex_x < -c1)
      {
        res = ln(x);
        times2pown(res,1); // res = 2*ln(x)
        res = ln(4+sqr(x)) - res; 
      }
      else
        if (ex_x > c1)
        {
          R = Inf(x);
          exx = expo(R);
          exxl = expo_gr(lr_part(R));
          if (514-exxl < -c3+exx)
            res = lx_interval(lx_real(0),S);
          else
            if (3 - exxl >= -c3 + exx) // dbl >= -c2
            {
              dbl = -2*(exx+exxl-3); // without any rounding!
              res = lx_interval(lx_real(0),lx_real(dbl,l_real(1)));
            }
            else // dbl < -c2
            {
              s = c2-exx;  s = (s - exx) -2*exxl + 6;
              r = (int) _double(s);
              res = lx_interval(lx_real(0),lx_real(-c2,comp(0.5,r+1)));
            }
        }
        else // x^2 now without integer overflow
          res = lnp1(sqr(2/x));
   }
   else // Now:  y>0  and  y != [1,1]  and  0<=x1<=x2;
   {
     if (ex_x>c1 || ex_y>c1 )
     {
       R = Inf(x);
       exx = expo(R); exxl = expo_gr(lr_part(R));
       R = Inf(y-1.0); // R != 0:
       exy = expo(R); exyl = expo_gr(lr_part(R));

       if (exxl<-100000) // x1 = 0;
         dbli = 2*( interval(exy) + (exyl-2) );
       else 
       {
         z = 2*( interval(exx) + (exxl-2) );
         dbli = 2*( interval(exy) + (exyl-2) );
         if (Sup(z) > Sup(dbli)) dbli = z;
       }
       dbli = (interval(exy) - dbli) + (exyl + 3);
       dbl = Sup(dbli);
       // Now it holds:  Sup{4y/(x^2+(1-y)^2)} < 2^dbl;
       // and this upper bound 2^dbl is guaranteed,
       // even if dbl is not an integer value!
       up = Inf( interval(-c2) - 1022);
       if (dbl<up)
         res = lx_interval(lx_real(0),S);
       else // dbl >= up
         if (dbl>=-c2)
         {
           if (!(Is_Integer(dbl)))
           {  // rounding to the next greater integer value.
             Dbl = floor(_double(dbl)) + 1; // Dbl: integer value
             dbl = (real) Dbl;
           }
           res = lx_interval(lx_real(0),lx_real(dbl,l_real(1)));
         }
         else
         {
           s = Sup(interval(c2) + dbl);
           if ( !(Is_Integer(s)) )
           {  // rounding to the next greater integer value.
             Dbl = floor(_double(s)) + 1;
             r = (int) Dbl;
           }
           else r = (int) _double(s);
           res = lx_interval(lx_real(0),lx_real(-c2,comp(0.5,r+1)));
         }
      }
      else // Now it holds: y>0  and  y!=1  and  0<=x1<=x2;
           // and:  x^2, (1-y)^2 produce no overflow for
           //       x2, |1-y| --> +infty;
           // furthermore |1-y| produces no integer overflow
           // for y --> +1.
      {
        res = y;
        times2pown(res,2); // res = 4*y;
        res = res/(sqr(x) + sqr(1-y));
        res = lnp1(res); // res: inclusion of Q(x,y);
      }
   } 
  }   
  return res;
} // Q_atan

lx_cinterval atan( const lx_cinterval& z ) throw()
{
	int stagsave = stagprec,
  	stagmax = 30;
  	if (stagprec>stagmax) stagprec = stagmax;

  	lx_cinterval res;
  	lx_interval
		rez = Re(z),
		imz = Im(z);

	lx_real
		irez = Inf(rez),
		srez = Sup(rez),
		iimz = Inf(imz),
		simz = Sup(imz);

	lx_interval  // all theses intervals are point intervals!
		hxl(irez), hxu(srez), hyl(iimz), hyu(simz); 

	lx_real
		resxl, resxu, resyl, resyu;
	//
   //  1st: check for singularities
	//
	if( (se_zero(irez) && ge_zero(srez)) && ( iimz <= -1.0 || simz >= 1.0 ) )
		cxscthrow(STD_FKT_OUT_OF_DEF("lx_cinterval atan( const lx_cinterval& z ); points of the branch cuts are not allowed in z."));
	//
   //  2nd: real part
   //  Re( atan( z ) ) = Arg( w ) / 2, where w = 1 - x^2 - y^2 + i * 2x )
	//
   //  evaluate atan on vertical boundaries
	//
	lx_interval
	// y_sqr = sqr( imz ),
	// rew_l = (1 - y_sqr) - sqr( hxl ),  
	// Blomquist; before: rew_l = 1 - sqr(hxl) - y_sqr, 
   //      rew_u = (1 - y_sqr) - sqr( hxu );  
	// Blomquist; before: rew_u = 1 - sqr(hxu) - y_sqr; 
		rew_l, rew_u;
//  ------------------------------ Blomquist ---------------------------------
//  ----------       Improvements for Im(z) = [1,1]  or  Im(z) = [-1,-1]  
	bool sqrImz_1 = (iimz==simz) && (iimz==1.0 || iimz==-1.0); 
	                   // Test for Im(z) = [1,1] or [-1,-1]

	if (sqrImz_1) 
	{ 
		rew_l = -abs(hxl);  hxl = lx_interval(0,l_interval(sign(irez))); 
		rew_u = -abs(hxu);  hxu = lx_interval(0,l_interval(sign(srez)));
	}
	else 
	{
	// rew_l = (1 - sqr( imz )) - sqr( hxl );
		re_atan(imz,hxl,rew_l); 
	// rew_u = (1 - sqr( imz )) - sqr( hxu );
		re_atan(imz,hxu,rew_u);
	}
//  ------------------------------ Blomquist; 22.02.05; --------------------  

	//
   //  left boundary
	//
	lx_real rew_inf = Inf( rew_l );
	lx_real rew_sup = Sup( rew_l );
	re_vert( irez, hxl, rew_inf, rew_sup, resxl, resxu );

	//
   //  right boundary
	//
	rew_inf = Inf( rew_u );
	rew_sup = Sup( rew_u );
	lx_real res_l, res_u;
	re_vert( srez, hxu, rew_inf, rew_sup, res_l, res_u );

	if (res_l<resxl) resxl = res_l;
	if (res_u>resxu) resxu = res_u;

	//
   //  look for extremal values on horizontal boundaries
   //  since atan( x+iy ) = atan( x-iy ),
   //  intersections can be considered in the upper half plane
	//
	lx_real abs_y_min = Inf( abs( imz ) );
	if( abs_y_min > 1.0 )
	{
		lx_interval
			abs_hyl = lx_interval( abs_y_min ),
         // abs_hxl = sqrt( sqr( abs_hyl ) - 1.0 );
			abs_hxl = sqrtx2m1(abs_hyl);  // Blomquist;

		if( Sup( abs_hxl ) > irez && Inf( abs_hxl ) < srez )
	   //  extremal curve intersects lower boundary of x+i|y| in quadrant I
	   //  intersection in Q I or Q IV: update minimum
	   //	resxl = inf( atan( abs_y_min / abs_hxl ) ) / 2.0;
			resxl = Inf( (Pi_lx_interval() - atan( 1.0 / abs_hxl ))/2.0 );
		else 
			if ( -Inf( abs_hxl ) > irez && -Sup( abs_hxl ) < srez )
	      //  extremal curve intersects lower boundary of x+i|y| in quadrant II
	      //  intersection in Q II or Q III: update maximum
				resxu = Sup( (atan( 1.0 / abs_hxl ) - Pi_lx_interval())/2.0 );
	}
   //  3rd: imaginary part
   //  Im( atan( z ) ) = +/- Ln( 1 +/- 4y/( x^2 + (1 -/+ y)^2 ) ) / 4
	//
   //  evaluate atan on horizontal boundaries
	lx_interval
		abs_rez = abs(rez), 
		im_atan_l, im_atan_u;

	if ( sm_zero(iimz) )
//    im_atan_l = -ln( 1 - 4 * hyl / ( x_sqr + sqr( 1 + hyl ) ) ) / 4.0;
//    im_atan_l = -lnp1(-4 * hyl / ( x_sqr + sqr( 1 + hyl ) )) / 4.0;  // Blomquist
		im_atan_l = -Q_atan(abs_rez,-hyl);  // Blomquist 
	else
//    im_atan_l = ln( 1 + 4 * hyl / ( x_sqr + sqr( 1 - hyl ) ) ) / 4.0;
		im_atan_l = Q_atan(abs_rez,hyl);  // Blomquist
	times2pown(im_atan_l,-2); // Division by 4
	if ( sm_zero(simz) )
//    im_atan_u = -ln( 1 - 4 * hyu / ( x_sqr + sqr( 1 + hyu ) ) ) / 4.0;
//    im_atan_u = -lnp1(-4 * hyu / ( x_sqr + sqr( 1 + hyu ) ) ) / 4.0; // Blomquist
		im_atan_u = -Q_atan(abs_rez,-hyu);  // Blomquist
	else
//    im_atan_u = ln( 1 + 4 * hyu / ( x_sqr + sqr( 1 - hyu ) ) ) / 4.0;
		im_atan_u = Q_atan(abs_rez,hyu);  // Blomquist
	times2pown(im_atan_u,-2); // Division by 4
	resyl = min( Inf( im_atan_l ), Inf( im_atan_u ) );
	resyu = max( Sup( im_atan_l ), Sup( im_atan_u ) );
	//
   //  look for extremal values on vertical boundaries,
   //  if vertical boundaries intersect extremal curves
	//
	lx_real abs_x_min = Inf( abs( rez ) );
	lx_interval
		x_extr = lx_interval( abs_x_min ),
//    y_extr = sqrt( 1.0 + sqr( x_extr ) );
		y_extr = sqrt1px2(x_extr); // Blomquist;

	if( Inf( y_extr ) < simz && Sup( y_extr ) > iimz )
    //  extremal curve intersects left boundary of |x|+iy in quadrant I
    //  update maximum
    //  resyu = Sup( ln( 1 + 4 * y_extr / ( sqr( x_extr ) + sqr( 1 - y_extr ) ) ) ) / 4.0;
    //	resyu = Sup( T_atan(abs_x_min)/4.0 );    // Blomquist
	{
		rez = T_atan(abs_x_min);
		times2pown(rez,-2);
		resyu = Sup(rez);
	}

	if( -Sup(y_extr) < simz && -Inf(y_extr) > iimz )
    //  extremal curve intersects left boundary of |x|+iy in quadrant IV
    //  update minimum
    //  resyl = -Sup( ln( 1 + 4 * y_extr / ( sqr( x_extr ) + sqr( 1 - y_extr ) ) ) ) / 4.0;
    //	resyl = -Sup( T_atan(abs_x_min)/4.0 );  // Blomquist
	{
		rez = T_atan(abs_x_min);
		times2pown(rez,-2);
		resyl = -Sup(rez);
	}

	res = lx_cinterval( lx_interval(resxl,resxu),lx_interval(resyl,resyu) );
	stagprec = stagsave;
	res = adjust(res);
	return res;
}
//
//-- end atan -----------------------------------------------------------------

lx_cinterval acot( const lx_cinterval& z ) throw()
{
	int stagsave = stagprec,
  	stagmax = 30;
  	if (stagprec>stagmax) stagprec = stagmax;

  	lx_cinterval res;
  	lx_interval
		rez = Re(z),
		imz = Im(z);

	lx_real
		irez = Inf(rez),
		srez = Sup(rez),
		iimz = Inf(imz),
		simz = Sup(imz);

	lx_interval  // all theses intervals are point intervals!
		hxl(irez), hxu(srez), hyl(iimz), hyu(simz);

	lx_real
		resxl, resxu, resyl, resyu;
	//
   //  1st: check for singularities
	//
	if( (se_zero(irez) && ge_zero(srez)) && ( iimz <= 1.0 || simz >= -1.0 ) )
		cxscthrow(STD_FKT_OUT_OF_DEF("lx_cinterval acot( const lx_cinterval& z ); points of the branch cuts are not allowed in z."));
	//
   //  2nd: real part
   //  Re( atan( z ) ) = Arg( w ) / 2, where w = 1 - x^2 - y^2 + i * 2x )
	//
   //  evaluate atan on vertical boundaries
	//
	lx_interval
//      y_sqr = sqr( imz ),
//      rew_l = (1 - y_sqr) - sqr( hxl ),  // Blomquist; before: rew_l = 1 - sqr(hxl) - y_sqr, 
//      rew_u = (1 - y_sqr) - sqr( hxu );  // Blomquist; before: rew_u = 1 - sqr(hxu) - y_sqr; 
		rew_l, rew_u;

//  ------------------------------ Blomquist ---------------------------------
//  ----------       Improvements for Im(z) = [1,1]  or  Im(z) = [-1,-1]  
	bool sqrImz_1 = (iimz==simz) && (iimz==1.0 || iimz==-1.0); 
	              // Test for Im(z) = [1,1] or [-1,-1]

	if (sqrImz_1) 
	{ 
		rew_l = abs(hxl);  hxl = lx_interval(0,l_interval(sign(irez))); 
		rew_u = abs(hxu);  hxu = lx_interval(0,l_interval(sign(srez)));
	}
	else 
	{
	// rew_l = (1 - sqr( imz )) - sqr( hxl );
		re_atan(imz,hxl,rew_l);
		rew_l = -rew_l; 
	// rew_u = (1 - sqr( imz )) - sqr( hxu );
		re_atan(imz,hxu,rew_u);
		rew_u = -rew_u;
	}
//  ------------------------------ Blomquist; 22.02.05; --------------------  

	//
   //  left boundary
	//
	lx_real rew_inf = Inf( rew_l );
	lx_real rew_sup = Sup( rew_l );
	re_vert( irez, hxl, rew_inf, rew_sup, resxl, resxu );

	//
   //  right boundary
	//
	rew_inf = Inf( rew_u );
	rew_sup = Sup( rew_u );
	lx_real res_l, res_u;
	re_vert( srez, hxu, rew_inf, rew_sup, res_l, res_u );

	if (res_l<resxl) resxl = res_l;
	if (res_u>resxu) resxu = res_u;

	//
   //  look for extremal values on horizontal boundaries
   //  since atan( x+iy ) = atan( x-iy ),
   //  intersections can be considered in the upper half plane
	//
	lx_real abs_y_min = Inf( abs( imz ) );

	if( abs_y_min > 1.0 )
	{
		lx_interval
			abs_hyl = lx_interval( abs_y_min ),
//      abs_hxl = sqrt( sqr( abs_hyl ) - 1.0 );
		abs_hxl = sqrtx2m1(abs_hyl);  // Blomquist;

		if( Sup( abs_hxl ) > irez && Inf( abs_hxl ) < srez )
	//  extremal curve intersects lower boundary of x+i|y| in quadrant I
	//  intersection in Q I or Q IV: update maximum
	//	resxl = inf( atan( abs_y_min / abs_hxl ) ) / 2.0;
			resxu = Sup( atan( 1.0 / abs_hxl ) / 2.0 );
		if ( -Inf( abs_hxl ) > irez && -Sup( abs_hxl ) < srez )
	//  extremal curve intersects lower boundary of x+i|y| in quadrant II
	//  intersection in Q II or Q III: update minimum
			resxl = -Sup( atan( 1.0 / abs_hxl ) /2.0 );
	}

   //  3rd: imaginary part
   //  Im( atan( z ) ) = +/- Ln( 1 +/- 4y/( x^2 + (1 -/+ y)^2 ) ) / 4
	//
   //  evaluate atan on horizontal boundaries
	lx_interval
		abs_rez = abs(rez), 
		im_atan_l, im_atan_u;

	if ( sm_zero(iimz) )
//    im_atan_l = -ln( 1 - 4 * hyl / ( x_sqr + sqr( 1 + hyl ) ) ) / 4.0;
//    im_atan_l = -lnp1(-4 * hyl / ( x_sqr + sqr( 1 + hyl ) )) / 4.0;  // Blomquist
		im_atan_l = -Q_atan(abs_rez,-hyl);  // Blomquist 
	else
//    im_atan_l = ln( 1 + 4 * hyl / ( x_sqr + sqr( 1 - hyl ) ) ) / 4.0;
		im_atan_l = Q_atan(abs_rez,hyl);  // Blomquist
	times2pown(im_atan_l,-2); // Division by 4

	if ( sm_zero(simz) )
//    im_atan_u = -ln( 1 - 4 * hyu / ( x_sqr + sqr( 1 + hyu ) ) ) / 4.0;
//    im_atan_u = -lnp1(-4 * hyu / ( x_sqr + sqr( 1 + hyu ) ) ) / 4.0; // Blomquist
		im_atan_u = -Q_atan(abs_rez,-hyu);  // Blomquist
	else
//    im_atan_u = ln( 1 + 4 * hyu / ( x_sqr + sqr( 1 - hyu ) ) ) / 4.0;
		im_atan_u = Q_atan(abs_rez,hyu);  // Blomquist
	times2pown(im_atan_u,-2); // Division by 4

	resyl = min( Inf( im_atan_l ), Inf( im_atan_u ) );
	resyu = max( Sup( im_atan_l ), Sup( im_atan_u ) );
	//
   //  look for extremal values on vertical boundaries,
   //  if vertical boundaries intersect extremal curves
	//
	lx_real abs_x_min = Inf( abs( rez ) );
	lx_interval
		x_extr = lx_interval( abs_x_min ),
//    y_extr = sqrt( 1.0 + sqr( x_extr ) );
		y_extr = sqrt1px2(x_extr);  // Blomquist;

	if( Inf( y_extr ) < simz && Sup( y_extr ) > iimz )
    //  extremal curve intersects left boundary of |x|+iy in quadrant I
    //  update maximum
    //  resyu = Sup( ln( 1 + 4 * y_extr / ( sqr( x_extr ) + sqr( 1 - y_extr ) ) ) ) / 4.0;
    //	resyu = Sup( T_atan(abs_x_min)/4.0 );    // Blomquist
	{
		rez = T_atan(abs_x_min);
		times2pown(rez,-2);
		resyu = Sup(rez);
	}

	if( -Sup(y_extr) < simz && -Inf(y_extr) > iimz )
    //  extremal curve intersects left boundary of |x|+iy in quadrant IV
    //  update minimum
    //  resyl = -Sup( ln( 1 + 4 * y_extr / ( sqr( x_extr ) + sqr( 1 - y_extr ) ) ) ) / 4.0;
    //	resyl = -Sup( T_atan(abs_x_min)/4.0 );  // Blomquist
	{
		rez = T_atan(abs_x_min);
		times2pown(rez,-2);
		resyl = -Sup(rez);
	}

	res = lx_cinterval( lx_interval(resxl,resxu),lx_interval(-resyu,-resyl) );
	stagprec = stagsave;
	res = adjust(res);
	return res;
}
//
//-- end acot -----------------------------------------------------------------

//-- atanh --------------------------------------------------------------------
//
lx_cinterval atanh( const lx_cinterval& z ) throw()
//
//  atanh( z ) = - i * atan( i * z )
//
{
	lx_cinterval res = atan( lx_cinterval( -Im(z), Re(z) ) );
	return lx_cinterval( Im(res), -Re(res) );
}
//
//-- end atanh ----------------------------------------------------------------

//-- acoth --------------------------------------------------------------------
//
lx_cinterval acoth( const lx_cinterval& z ) throw()
//
//  acoth( z ) = i * acot( i * z )
//
{
	lx_cinterval res = acot( lx_cinterval( -Im(z), Re(z) ) );
	return lx_cinterval( -Im(res), Re(res) );
}
//
//-- end acoth ----------------------------------------------------------------

// ---- sqrt1px2 --------------------------------------------------------------
//
lx_cinterval sqrt1px2(const lx_cinterval& z) throw()
// sqrt(1+z^2);
// Blomquist, 03.07.2008;
{
		const lx_real c = lx_real(1600,l_real(1));
		int stagsave = stagprec,
  stagmax = 30; // l_imath.cpp: sqrt() uses stagmax = 30;
  if (stagprec>stagmax) stagprec = stagmax;
	
  lx_cinterval res;
  lx_interval absz(abs(z));
  lx_real Inf_absz(Inf(absz));
	
  if (Inf_absz > c)
  {
	  absz = 1 / lx_interval(Inf_absz);
	  Inf_absz = Sup(absz);
	  res = lx_cinterval(lx_interval(-Inf_absz,Inf_absz),
								lx_interval(-Inf_absz,Inf_absz));
		 // res is the correcture interval, i.e.
		 // z + res or -z + res is the inclusionof sqrt{1+z^2}
	  res = (Inf(Re(z))>=0)? z + res : -z + res;
  }
  else 
  {
	  res =  lx_cinterval(lx_interval(0,l_interval(0)),
								 lx_interval(0,l_interval(1))); // res = i
	  if ( Sup(abs(z-res))<0.5 || Sup(abs(z+res))<0.5 )
	  {
		  res = lx_cinterval(-Im(z),Re(z)); // Res = i*z;
		    // (1 - i*z)*(1 + i*z) = 1+z^2;
		  res = sqrt( (1-res)*(1+res) );
	  }
	  else
		  res = sqrt(1+sqr(z));
  }
  if (sm_zero(Inf(Re(res))))
	  res = lx_cinterval(lx_interval(lx_real(0),Sup(Re(res))),Im(res));
  stagprec = stagsave;
  res = adjust(res);	
  return res;
}
// -- end sqrt1px2 ------------------------------------------------------------

lx_cinterval sqrt1mx2(const lx_cinterval& z) throw()
// sqrt(1-z^2);
// Blomquist, 03.07.2008;
{
		const lx_real c = lx_real(1600,l_real(1));
		int stagsave = stagprec,
  stagmax = 30; // l_imath.cpp: sqrt() uses stagmax = 30;
  if (stagprec>stagmax) stagprec = stagmax;
	
  lx_cinterval res,u;
  lx_interval absz(abs(z));
  lx_real Inf_absz(Inf(absz));
	
  if (Inf_absz > c)
  {
	  absz = 1 / lx_interval(Inf_absz);
	  Inf_absz = Sup(absz);
	  res = lx_cinterval(lx_interval(-Inf_absz,Inf_absz),
								lx_interval(-Inf_absz,Inf_absz)); // res = Delta
	  u = lx_cinterval(-Im(z),Re(z)); // u = i*z;
		// res is the correcture interval, i.e.
		// i*z + res or -i*z + res is the inclusion of sqrt{1-z^2}
	  res = (Inf(Im(z))>=0)? -u + res : u + res;
  }
  else 
  {
	  res = 1-z;  u = 1+z;
	  res = (Sup(abs(res))<0.5 || Sup(abs(u))<0.5)? sqrt(res*u) : sqrt(1-sqr(z));
  }
  if (sm_zero(Inf(Re(res))))
	  res = lx_cinterval(lx_interval(lx_real(0),Sup(Re(res))),Im(res));
  stagprec = stagsave;
  res = adjust(res);	
  return res;
}

lx_cinterval sqrtx2m1(const lx_cinterval& z) throw()
// sqrt(z^2-1);
// Blomquist, 03.07.2008;
{
		const lx_real c = lx_real(1600,l_real(1));
		int stagsave = stagprec,
  stagmax = 30; // l_imath.cpp: sqrt() uses stagmax = 30;
  if (stagprec>stagmax) stagprec = stagmax;
	
  lx_cinterval res,u;
  lx_interval absz(abs(z));
  lx_real Inf_absz(Inf(absz));
	
  if (Inf_absz > c)
  {
	  absz = 1 / lx_interval(Inf_absz);
	  Inf_absz = Sup(absz);
	  res = lx_cinterval(lx_interval(-Inf_absz,Inf_absz),
								lx_interval(-Inf_absz,Inf_absz)); // res = Delta
	   // res is the correcture interval, i.e.
	  res = (Inf(Re(z))>=0)? z + res : -z + res;
  }
  else 
  {
	  res = z-1;  u = z+1;
	  res = (Sup(abs(res))<0.5 || Sup(abs(u))<0.5)? sqrt(res*u) : sqrt(sqr(z)-1);
  }
	
  if (sm_zero(Inf(Re(res))))
	  res = lx_cinterval(lx_interval(lx_real(0),Sup(Re(res))),Im(res));
 
  stagprec = stagsave;
  res = adjust(res);	
  return res;
}

lx_cinterval sqrtp1m1(const lx_cinterval& z) throw()
// sqrt(1+z)-1;
// Blomquist, 08.07.2008;
{
		const real c = 0.125;
		int stagsave = stagprec,
  stagmax = 30; // l_imath.cpp: sqrt() uses stagmax = 30;
  if (stagprec>stagmax) stagprec = stagmax;
	
  lx_cinterval res;
  lx_interval absz(abs(z));
  lx_real Sup_absz(Sup(absz));
	
  if (Sup_absz < c)
	  res = z / (sqrt(z+1) + 1);
  else 
	  res = sqrt(z+1) - 1;
	
  stagprec = stagsave;
  res = adjust(res);	
  return res;
}

lx_cinterval expm1(const lx_cinterval& z) throw()
// exp(z) - 1;
// Blomquist, 09.08.2008;
{
		int stagsave = stagprec,
  stagmax = 30; // l_imath.cpp: sqrt() uses stagmax = 30;
  if (stagprec>stagmax) stagprec = stagmax;
	
  const lx_interval cancl_test = lx_interval(0,l_interval(0.995,1.005));
  lx_interval rez(Re(z)), imz(Im(z));
  lx_interval exp_x, sin_y, cos_y, h, Xt;
  lx_cinterval res;
	
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
  res = lx_cinterval(h,imz);
	
  stagprec = stagsave;
  res = adjust(res);
		
  return res;
}

lx_cinterval lnp1(const lx_cinterval& z) throw()
{// ln(1+z);
 // Blomquist, 11.08.2008;
	int stagsave = stagprec,
  	stagmax = 30; 
  	if (stagprec > stagmax) stagprec = stagmax;
  	const real c1 = 1e-128;
  	lx_cinterval y;
  	lx_interval abs_z(abs(z));
  	lx_real
		srez = Sup( Re(z) ),
		simz = Sup( Im(z) ),
		iimz = Inf( Im(z) );
	
	if (-1 <= z) //  z contains -1
		cxscthrow(STD_FKT_OUT_OF_DEF(
			"lx_cinterval lnp1(const lx_cinterval& z); z contains -1"));
	if ( srez<-1 && iimz<0 && simz>=0 ) 
		cxscthrow(STD_FKT_OUT_OF_DEF(
			"lx_cinterval lnp1(const lx_cinterval& z); z not allowed"));
	if (Sup(abs_z) < c1)
	{
		abs_z = Re(z);
		abs_z = lnp1( abs_z*(2+abs_z) + sqr(Im(z)) );
		times2pown(abs_z,-1);
		y = lx_cinterval( abs_z, arg(1+z) );
	}
	else
		y = Ln(1+z);
	stagprec = stagsave;
	y = adjust(y);

	return y;
}

//-- pow_all ------
//
//  Non-analytic power function for real l_interval exponent.
//
//  If 0 \not\in z, then compute four rectangular intervals that comprehend
//  an annulus which contains all values  zeta^phi, zeta in z, phi in p.
//
//  If 0 in z and negative reals in p, then abort execution
//  (potential modification: return the entire complex plane).
//
std::list<lx_cinterval> pow_all( const lx_cinterval& z, 
											const lx_interval& p ) throw()
{
	lx_interval abs_z = abs(z);

	if( 0.0 < Inf( abs_z ) )
	{
		lx_interval abs_z_p = exp( p * ln( abs_z ) );

      //  Inner and outer radii of the annulus are inf/sup( abs_z_n )
      //  Inscribed square has side length sqrt( 2 ) * rad_1
		lx_interval rad_1 = Sqrt2r_l_interval() * lx_interval(Inf( abs_z_p ));
		lx_interval rad_2 = lx_interval(Sup( abs_z_p ));

		std::list<lx_cinterval> res;

      //  4 intervals covering the annulus, counter-clockwise
		res.push_back(lx_cinterval(lx_interval( Inf( rad_1 ), Sup( rad_2 ) ),
									     lx_interval( -Sup( rad_1 ), Sup( rad_2 ) )));
		res.push_back(lx_cinterval(lx_interval( -Sup( rad_2 ), Sup( rad_1 ) ),
											lx_interval( Inf( rad_1 ), Sup( rad_2 ) ) ));
		res.push_back(lx_cinterval(lx_interval( -Sup( rad_2 ), -Inf( rad_1 ) ),
											lx_interval( -Sup( rad_2 ), Sup(rad_1) ) ) );
		res.push_back(lx_cinterval(lx_interval( -Sup( rad_1 ), Sup( rad_2 ) ),
											lx_interval( -Sup( rad_2 ), -Inf(rad_1) ) ) );

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
			lx_interval abs_z_p = exp( p * ln( lx_interval( Sup( abs_z ) ) ) );
			lx_real rad_p = Sup( abs_z_p );

			std::list<lx_cinterval> res;

			res.push_back( lx_cinterval( lx_interval( -rad_p, rad_p ),
												  lx_interval( -rad_p, rad_p ) ) );
			return res;
		}
		else
		{
		//
		//  The set   zeta^phi, zeta in z, phi in p   is unbounded
		//  if inf( p ) < 0.  0^p is undefined for p <= 0.
		//
			cxscthrow(STD_FKT_OUT_OF_DEF(
			"pow_all(lx_cinterval, lx_interval); 0^p is undefined for p <= 0."));
			std::list<lx_cinterval> res;
			return res;
		} 
	}
}
//
//-- end pow_all --------------------------------------------------------------

} // namespace cxsc
