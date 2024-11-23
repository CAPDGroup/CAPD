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

/* CVS $Id: testmatr.hpp,v 1.24 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_TESTMATR_HPP_INCLUDED
#define _CXSC_TESTMATR_HPP_INCLUDED
#include <string>
// -------------------------- Matrix - Tests

namespace cxsc {

template <class M,class V,class M2,class S,int n>
class matrixconstr : public testclass
{
	public:
	matrixconstr(const M &m, const V &v, const M2 &m2)
	{
		S s;
		for(int i=1;i<=n;i++)
			v[i]=(-5)*i;
		
		testing(nameof(m),nameof(m)+"("+nameof(s)+")");
		if(!tested())
		{
			s=3234;
			if((M(s)[1][1]!=s)||(Lb(M(s),ROW)!=1)||(Lb(M(s),COL)!=1)||(Ub(M(s),ROW)!=1)||(Ub(M(s),COL)!=1))
				error();
			else
				ok();
		}
		testing(nameof(m),nameof(m)+"("+nameof(v)+")");
		if(!tested())
		{
			int err=0;
			for(int i=1;i<=n;i++)
			{
				if((M(v)[i][1]!=(-5)*i)||(Lb(M(v),ROW)!=1)||(Lb(M(v),COL)!=1)||(Ub(M(v),ROW)!=n)||(Ub(M(v),COL)!=1))
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(m),nameof(m)+"("+nameof(v)+"())");
		if(!tested())
		{
			int err=0;
			for(int i=1;i<=n;i++)
			{
				if((M(v(n))[i][1]!=(-5)*i)||(Lb(M(v(n)),ROW)!=1)||(Lb(M(v(n)),COL)!=1)||(Ub(M(v(n)),ROW)!=n)||(Ub(M(v(n)),COL)!=1))
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(m),nameof(m)+"("+nameof(m2)+")");
		for(int i=1;i<=n;i++)
			for(int j=1;j<=n;j++)
				m2[i][j]=i+n*j;
		if(!tested())
		{
			int err=0;
			for(int i=1;i<=n&& !err;i++)
			{
				for(int j=1;j<=n&& !err;j++)
					if((M(m2)[i][j]!=i+n*j)||(Lb(M(m2),ROW)!=1)||(Lb(M(m2),COL)!=1)||(Ub(M(m2),ROW)!=n)||(Ub(M(m2),COL)!=n))
						err=1;
				cout<<"."<<flush;
			}
			cout<<endl;
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(m),nameof(m)+"("+nameof(m2)+"())");
		if(!tested())
		{
			int err=0;
			for(int i=1;i<=n&& !err;i++)
			{
				for(int j=1;j<=n&& !err;j++)
					if((M(m2(n,n))[i][j]!=i+n*j)||(Lb(M(m2(n,n)),ROW)!=1)||(Lb(M(m2(n,n)),COL)!=1)||(Ub(M(m2(n,n)),ROW)!=n)||(Ub(M(m2(n,n)),COL)!=n))
						err=1;
				cout<<"."<<flush;
			}
			cout<<endl;
			if(err)
				error();
			else
				ok();
		}
	}
};


template <class M,class V,class M2,class S,int n>
class matrixassign : public testclass
{
   // Tests various assignments to matrices
   public:
      matrixassign(const M &x, const V &v, const M2 &y)
      {
         S s;
         
         testing(nameof(x),"operator [][].operator =("+nameof(s)+")");
         if(!tested())
         {
            int i,j;
            for(i=1;i<=n;i++)
            {
					for(j=1;j<=n;j++)
					{
               	s=1234*j-i;
               	x[i][j]=s;
					}
            }
            int err=0;
            for(i=1;i<=n;i++)
            {
					for(j=1;j<=n;j++)
					{
               	s=1234*j-i;
               	if(!(x[i][j]==s && !(x[i][j]!=s)))
               	   err=1;
					}
            }
            if(err)
               error();
            else
               ok();
         }
         testing(nameof(x),"operator [Col()][].operator =("+nameof(s)+")");
         if(!tested())
         {
            int i,j;
            for(i=1;i<=n;i++)
            {
					for(j=1;j<=n;j++)
					{
               	s=1234*j-i;
               	x[Col(j)][i]=s;
					}
            }
            int err=0;
            for(i=1;i<=n;i++)
            {
					for(j=1;j<=n;j++)
					{
               	s=1234*j-i;
               	if(!(x[i][j]==s && !(x[i][j]!=s)))
               	   err=1;
					}
            }
            if(err)
               error();
            else
               ok();
         }

         testing(nameof(x),"operator =("+nameof(y)+")");
         if(!tested())
         {
            int i,j;
				x=S(0.0);
            for(i=1;i<=n;i++)
            {
					for(j=1;j<=n;j++)
					{
						s=1234*j-i;
               	y[i][j]=s;
					}
            }
            int err=0;
				x=y;
            for(i=1;i<=n;i++)
            {
					for(j=1;j<=n;j++)
					{
               	s=1234*j-i;
               	if(x[i][j]!=s || !(x[i][j]==s))
               	   err=1;
					}
            }
            if(err)
               error();
            else
               ok();
         }

         testing(nameof(x),"operator (int,int).[][]=("+nameof(s)+")");
         if(!tested())
         {
            int i,j;
				for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
						x[i][j]=1234*j+i;
         
            int err=0;
            for(i=1;i<=n-10;i++)
            {
					for(j=1;j<=n-10;j++)
					{
               	s=1234*j+i;
               	if(!(x(n-10,n-10)[i][j]==s && !(x(n-10,n-10)[i][j]!=s)))
               	   err=1;
					}
            }
            if(err)
               error();
            else
               ok();
         }
         testing(nameof(x),"operator (int,int,int,int).[][]=("+nameof(s)+")");
         if(!tested())
         {
            int i,j;
				for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
						x[i][j]=1234*j+i;
         
            int err=0;
            for(i=10;i<=n-10;i++)
            {
					for(j=10;j<=n-10;j++)
					{
               	s=1234*j+i;
               	if(!(x(10,n-10,10,n-10)[i][j]==s && !(x(10,n-10,10,n-10)[i][j]!=s)))
               	   err=1;
					}
            }
            if(err)
               error();
            else
               ok();
         }
         testing(nameof(x),"operator (int,int,int,int).[][]=("+nameof(s)+")");
         if(!tested())
         {
            int i,j;
				for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
						x[i][j]=1234*j+i;
         
            int err=0;
            for(i=10;i<=n-10;i++)
            {
					for(j=10;j<=n-10;j++)
					{
						s=123*j-i;
               	x(10,n-10,10,n-10)[i][j]=s;
               	if(!(x[i][j]==s && !(x[i][j]!=s)))
               	   err=1;
					}
            }
            if(err)
               error();
            else
               ok();
         }
         testing(nameof(x),"operator =("+nameof(y)+"())");
         if(!tested())
         {
				x=S(0.0);
            int i,j;
				for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
						y[i][j]=1234*j+i;
         
            int err=0;
				x=y(1,n,1,n);
            for(i=1;i<=n;i++)
            {
					for(j=1;j<=n;j++)
					{
               	s=1234*j+i;
               	if(!(x[i][j]==s && !(x[i][j]!=s)))
               	   err=1;
					}
            }
            if(err)
               error();
            else
               ok();
         }
         testing(nameof(x),"operator ().=("+nameof(y)+")");
         if(!tested())
         {
            int i,j;
				for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
						x[i][j]=1234*j+i;
         
            for(i=10;i<=n-10;i++)
            {
					for(j=10;j<=n-10;j++)
					{
               	s=123*j-i;
               	y[i][j]=s;
					}
            }
            int err=0;
				x(1,n,1,n)=y;
            for(i=10;i<=n-10;i++)
            {
					for(j=10;j<=n-10;j++)
					{
               	s=123*j-i;
               	if(!(x[i][j]==s && !(x[i][j]!=s)))
               	   err=1;
					}
            }
            if(err)
               error();
            else
               ok();
         }
         testing(nameof(x),"operator ().=("+nameof(y)+"())");
         if(!tested())
         {
            int i,j;
				for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
						x[i][j]=1234*j+i;
         
            for(i=1;i<=n;i++)
            {
					for(j=1;j<=n;j++)
					{
               	s=123*j-i;
               	y[i][j]=s;
					}
            }
            int err=0;
				x(10,n-10,10,n-10)=y(10,n-10,10,n-10);
            for(i=10;i<=n-10;i++)
            {
					for(j=10;j<=n-10;j++)
					{
               	s=123*j-i;
               	if(!(x[i][j]==s && !(x[i][j]!=s)))
               	   err=1;
					}
            }
            if(err)
               error();
            else
               ok();
         }
			testing(nameof(x),"operator =("+nameof(v)+")");
			if(!tested())
			{
				int err=0;
				for(int i=1;i<=n;i++)
					v[i]=3*i+2;
				if(nameof(x).substr(nameof(x).length()-5,5)=="slice")
					x(n,1)=v;
				else
					x=v;
				for(int i=1;i<=n;i++)
				{
					if(x[i][1]!=3*i+2)
						err=1;
					if((nameof(x).substr(nameof(x).length()-5,5)!="slice")&&((Lb(x,ROW)!=1)||(Lb(x,COL)!=1)||(Ub(x,ROW)!=n)||(Ub(x,COL)!=1)))
						err=1;
				}
				if(err)
					error();
				else
					ok();
				x=y;
			}
			testing(nameof(x),"operator =("+nameof(v)+"())");
			if(!tested())
			{
				x=S(0.0);
				int err=0;
				for(int i=1;i<=n;i++)
					v[i]=3*i+2;
				if(nameof(x).substr(nameof(x).length()-5,5)=="slice")
					x(n,1)=v;
				else
					x=v(n);
//				cout << "v"<<endl<<v<<endl;
//				cout<<"x"<<endl<<x<<endl;
				for(int i=1;i<=n;i++)
				{
//					cout<<"i "<<i<<endl;
					if(x[i][1]!=3*i+2)
						err=1;
					if((nameof(x).substr(nameof(x).length()-5,5)!="slice")&&((Lb(x,ROW)!=1)||(Lb(x,COL)!=1)||(Ub(x,ROW)!=n)||(Ub(x,COL)!=1)))
						err=1;
				}
				if(err)
					error();
				else
					ok();
				x=y;
			}
      }
};

template <class A,class B,class E,class DP,int n>
class matrixmult: public testclass
{
	public:
		matrixmult(const A &a, const B &b)
		{
			E e(n,n);
			DP dp;
//         int i,j,err;
         
//         matrixassign<A,n,m,int> ai;
         
//         setfail(!ai);
         
      	testing(nameof(e),"operator *("+nameof(a)+","+nameof(b)+")");
			if(!tested())
			{
				int err=0;
				E e(n,n);
				for(int i=1;i<=n;i++)
				{
					for(int j=1;j<=n;j++)
					{
						a[i][j]=10+2*i-j;
						b[i][j]=100-3*i+4*j;
					}
				}
				e=a*b;
				for(int i=1;i<=n;i++)
				{
					for(int j=1;j<=n;j++)
					{
						dp=0.0;
						for(int k=1;k<=n;k++)
							accumulate(dp,10+2*i-k,100-3*k+4*j);
						if(e[i][j]!=rnd(dp))
							err=1;
					}
				}
				if(err)
					error();
				else
					ok();		
			}
      	testing(nameof(e),"operator *("+nameof(b)+","+nameof(a)+")");
			if(!tested())
			{
				int err=0;
				E e(n,n);
				for(int i=1;i<=n;i++)
				{
					for(int j=1;j<=n;j++)
					{
						a[i][j]=10+2*i-j;
						b[i][j]=100-3*i+4*j;
					}
				}
				e=b*a;
				for(int i=1;i<=n;i++)
				{
					for(int j=1;j<=n;j++)
					{
						dp=0.0;
						for(int k=1;k<=n;k++)
							accumulate(dp,10+2*k-j,100-3*i+4*k);
						if(e[i][j]!=rnd(dp))
							err=1;
					}
				}
				if(err)
					error();
				else
					ok();		
			}
      	testing(nameof(e),"operator *=("+nameof(a)+","+nameof(b)+")");
			if(!tested())
			{
				int err=0;
				E e(n,n);
				for(int i=1;i<=n;i++)
				{
					for(int j=1;j<=n;j++)
					{
						a[i][j]=10+2*i-j;
						b[i][j]=100-3*i+4*j;
					}
				}
				e=(a*=b);
				for(int i=1;i<=n;i++)
				{
					for(int j=1;j<=n;j++)
					{
						dp=0.0;
						for(int k=1;k<=n;k++)
							accumulate(dp,10+2*i-k,100-3*k+4*j);
						if((e[i][j]!=rnd(dp))||(a[i][j]!=rnd(dp)))
							err=1;
					}
				}
				if(err)
					error();
				else
					ok();		
			}
		}
};

template <class A,class S,class E,int n>
class matrixscalarmult : public testclass
{
	public:
		matrixscalarmult(const A &a)
		{
			E e(n,n);
			S s;
			int i,j,err;
//			matrixassign<A,n,n,S> ai;
//			setfail(!ai);

			testing(nameof(a),"operator *("+nameof(s)+")");
			if(!tested())
			{
				err=0;
				s=21;
				for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
					{
						a[i][j]=123*i+4*j;
						e[i][j]=(123*i+4*j)*s;
					}
				}
				E c=a*s;
				for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
   	            if(!(c[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
						{
	                  err=1;
						}
            if(err)
               error();
            else
               ok();
         }

			testing(nameof(a),"operator /("+nameof(s)+")");
			if(!tested())
			{
				err=0;
				s=21;
				for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
					{
						a[i][j]=123*i+4*j,e[i][j]=(123*i+4*j)/s;
					}
				}
				E c=a/s;
				for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
   	            if(!(c[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
	                  err=1;
				}
            if(err)
               error();
            else
               ok();
         }
			testing(nameof(a),"operator *=("+nameof(s)+")");
			if(!tested())
			{
				err=0;
				s=21;
				for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
					{
						a[i][j]=123*i+4*j;
						e[i][j]=(123*i+4*j)*s;
					}
				}
				a*=s;
				for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
   	            if(!(a[i][j]==e[i][j] && !(a[i][j]!=e[i][j])))
						{
	                  err=1;
						}
            if(err)
               error();
            else
               ok();
         }

			testing(nameof(a),"operator /("+nameof(s)+")");
			if(!tested())
			{
				err=0;
				s=21;
				for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
					{
						a[i][j]=123*i+4*j,e[i][j]=(123*i+4*j)/s;
					}
				}
				a/=s;
				for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
   	            if(!(a[i][j]==e[i][j] && !(a[i][j]!=e[i][j])))
	                  err=1;
				}
            if(err)
               error();
            else
               ok();
         }
		}
};
			

template <class A,class B,class E,int n>
class matrixaddsub : public testclass
{
   // Tests adding to matrix of same type
   // (soon to diffrent, but there is yet only one... ;)
   public:
      matrixaddsub(const A &a, const B &b)
      {
         int i,j;
         E c,e(n,n);
         testing(nameof(e),"operator +("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
	               a[i][j]= 123-i+j,b[i][j]=1023-2*i-3*j,e[i][j]=123+1023-3*i-2*j;

            int err=0;
				c=a+b;

            for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
   	            if(!(c[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
	                  err=1;
				}
            if(err)
               error();
            else
               ok();
         }
         
         testing(nameof(e),"operator +=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
	               a[i][j]= 123-i+j,b[i][j]=1023-2*i-3*j,e[i][j]=123+1023-3*i-2*j;


            c=(a+=b);

            int err=0;

            for(i=1;i<=n;i++)
            {
               //cout << c[i][j] << " " << a[i][j] << " " << e[i][j] << endl;
					for(j=1;j<=n;j++)
	               if(c[i][j]!=e[i][j] || a[i][j]!=e[i][j])
  	                err=1;
            }
            if(err)
               error();
            else
               ok();
         }
         
         testing(nameof(e),"operator -("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
	               a[i][j]= 123-i+j,b[i][j]=1023-2*i-3*j,e[i][j]=123-1023+i+4*j;

            int err=0;
				c=(a-b);

            for(i=1;i<=n;i++)
            {
					for(j=1;j<=n;j++)
	               if(c[i][j]!=e[i][j] || !(c[i][j]==e[i][j]))
  	                err=1;
            }

            if(err)
               error();
            else
               ok();
         }

         testing(nameof(e),"operator -=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
	               a[i][j]= 123-i+j,b[i][j]=1023-2*i-3*j,e[i][j]=123-1023+i+4*j;

            c=(a-=b);

            int err=0;

            for(i=1;i<=n;i++)
            {
					for(j=1;j<=n;j++)
						if(c[i][j]!=e[i][j] || a[i][j]!=e[i][j])
  	                err=1;
            }

            if(err)
               error();
            else
               ok();
         }

         testing(nameof(e),"operator -("+nameof(a)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
	               a[i][j]= 123-i+j,e[i][j]=-(123-i+j);

            int err=0;

				b= -a;
            for(i=1;i<=n;i++)
            {
					for(j=1;j<=n;j++)
						if( b[i][j]!=e[i][j] || !(b[i][j]==e[i][j]))
  	                err=1;
            }

            if(err)
               error();
            else
               ok();
         }
         testing(nameof(e),"operator +("+nameof(a)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
	               a[i][j]= 123-i+j,e[i][j]=+(123-i+j);

            int err=0;
				b= +a;

            for(i=1;i<=n;i++)
            {
					for(j=1;j<=n;j++)
						if( b[i][j]!=e[i][j] || !(b[i][j]==e[i][j]))
  	                err=1;
            }
            if(err)
               error();
            else
               ok();
         }
         
         testing(nameof(e),"operator +("+nameof(a)+"(),"+nameof(b)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
	          	{
						a[i][j]= 123-i+j;
						if((i>=10)&&(j>=10)&&(i<=n-10)&&(j<=n-10))
						{
							b[i][j]=e[i][j]=1023-2*i-3*j;
							e[i][j]+=123-i+j;
						}
					}
				}
            int err=0;
				c=a(10,n-10,10,n-10)+B(b(10,n-10,10,n-10));

            for(i=10;i<=n-10;i++)
				{
					for(j=10;j<=n-10;j++)
   	            if(!(c[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
	                  err=1;
				}
            if(err)
               error();
            else
               ok();
         }
         
         testing(nameof(e),"operator +("+nameof(a)+","+nameof(b)+"() )");
         if(!tested())
         {
            for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
	          	{
						a[i][j]= 123-i+j;
						if((i>=10)&&(j>=10)&&(i<=n-10)&&(j<=n-10))
						{
							b[i][j]=e[i][j]=1023-2*i-3*j;
							e[i][j]+=123-i+j;
						}
					}
				}
            int err=0;
				c=B(b(10,n-10,10,n-10))+a(10,n-10,10,n-10);

            for(i=10;i<=n-10;i++)
				{
					for(j=10;j<=n-10;j++)
   	            if(!(c[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
	                  err=1;
				}
            if(err)
               error();
            else
               ok();
         }
         
         testing(nameof(e),"operator +("+nameof(a)+"() ,"+nameof(b)+"() )");
         if(!tested())
         {
            for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
	          	{
						a[i][j]= 123-i+j;
						b[i][j]=1023-2*i-3*j;
						if((i>=10)&&(j>=10)&&(i<=n-10)&&(j<=n-10))
						{
							e[i][j]=1023-2*i-3*j+123-i+j;
						}
					}
				}
            int err=0;
				c=b(10,n-10,10,n-10)+a(10,n-10,10,n-10);

            for(i=10;i<=n-10;i++)
				{
					for(j=10;j<=n-10;j++)
   	            if(!(c[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
	                  err=1;
				}
            if(err)
               error();
            else
               ok();
         }
         
         testing(nameof(e),"operator ().+=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
	          	{
						a[i][j]=e[i][j]=123-i+j;
						if((i>=10)&&(j>=10)&&(i<=n-10)&&(j<=n-10))
						{
							e[i][j]+=b[i][j]=1023-2*i-3*j;
						}
					}
				}
            int err=0;
				c=(a(10,n-10,10,n-10)+=B(b(10,n-10,10,n-10)));

            for(i=10;i<=n-10;i++)
				{
					for(j=10;j<=n-10;j++)
   	            if(!(a[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
	                  err=1;
				}
            if(err)
               error();
            else
               ok();
         }
			
         testing(nameof(e),"operator ().+=("+nameof(a)+","+nameof(b)+"() )");
         if(!tested())
         {
            for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
	          	{
						a[i][j]=e[i][j]=123-i+j;
						b[i][j]=1023-2*i-3*j;
						if((i>=10)&&(j>=10)&&(i<=n-10)&&(j<=n-10))
						{
							e[i][j]+=b[i][j];
						}
					}
				}
            int err=0;
				c=(a(10,n-10,10,n-10)+=b(10,n-10,10,n-10));

            for(i=10;i<=n-10;i++)
				{
					for(j=10;j<=n-10;j++)
   	            if(!(a[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
	                  err=1;
				}
            if(err)
               error();
            else
               ok();
         }
			
		testing(nameof(e),"operator +=("+nameof(a)+","+nameof(b)+"() )");
		if(!tested())
		{
			for(i=1;i<=n;i++)
			{
				for(j=1;j<=n;j++)
				{
					e[i][j]= 123-i+j+1023-2*i-3*j;
					b[i][j]=1023-2*i-3*j;
					a[i][j]=123-i+j;
				}
			}
			int err=0;
			c=(a+=b(n,n));

			for(i=1;i<=n;i++)
			{
				for(j=1;j<=n;j++)
					if(!(a[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
						err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		
         testing(nameof(e),"operator -("+nameof(a)+"(),"+nameof(b)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
	          	{
						a[i][j]= 123-i+j;
						if((i>=10)&&(j>=10)&&(i<=n-10)&&(j<=n-10))
						{
							e[i][j]=123-i+j;
							e[i][j]-=b[i][j]=1023-2*i-3*j;
						}
					}
				}
            int err=0;
				c=a(10,n-10,10,n-10)-B(b(10,n-10,10,n-10));

            for(i=10;i<=n-10;i++)
				{
					for(j=10;j<=n-10;j++)
   	            if(!(c[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
	                  err=1;
				}
            if(err)
               error();
            else
               ok();
         }
         
         testing(nameof(e),"operator -("+nameof(a)+","+nameof(b)+"() )");
         if(!tested())
         {
            for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
	          	{
						a[i][j]= 123-i+j;
						if((i>=10)&&(j>=10)&&(i<=n-10)&&(j<=n-10))
						{
							e[i][j]=b[i][j]=1023-2*i-3*j;
							e[i][j]-=123-i+j;
						}
					}
				}
            int err=0;
				c=B(b(10,n-10,10,n-10))-a(10,n-10,10,n-10);

            for(i=10;i<=n-10;i++)
				{
					for(j=10;j<=n-10;j++)
   	            if(!(c[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
	                  err=1;
				}
            if(err)
               error();
            else
               ok();
         }
         
         testing(nameof(e),"operator -("+nameof(a)+"() ,"+nameof(b)+"() )");
         if(!tested())
         {
            for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
	          	{
						a[i][j]= 123-i+j;
						b[i][j]=1023-2*i-3*j;
						if((i>=10)&&(j>=10)&&(i<=n-10)&&(j<=n-10))
						{
							e[i][j]=1023-2*i-3*j-(123-i+j);
						}
					}
				}
            int err=0;
				c=b(10,n-10,10,n-10)-a(10,n-10,10,n-10);

            for(i=10;i<=n-10;i++)
				{
					for(j=10;j<=n-10;j++)
   	            if(!(c[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
	                  err=1;
				}
            if(err)
               error();
            else
               ok();
         }
         
		testing(nameof(e),"operator ().-=("+nameof(a)+","+nameof(b)+")");
		if(!tested())
		{
			for(i=1;i<=n;i++)
			{
				for(j=1;j<=n;j++)
				{
					a[i][j]=e[i][j]= 123-i+j;
					if((i>=10)&&(j>=10)&&(i<=n-10)&&(j<=n-10))
					{
						e[i][j]-=b[i][j]=1023-2*i-3*j;
					}
				}
			}
			int err=0;
			c=(a(10,n-10,10,n-10)-=B(b(10,n-10,10,n-10)));

			for(i=10;i<=n-10;i++)
			{
				for(j=10;j<=n-10;j++)
					if(!(a[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
						err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		
		testing(nameof(e),"operator ().-=("+nameof(a)+","+nameof(b)+"() )");
		if(!tested())
		{
			for(i=1;i<=n;i++)
			{
				for(j=1;j<=n;j++)
				{
					a[i][j]=e[i][j]= 123-i+j;
					b[i][j]=1023-2*i-3*j;
					if((i>=10)&&(j>=10)&&(i<=n-10)&&(j<=n-10))
					{
						e[i][j]-=b[i][j];
					}
				}
			}
			int err=0;
			c=(a(10,n-10,10,n-10)-=b(10,n-10,10,n-10));

			for(i=10;i<=n-10;i++)
			{
				for(j=10;j<=n-10;j++)
					if(!(c[i][j]==e[i][j] && !(a[i][j]!=e[i][j])))
						err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		
		testing(nameof(e),"operator -=("+nameof(a)+","+nameof(b)+"() )");
		if(!tested())
		{
			for(i=1;i<=n;i++)
			{
				for(j=1;j<=n;j++)
				{
					e[i][j]= 123-i+j-(1023-2*i-3*j);
					b[i][j]=1023-2*i-3*j;
					a[i][j]=123-i+j;
				}
			}
			int err=0;
			c=(a-=b(1,n,1,n));

			for(i=1;i<=n;i++)
			{
				for(j=1;j<=n;j++)
					if(!(c[i][j]==e[i][j] && !(a[i][j]!=e[i][j])))
						err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		
         testing(nameof(e),"operator -("+nameof(a)+"() )");
         if(!tested())
         {
            for(i=1;i<=n;i++)
				{
					for(j=1;j<=n;j++)
	          	{
						a[i][j]= 123-i+j;
					}
				}
            int err=0;

				c= -a(10,n-10,10,n-10);
            for(i=10;i<=n-10;i++)
            {
					for(j=10;j<=n-10;j++)
						if( c[i][j]!=-a[i][j] || !(c[i][j]==-a[i][j]))
  	                err=1;
            }

            if(err)
               error();
            else
               ok();
         }
         testing(nameof(e),"operator +("+nameof(a)+"() )");
         if(!tested())
         {
            for(i=1;i<=n;i++)
					for(j=1;j<=n;j++)
	               a[i][j]= 123-i+j;

            int err=0;
				c= +a;

            for(i=10;i<=n-10;i++)
            {
					for(j=10;j<=n-10;j++)
						if( c[i][j]!=a[i][j] || !(c[i][j]==a[i][j]))
  	                err=1;
            }
            if(err)
               error();
            else
               ok();
         }
         
      }
};

template <class M1,class M2,int n,class E>
class matrixconv: public testclass
{
	public:
	matrixconv(const M1 &a, const M2 &b)
	{
		int i,j;
		E c,e(n,n);
		testing(nameof(e),"operator |("+nameof(a)+","+nameof(b)+")");
		if(!tested())
		{
			for(i=1;i<=n;i++)
				for(j=1;j<=n;j++)
					a[i][j]= 23-i+j,b[i][j]=1023-2*i-3*j,e[i][j]=(23-i+j)|(1023-2*i-3*j);

			int err=0;
			c=a|b;

			for(i=1;i<=n;i++)
			{
				for(j=1;j<=n;j++)
					if(!(c[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
						err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(e),"operator |("+nameof(b)+","+nameof(a)+")");
		if(!tested())
		{
			for(i=1;i<=n;i++)
				for(j=1;j<=n;j++)
					a[i][j]= 23-i+j,b[i][j]=1023-2*i-3*j,e[i][j]=(23-i+j)|(1023-2*i-3*j);

			int err=0;
			c=b|a;

			for(i=1;i<=n;i++)
			{
				for(j=1;j<=n;j++)
					if(!(c[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
						err=1;
			}
			if(err)
				error();
			else
				ok();
		}
	}
};

template <class M1,class M2,int n,class E>
class matrixconvassign: public testclass
{
	public:
	matrixconvassign(const M1 &a, const M2 &b)
	{
		int i,j;
		E c,e(n,n);
		testing(nameof(e),"operator |=("+nameof(a)+","+nameof(b)+")");
		if(!tested())
		{
			for(i=1;i<=n;i++)
				for(j=1;j<=n;j++)
					a[i][j]= 23-i+j,b[i][j]=1023-2*i-3*j,e[i][j]=(123-i+j)|(1023-2*i-3*j);


			c=(a|=b);

			int err=0;

			for(i=1;i<=n;i++)
			{
				//cout << c[i][j] << " " << a[i][j] << " " << e[i][j] << endl;
				for(j=1;j<=n;j++)
					if(c[i][j]!=e[i][j] || a[i][j]!=e[i][j])
					 err=1;
			}
			if(err)
				error();
			else
				ok();
		}
	}
};
template <class M1,class M2,int n,class E>
class matrixsect: public testclass
{
	public:
	matrixsect(const M1 &a, const M2 &b)
	{
		int i,j;
		E c,e(n,n);
		testing(nameof(e),"operator &("+nameof(a)+","+nameof(b)+")");
		if(!tested())
		{
			for(i=1;i<=n;i++)
				for(j=1;j<=n;j++)
					a[i][j]= 23-i+j,b[i][j]=1023-2*i-3*j,e[i][j]=(23-i+j)&(1023-2*i-3*j);

			int err=0;
			c=a&b;

			for(i=1;i<=n;i++)
			{
				for(j=1;j<=n;j++)
					if(!(c[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
						err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(e),"operator &("+nameof(b)+","+nameof(a)+")");
		if(!tested())
		{
			for(i=1;i<=n;i++)
				for(j=1;j<=n;j++)
					a[i][j]= 23-i+j,b[i][j]=1023-2*i-3*j,e[i][j]=(23-i+j)&(1023-2*i-3*j);

			int err=0;
			c=b&a;

			for(i=1;i<=n;i++)
			{
				for(j=1;j<=n;j++)
					if(!(c[i][j]==e[i][j] && !(c[i][j]!=e[i][j])))
						err=1;
			}
			if(err)
				error();
			else
				ok();
		}
	}
};

template <class M1,class M2,int n,class E>
class matrixsectassign: public testclass
{
	public:
	matrixsectassign(const M1 &a, const M2 &b)
	{
		int i,j;
		E c,e(n,n);
		testing(nameof(e),"operator &=("+nameof(a)+","+nameof(b)+")");
		if(!tested())
		{
			for(i=1;i<=n;i++)
				for(j=1;j<=n;j++)
					a[i][j]= 23-i+j,b[i][j]=1023-2*i-3*j,e[i][j]=(123-i+j)&(1023-2*i-3*j);


			c=(a&=b);

			int err=0;

			for(i=1;i<=n;i++)
			{
				//cout << c[i][j] << " " << a[i][j] << " " << e[i][j] << endl;
				for(j=1;j<=n;j++)
					if(c[i][j]!=e[i][j] || a[i][j]!=e[i][j])
					 err=1;
			}
			if(err)
				error();
			else
				ok();
		}
	}
};
     
} // namespace cxsc 
   
#endif // _CXSC_TESTMATR_HPP_INCLUDED
