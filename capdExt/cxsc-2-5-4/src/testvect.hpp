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

/* CVS $Id: testvect.hpp,v 1.24 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_TESTVECT_HPP_INCLUDED
#define _CXSC_TESTVECT_HPP_INCLUDED
// -------------------------- Vector - Tests

namespace cxsc {


template <class V,class W,class M,class S,int n>
class vectorconstr : public testclass
{
	public:
	vectorconstr( V &v, W &w,  M &m)
	{
		for(int i=1;i<=n;i++)
			w[i]=5*i+4;
		S s=S(1345.0);

		v=S(0.0);
		
		testing(nameof(v),nameof(v)+"("+nameof(s)+")");
		if(!tested())
		{
			if(V(s)[1]!=s)
				error();
			else
				ok();
		}
	
		testing(nameof(v),nameof(v)+"("+nameof(w)+")");
		if(!tested())
		{
			int err=0;
			for(int i=1;i<=n;i++)
			{
//				cout<<"V(w)[i]="<<V(w)[i]<<" ("<<5*i+4<<endl;
				if(V(w)[i]!=5*i+4)
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(v),nameof(v)+"("+nameof(w)+"())");
		if(!tested())
		{
			int err=0;
			for(int i=1;i<=n;i++)
				if(V(w(n))[i]!=5*i+4)
					err=1;
			if(err)
				error();
			else
				ok();
		}
		for(int i=1;i<=n;i++)
			for(int j=1;j<=n;j++)
				m[i][j]=50*i+3*j;
		testing(nameof(v),nameof(v)+"("+nameof(m)+")");
		if(!tested())
		{
			int err=0;
			for(int i=1;i<=n;i++)
			{
				if(V(M(m(1,1,1,n)))[i]!=50+i*3)
					err=1;
			}
			for(int i=1;i<=n;i++)
			{
				if(V(M(m(1,n,1,1)))[i]!=50*i+3)
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(v),nameof(v)+"("+nameof(m)+"())");
		if(!tested())
		{
			int err=0;
			for(int i=1;i<=n;i++)
			{
				if(V(m(1,1,1,n))[i]!=50+i*3)
					err=1;
			}
			for(int i=1;i<=n;i++)
			{
				if(V(m(1,n,1,1))[i]!=50*i+3)
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(v),nameof(v)+"("+nameof(m)+"[])");
		if(!tested())
		{
			int err=0;
			for(int i=1;i<=n;i++)
			{
				if(V(m[4])[i]!=200+3*i)
					err=1;
			}
			for(int i=1;i<=n;i++)
			{
				if(V(m[Col(4)])[i]!=50*i+12)
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
	}
};

template <class T,class V,class M,int n,class S>
class vectorassign : public testclass
{
   // Tests various assignments to vectors
   public:
      vectorassign(const T &v, const V &w, const M &m)
      {
         S s;
  //       scalarassign<S,int>    assign_from_int;
         
         
			for(int i=1;i<=n;i++)
			{
				for(int j=1;j<=n;j++)
					m[i][j]=50*i+3*j;
				w[i]=5*i+4;
			}

         testing(nameof(v),"operator []");
         if(!tested())
         {
            int i;
            for(i=1;i<=n;i++)
            {
               s=1234-i;
               v[i]=s;
            }
            int err=0;
            for(i=1;i<=n;i++)
            {
               s=1234-i;
               if(!(v[i]==s && !(v[i]!=s)))
                  err=1;
            }
            if(err)
               error();
            else
               ok();
         }
         testing(nameof(v),"operator ()(int).[]");
         if(!tested())
         {
            int i;
            for(i=1;i<=n;i++)
            {
               s=1234-i;
               v[i]=s;
            }
            int err=0;
            for(i=1;i<=n/2;i++)
            {
               s=1234-i;
               if(!(v(n/2)[i]==s && !(v(n/2)[i]!=s)))
                  err=1;
            }
            if(err)
               error();
            else
               ok();
         }
         testing(nameof(v),"operator ()(int,int).[]");
         if(!tested())
         {
            int i;
            for(i=1;i<=n;i++)
            {
               s=1234-i;
               v[i]=s;
            }
            int err=0;
            for(i=2;i<=n/2;i++)
            {
               s=1234-i;
               if(!(v(2,n/2)[i]==s && !(v(2,n/2)[i]!=s)))
                  err=1;
            }
            if(err)
               error();
            else
               ok();
         }
         testing(nameof(v),"operator =("+nameof(s)+")");
         if(!tested())
         {
        		s=2007;
				int i;
            v=s;
            int err=0;
            for(i=1;i<=n;i++)
            {
               if(!(v[i]==2007 && !(v[i]!=2007)))
                  err=1;
            }
            if(err)
               error();
            else
               ok();
         }
			testing(nameof(v),"operator =("+nameof(w)+")");
			if(!tested())
			{
				v=w;
				int err=0;
				for(int i=1;i<=n;i++)
					if(v[i]!=5*i+4)
						err=1;
				if(err)
					error();
				else
					ok();
				v=S(0);
			}
		testing(nameof(v),nameof(v)+"("+nameof(m)+")");
		if(!tested())
		{
			int err=0;
			v=M(m(1,1,1,n));
			for(int i=1;i<=n;i++)
			{
				if(v[i]!=50+3*i)
					err=1;
			}
			v=M(m(1,n,1,1));
			for(int i=1;i<=n;i++)
			{
				if(v[i]!=50*i+3)
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}

			testing(nameof(v),"operator =("+nameof(w)+"())");
			if(!tested())
			{
				v=w(50);
				int err=0;
				for(int i=1;i<=n;i++)
					if(v[i]!=5*i+4)
						err=1;
				if(err)
					error();
				else
					ok();
				v=S(0);
			}
		testing(nameof(v),nameof(v)+"("+nameof(m)+"())");
		if(!tested())
		{
			int err=0;
			v=m(1,1,1,n);
			for(int i=1;i<=n;i++)
			{
				if(v[i]!=50+3*i)
					err=1;
			}
			v=m(1,n,1,1);
			for(int i=1;i<=n;i++)
			{
				if(v[i]!=50*i+3)
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
			testing(nameof(v),"operator =("+nameof(m)+"[])");
			if(!tested())
			{
				v=m[3];
				int err=0;
				for(int i=1;i<=n;i++)
					if(v[i]!=150+3*i)
						err=1;
				if(err)
					error();
				else
					ok();
				v=S(0);
			}

      }
};

template <class A, class B, int n, class ERG>
class vectoraddsub : public testclass
{
   // Tests adding to vector of same type
   // (soon to different, but there is yet only one... ;)
   public:
      vectoraddsub(const A &a, const B &b)
      {
			ERG e(n);
         int i;
         
//         vectorassign<A,n,int> ai(a);
         
//         setfail(!ai);
         
         testing(nameof(e),"operator +("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
               a[i]= 123-i,b[i]=1023-2*i,e[i]=123+1023-3*i;

            int err=0;

            for(i=1;i<=n;i++)
               if(!((a+b)[i]==e[i] && !((a+b)[i]!=e[i])))
                  err=1;
            if(err)
               error();
            else
               ok();
         }
         
         testing(nameof(e),"operator -("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
               a[i]= 123-i,b[i]=1023-2*i,e[i]=123-1023+i;

            int err=0;

            for(i=1;i<=n;i++)
               if(!((a-b)[i]==e[i] && !((a-b)[i]!=e[i])))
                  err=1;
            if(err)
               error();
            else
               ok();
         }

         testing(nameof(e),"operator -("+nameof(a)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
               a[i]= 123-i,e[i]=-(123-i);

            int err=0;

            for(i=1;i<=n;i++)
               if( (-a)[i]!=e[i] || !((-a)[i]==e[i]))
                  err=1;
            if(err)
               error();
            else
               ok();
         }
         testing(nameof(e),"operator +("+nameof(a)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
               a[i]= 123-i,e[i]=+(123-i);

            int err=0;

            for(i=1;i<=n;i++)
               if( (+a)[i]!=e[i] || !((+a)[i]==e[i]))
                  err=1;
            if(err)
               error();
            else
               ok();
         }
		}
};

template <class A, class B, int n, class ERG>
class vectoraddsubassign : public testclass
{
   // Tests adding to vector of same type
   // (soon to different, but there is yet only one... ;)
   public:
      vectoraddsubassign(const A &a, const B &b)
      {
			ERG e(n);
         int i;
         
//         vectorassign<A,n,int> ai(a);
         
//         setfail(!ai);
         
         testing(nameof(e),"operator +=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            for(i=1;i<=n;i++)
               a[i]= 123-i,b[i]=1023-2*i,e[i]=123+1023-3*i;

            ERG c=(a+=b);

            int err=0;

            for(i=1;i<=n;i++)
            {
               //cout << c[i] << " " << a[i] << " " << e[i] << endl;
               if(c[i]!=e[i] || a[i]!=e[i])
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
               a[i]= 123-i,b[i]=1023-2*i,e[i]=123-1023+i;

            ERG c=(a-=b);

            int err=0;

            for(i=1;i<=n;i++)
            {
               if(c[i]!=e[i] || a[i]!=e[i])
                  err=1;
            }
            if(err)
               error();
            else
               ok();
         }

         
      }
};

template <class V1,class V2,int n,class E>
class vectorconv : public testclass
{
	public:
	vectorconv(const V1 &a, const V2 &b)
	{
		E e(n);
		testing(nameof(e),"operator |("+nameof(a)+","+nameof(b)+")");
		if(!tested())
		{
			for(int i=1;i<=n;i++)
				a[i]= 123-i,b[i]=3*i-20,e[i]=real(123-i)|real(3*i-20);

			int err=0;

			for(int i=1;i<=n;i++)
				if(!((a|b)[i]==e[i] && !((a|b)[i]!=e[i])))
					err=1;
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(e),"operator |("+nameof(b)+","+nameof(a)+")");
		if(!tested())
		{
			for(int i=1;i<=n;i++)
				a[i]= 123-i,b[i]=3*i-20,e[i]=real(123-i)|real(3*i-20);

			int err=0;

			for(int i=1;i<=n;i++)
				if(!((b|a)[i]==e[i] && !((b|a)[i]!=e[i])))
					err=1;
			if(err)
				error();
			else
				ok();
		}
	}
};

template <class V1,class V2,int n,class ERG>
class vectorconvassign : public testclass
{
	public:
	vectorconvassign(const V1 &a, const V2 &b)
	{
		ERG e(n);
		testing(nameof(e),"operator |=("+nameof(a)+","+nameof(b)+")");
		if(!tested())
		{
			for(int i=1;i<=n;i++)
				a[i]= 123-i,b[i]=3*i-20,e[i]=real(123-i)|real(3*i-20);

			ERG c=(a|=b);

			int err=0;

			for(int i=1;i<=n;i++)
			{
				if(c[i]!=e[i] || a[i]!=e[i])
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
	}
};

template <class V1,class V2,int n,class E,class S>
class vectorsect : public testclass
{
	public:
	vectorsect(const V1 &a, const V2 &b)
	{
		E e(n);
		testing(nameof(e),"operator &("+nameof(a)+","+nameof(b)+")");
		if(!tested())
		{
			for(int i=1;i<=n;i++)
			{
				a[i]=223-i,e[i]=b[i]=2*i;
			}
			SetInf(a,S(-100));
			int err=0;

			for(int i=1;i<=n;i++)
				if(!((a&b)[i]==e[i] && !((a&b)[i]!=e[i])))
					err=1;
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(e),"operator &("+nameof(b)+","+nameof(a)+")");
		if(!tested())
		{
			for(int i=1;i<=n;i++)
			{
				a[i]= 223-i,e[i]=b[i]=2*i;
			}
			SetInf(a,S(-100));

			int err=0;

			for(int i=1;i<=n;i++)
				if(!((b&a)[i]==e[i] && !((b&a)[i]!=e[i])))
					err=1;
			if(err)
				error();
			else
				ok();
		}
	}
};

template <class V1,class V2,int n,class ERG,class S>
class vectorsectassign : public testclass
{
	public:
	vectorsectassign(const V1 &a, const V2 &b)
	{
		ERG e(n);
		testing(nameof(e),"operator &=("+nameof(a)+","+nameof(b)+")");
		if(!tested())
		{
			for(int i=1;i<=n;i++)
			{
				a[i]= 223-i,e[i]=b[i]=2*i;
			}
			SetInf(a,S(-100));

			ERG c=(a&=b);

			int err=0;

			for(int i=1;i<=n;i++)
			{
				//cout << c[i] << " " << a[i] << " " << e[i] << endl;
				if(c[i]!=e[i] || a[i]!=e[i])
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
	}
};


template <class V1,class S,int n>
class vectorscalar : public testclass
{
	public:
	vectorscalar(const V1 &v)
	{
		S s;

		testing(nameof(v),"operator *("+nameof(v)+","+nameof(s)+")");
		if(!tested())
		{
			int err=0;
			s=S(2341);
			v=S(0.0);
			for(int i=1;i<=n;i++)
			{
				v[i]=3+4*i;
				if((v*s)[i]!=((3+4*i)*s))
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(v),"operator *("+nameof(s)+","+nameof(v)+")");
		if(!tested())
		{
			int err=0;
			s=S(2341);
			v=S(0.0);
			for(int i=1;i<=n;i++)
			{
				v[i]=3+4*i;
				if((s*v)[i]!=((3+4*i)*s))
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(v),"operator /("+nameof(v)+","+nameof(s)+")");
		if(!tested())
		{
			int err=0;
			s=S(2341);
			v=S(0.0);
			for(int i=1;i<=n;i++)
			{
				v[i]=3+4*i;
				if((v/s)[i]!=((3+4*i)/s))
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
	}
};

template <class V1,class S,class E,int n>
class vectorscalarassign : public testclass
{
	public:
	vectorscalarassign(const V1 &v, const E &e)
	{
		S s;

		testing(nameof(e),"operator *=("+nameof(v)+","+nameof(s)+")");
		if(!tested())
		{
			int err=0;
			s=S(2341);
			E t(n);
			for(int i=1;i<=n;i++)
			{
				e[i]=3+4*i;
				t[i]=(3+4*i)*s;
			}
			E r=(e*=s);
			for(int i=1;i<=n;i++)
			{
				if((r[i]!=t[i])||(e[i]!=t[i]))
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(e),"operator /=("+nameof(v)+","+nameof(s)+")");
		if(!tested())
		{
			int err=0;
			s=S(2341);
			E t(n);
			for(int i=1;i<=n;i++)
			{
				e[i]=3+4*i;
				t[i]=(3+4*i)/s;
			}
			E r=(e/=s);
			for(int i=1;i<=n;i++)
			{
				if((r[i]!=t[i])||(e[i]!=t[i]))
					err=1;
			}
			if(err)
				error();
			else
				ok();
		}
	}
};

template <class V1,class V2,int n,class S>
class vectormult : public testclass
{
	public:
	vectormult(const V1 &v, const V2 &w)
	{
		S s;

		testing(nameof(v),nameof(s)+" operator *("+nameof(v)+","+nameof(w)+")");
		if(!tested())
		{
			int err=0;
			s=S(0);
			for(int i=1;i<=n;i++)
			{
				v[i]=3+4*i;
				w[i]=5*i;
				s+=(3+4*i)*(5*i);
			}
			if(v*w!=s)
				err=1;
			if(err)
				error();
			else
				ok();
		}
		testing(nameof(v),nameof(s)+" operator *("+nameof(w)+","+nameof(v)+")");
		if(!tested())
		{
			int err=0;
			s=S(0);
			for(int i=1;i<=n;i++)
			{
				v[i]=3+4*i;
				w[i]=5*i;
				s+=(3+4*i)*(5*i);
			}
			if(w*v!=s)
				err=1;
			if(err)
				error();
			else
				ok();
		}
	}
};


template <class V1,class V2,int n,class DP>
class vectoraccu : public testclass
{
	public:
	vectoraccu(const V1 &v, const V2 &w)
	{
		DP dp,e;

		testing(nameof(v),"accumulate ("+nameof(dp)+","+nameof(v)+","+nameof(w)+")");
		if(!tested())
		{
			dp=0.0;
			for(int i=1;i<=n;i++)
			{
				v[i]=3+4*i;
				w[i]=5*i;
				accumulate(dp,(3+4*i),(5*i));
			}
			e=0.0;
			accumulate(e,v,w);
			if(e!=dp)
				error();
			else
				ok();
		}
		testing(nameof(v),"accumulate ("+nameof(dp)+","+nameof(w)+","+nameof(v)+")");
		if(!tested())
		{
			dp=0.0;
			for(int i=1;i<=n;i++)
			{
				v[i]=3+4*i;
				w[i]=5*i;
				accumulate(dp,5*i,3+4*i);
			}
			e=0.0;
			accumulate(e,w,v);
			if(e!=dp)
				error();
			else
				ok();
		}
	}
};

} // namespace cxsc 

#endif // _CXSC_TESTVECT_HPP_INCLUDED
