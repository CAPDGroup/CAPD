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

/* CVS $Id: errtest.cpp,v 1.23 2014/01/30 17:23:45 cxsc Exp $ */

#include <iostream>
#include <memory>
#include "except.hpp"

namespace cxsc {

class err1701 : public ERROR_LRVECTOR, public WRONG_ROUNDING
{
	public:
	virtual int errnum() const throw() { return 1701; }
	err1701() throw() { }
};


using namespace std;

int main()
{
	try
	{
		for(;;) new char[10000];
		throw err1701();
		cout << "weiter"<<endl;
	}
	catch(const err1701 &e)
	{
		cout<<"Rundungsfehler: "<<e.errnum()<<endl;
	}
	catch(bad_alloc)
	{
		cout<<"Nicht genug Speicher!"<<endl;
	}
	catch(const ERROR_VECTOR &e)
	{
		cout <<"Vektor-Fehler: "<<e.errnum()<<endl;
	}
	catch(const ERROR_ALL &e)
	{
		cout <<"Fehler: "<<e.errnum()<<endl;
	}
	cout<<"ende"<<endl;
}

} // namespace cxsc

