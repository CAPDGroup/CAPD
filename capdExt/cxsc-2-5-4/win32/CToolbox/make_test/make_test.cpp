/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.2.0)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2006 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**	      Author: Boris Hoeffgen <hoeffgen@uni-wuppertal.de>
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

#include <iostream>
#include <fstream>
#include "cxscconf.h"
#include <conio.h>
#include <direct.h>
#include <string.h>

using namespace std;

int main(int argc, char* argv[]){
  char dummy[1024];
  char buffer[_MAX_PATH];
  char dir[1024];

  _getcwd( buffer, _MAX_PATH );
  sprintf(dir,"%s\\examples\\test",buffer);
  _chdir(dir);
  ifstream sum("make_test.txt");

  if (!sum){
	_getcwd( buffer, _MAX_PATH );
	cout << buffer << endl;
    cout << "Error trying to read make_test.txt" << endl;
#if WINDOWS_X86_32
	int taste=0;

	cout << endl;
	cout << "Press a key to close the window..." << endl;

	while(taste==0){
		taste = getch();
	}
#endif
    return 1;
  }

  while (!sum.eof()){
	sum.getline(dummy,100);
	system(dummy);
  }

  return 0;
}
