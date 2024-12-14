/******************************************************************************/
/*                                                                            */
/*                      sum up all partitional results                        */
/*                                                                            */
/******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "cxscconf.h"
#if WINDOWS_X86_32
#include <conio.h>
#endif

using namespace std;

int main ()
{
//  int i;
  int counter[4];
  counter[0] = 0;	                                            //equalities
  counter[1] = 0;	                                               //subsets
  counter[2] = 0;	                                             //supersets
  counter[3] = 0;	                                       //other relations

  int dummy;

  cout << endl;
  cout << endl;
  cout << endl;
  cout << "***********************************************************" << endl;
  cout << "*                      FINAL RESULT                       *" << endl;
  cout << "***********************************************************" << endl;

  ifstream sum("ergebnis.dat");

  if (!sum)
  {
    cout << "Error trying to read ergebnis.dat" << endl;
    return 1;
  }

  while (!sum.eof())
  {
    sum >> dummy;
    counter[0] = counter[0] + dummy;
    sum >> dummy;
    counter[1] = counter[1] + dummy;
    sum >> dummy;
    counter[2] = counter[2] + dummy;
    sum >> dummy;
    counter[3] = counter[3] + dummy;
  }

  sum.close();

  counter[0] = counter[0] - dummy;
  counter[1] = counter[1] - dummy;
  counter[2] = counter[2] - dummy;
  counter[3] = counter[3] - dummy;

  cout << endl;
  cout << "Number of equalities :" << counter[0] << endl;
  cout << "Number of subsets    :" << counter[1] << endl;
  cout << "Number of supersets  :" << counter[2] << endl;
  cout << "Other relations      :" << counter[3] << endl;

  if (counter[0] == 128)                  //128 = Number of intervals to compare
  {
    cout << endl;
    cout << "*********************************************************" << endl;
    cout << "*    All results are equal to the expected solutions.   *" << endl;
    cout << "*  It seems the libraries (C-XSC and CToolbox) are OK!  *" << endl;
    cout << "*********************************************************" << endl;
  }
  else if ((counter[0]+counter[2]) == 128) 
  {
    cout << endl;
    cout << "*********************************************************" << endl;
    cout << "*    All results are equal to or supersets of the       *" << endl;
    cout << "*                  expected solutions.                  *" << endl;
    cout << "*  It seems the libraries (C-XSC and CToolbox) are OK!  *" << endl;
    cout << "*********************************************************" << endl;
  }
  else
  {
     cout << endl;
     cout << "********************************************************" << endl;
     cout << "*     Warning: There are some unexpected solutions!    *" << endl;
     cout << "* It seems one library (C-XSC or CToolbox) in not OK!  *" << endl;
     cout << "*                                                      *" << endl;
     cout << "*                  !!!  CAUTION  !!!                   *" << endl;
     cout << "*                                                      *" << endl;
     cout << "* Do not use this library for interval arithmetic or   *" << endl;
     cout << "*        for numerics with result verification!        *" << endl;
     cout << "*                                                      *" << endl;
     cout << "*   Maybe you can try to compile the library without   *" << endl;
     cout << "*             any compiler optimization                *" << endl;
     cout << "*   (try 'make clean' and afterwards 'install_cxsc')   *" << endl;
     cout << "*                                                      *" << endl;
     cout << "*         If you have continuous problems:             *" << endl;
     cout << "*         e-mail: xsc@math.uni-wuppertal.de            *" << endl;    
     cout << "********************************************************" << endl;
     exit(1);
  }

#if WINDOWS_X86_32
	int taste=0;

	cout << endl;
	cout << "Press a key to close the window..." << endl;

	while(taste==0){
		taste = getch();
	}
#endif
  return 0;
}
