/******************************************************************************/
/*                                                                            */
/*                      compares two interval-data-files                      */
/*                                                                            */
/******************************************************************************/

#include <iostream>
#include <fstream>
#include <cstring>
#include <interval.hpp>                                      //interval-data-typ

using namespace std;
using namespace cxsc;

template<typename T> static T basename(T s)
{
	T p;
	size_t z = strlen(s);
	for (p = s + z - 1; p >= s && *p != '/'; --p)
		;
	return (p >= s) ? p : s;
}

int main (int argAnz, char* argFeld[])
{
  int i;
  int counter[4];
  counter[0] = 0;	                                            //equalities
  counter[1] = 0;	                                               //subsets
  counter[2] = 0;	                                             //supersets
  counter[3] = 0;	                                        //other relation

  interval data[100];	                                    //max number of data
  interval vgl[100];

  int anzdaten;
  int anzvgl;

  cout << Scientific << SetPrecision(23,15);

  cout << endl;
  cout << endl;
  cout << endl;
  cout << "***********************************************************" << endl;
  cout << " checking program " << argFeld[1] << endl;
  cout << "***********************************************************" << endl;

  ifstream eindat1(argFeld[1]);
  if (!eindat1)
  {
    cout << "Error trying to read new data" << endl;
    return 1;
  }

/******************** special cases (begin) ********************/
  if (strcmp(basename(argFeld[1]), "hess_ex.out") == 0)
  {
    while (eindat1.peek()!='I' && !eindat1.eof())
    {
      eindat1.get();
      if (eindat1.eof()) break;
    }
  }

  if (strcmp(basename(argFeld[1]), "lop_ex.out1") == 0)
  {
    while (eindat1.peek()!='I' && !eindat1.eof())
    {
      eindat1.get();
      if (eindat1.eof()) break;
    }
  }

  if (strcmp(basename(argFeld[1]), "lop_ex.out2") == 0)
  {
    while (eindat1.peek()!='I' && !eindat1.eof())
    {
      eindat1.get();
      if (eindat1.eof()) break;
    }
  }

/********************* special cases (end) *********************/

  anzdaten = 0;
  while (!eindat1.eof())
  {
	while (eindat1.peek()!='[' && !eindat1.eof()) eindat1.get();
    if (eindat1.eof()) break;

    eindat1 >> data[anzdaten];
    anzdaten++;
  }

  eindat1.close();

  ifstream eindat2(argFeld[2]);
  if (!eindat2)
  {
    cout << "Error trying to read compare data" << endl;
    return 1;
  }

/******************** special cases (begin) ********************/

  if (strcmp(basename(argFeld[2]), "hess_ex.vgl") == 0)
  {
    while (eindat2.peek()!='I' && !eindat2.eof())
    {
      eindat2.get();
      if (eindat2.eof()) break;
    }
  }
 
  if (strcmp(basename(argFeld[2]), "lop_ex.vgl1") == 0)
  {
    while (eindat2.peek()!='I' && !eindat2.eof())
    {
      eindat2.get();
      if (eindat2.eof()) break;
    }
  }
 
  if (strcmp(basename(argFeld[2]), "lop_ex.vgl2") == 0)
  {
    while (eindat2.peek()!='I' && !eindat2.eof())
    {
      eindat2.get();
      if (eindat2.eof()) break;
    }
  }

/********************* special cases (end) *********************/

  anzvgl = 0;
  while (!eindat2.eof())
  {
    while (eindat2.peek()!='[' && !eindat2.eof()) eindat2.get();
    if (eindat2.eof()) break;

    eindat2 >> vgl[anzvgl];
    anzvgl++;
  } 

  eindat2.close();

  cout << "***********************************************************" << endl;
  cout << "*                   start checking data                   *" << endl;
  cout << "***********************************************************" << endl;
  cout << endl;


  if (anzdaten == anzvgl)
  {
    for (i=0; i<anzdaten; i++)
    {
//      cout << data[i] << " " << vgl[i] << " ";
      cout << data[i] << endl;
      cout << vgl[i] << " ";
      if (data[i] == vgl[i])
      {
        counter[0]++;
        cout << "equal" << endl;
      }
      else
        if (data[i] <= vgl[i])
        {
          counter[1]++;
          cout << "subset" << endl;
        }
        else
          if (data[i] >= vgl[i])
          {
            counter[2]++;
            cout << "superset" << endl;
          }
          else
            {
              counter[3]++;
              cout << "other relation" << endl;
            }
    }

    cout << endl;
    cout << "Number of equalities                 :" << counter[0] << endl;
    cout << "Number of subsets to compare data    :" << counter[1] << endl;
    cout << "Number of supersets to compare data  :" << counter[2] << endl;
    cout << "Other relation                       :" << counter[3] << endl;

#if WINDOWS_X86_32 
    ofstream outfile("ergebnis.dat", ios::app);
#else
	ofstream outfile("./ergebnis.dat", ios::app);
#endif
    outfile << ' ' << counter[0] << ' ' << counter[1] << ' ' << counter[2] << ' ' << counter[3] << ' ' << endl;
    outfile.close();
  }
  else
  {
    cout << "The number of data is not equal in both files !!!" << endl;
    cout << "Number of new data: " << anzdaten << endl;
    cout << "Number of compare data: " << anzvgl << endl;
  }

  cout << endl;
  cout << "************************************************************" << endl;
  cout << "*                           End                            *" << endl;
  cout << "************************************************************" << endl;

  return 0;
}
