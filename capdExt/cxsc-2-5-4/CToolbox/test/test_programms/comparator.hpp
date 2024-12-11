/******************************************************************************/
/*                                                                            */
/*                      compares two interval-data-streams                    */
/*                      (based on tester.cpp from C-XSC 2.5.3)                */
/* ComparatorCTests CLASS for C-Tests for CMAKE                               */
/* Klaus-Peter Watzlawek, Bergische Universität Wuppertal, Version 2013-01-22 */
/******************************************************************************/

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <interval.hpp>

using namespace std;
using namespace cxsc;

// Defintion --------------------------------------------------------
class ComparatorCTests {
  private:
    int counter[4];
    interval output_interval[100];
    interval compare_interval[100];
    
    int count_data;
    int count_compare;  
    
    bool debug;
    
    void main_compare();
    void create_intervall(istream &in_istream, 
			  int &in_count_stream, 
			  interval *in_interval,
			  string special = ""
 			);
  public:
    ComparatorCTests();    
    void compare(istream &in_output_file,
		 istream &in_compare_file, 
		 string mode = "");    
    
    int* getCounterStats();
    
    int getResult(istream &in_expected);
    
    void debugmode(bool d);
};

// Implementation --------------------------------------------------

int* ComparatorCTests::getCounterStats() {
    return counter;
}

ComparatorCTests::ComparatorCTests() {
	for (int i = 0; i < 4; i++)
	  counter[i] = 0;
	
	count_data = 0;
	count_compare = 0;	
	debug = 0;
}

void ComparatorCTests::debugmode(bool d) {
	debug = d;
}

void ComparatorCTests::create_intervall(istream &in_istream, 
					int &in_count_stream, 
					interval *in_interval, 
					string special) {
 
  //Special Cases, can be added in future
  if (special.compare("hess_test") == 0 || special.compare("lop_test") == 0) {
    while (in_istream.peek()!='I' && !in_istream.eof())
    {
      in_istream.get();
      if (in_istream.eof()) break;
    }
    
    if (debug) cout << endl << "Special Case Run:" << special << endl;
  }
  
  in_count_stream = 0;
  while (!in_istream.eof()) {
    while (in_istream.peek()!='[' && !in_istream.eof())
      in_istream.get();
    
    if (in_istream.eof()) break;
  
    in_istream >> in_interval[in_count_stream];
    in_count_stream++;   
  }
}

void ComparatorCTests::compare(istream &in_output_file, 
				       istream &in_compare_file, string mode) {
  
  create_intervall(in_output_file,count_data,output_interval,mode);
  create_intervall(in_compare_file,count_compare,compare_interval,mode);
  
  if (debug) {
    cout << "count_data=" << count_data << endl;
    cout << "count_compare=" << count_compare << endl;
  }
  
  main_compare();
}

int ComparatorCTests::getResult(istream &in_expected) {   
  int tmp_counter[4];
  
  //init 0
  for (int i = 0; i < 4; i++)
      tmp_counter[i] = 0;
  
  int i = 0;
  //read from file with expected data
  while (!in_expected.eof())  {    
    in_expected >> tmp_counter[i]; 
    i++;
  } 
  
  //compare calculated data with expected data
  for (int i = 0; i < 4; i++)
      if (tmp_counter[i] != counter[i])
	return 1;
  
    
  return 0;
}

void ComparatorCTests::main_compare() {
  ostringstream debugger;  
  debugger << "***********************************************************" << endl;
  debugger << "*                   start checking data                   *" << endl;
  debugger << "***********************************************************" << endl;
  debugger << endl;


  if (count_data == count_compare)
  {
    for (int i=0; i<count_data; i++)
    {
//      cout << data[i] << " " << vgl[i] << " ";
      debugger << output_interval[i] << endl;
      debugger << compare_interval[i] << " ";
      if (output_interval[i] == compare_interval[i])
      {
        counter[0]++;
        debugger << "equal" << endl;
      }
      else
        if (output_interval[i] <= compare_interval[i])
        {
          counter[1]++;
          debugger << "subset" << endl;
        }
        else
          if (output_interval[i] >= compare_interval[i])
          {
            counter[2]++;
            debugger << "superset" << endl;
          }
          else
            {
              counter[3]++;
              debugger << "other relation" << endl;
            }
    }

    debugger << endl;
    debugger << "Number of equalities                 :" << counter[0] << endl;
    debugger << "Number of subsets to compare data    :" << counter[1] << endl;
    debugger << "Number of supersets to compare data  :" << counter[2] << endl;
    debugger << "Other relation                       :" << counter[3] << endl;


  }
  else
  {
    debugger << "The number of data is not equal in both files !!!" << endl;
    debugger << "Number of new data: " << count_data << endl;
    debugger << "Number of compare data: " << count_compare << endl;
  }

  debugger << endl;
  debugger << "************************************************************" << endl;
  debugger << "*                           End                            *" << endl;
  debugger << "************************************************************" << endl;
  
  if (debug)
    cout << debugger.str();
}

