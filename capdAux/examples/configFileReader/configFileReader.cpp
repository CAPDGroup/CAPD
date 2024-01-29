
/////////////////////////////////////////////////////////////////////////////
///
/// @file configFileReader.cpp
///
/// @author Tomasz Kapela
///
/////////////////////////////////////////////////////////////////////////////

// This file constitutes a part of the Homology Library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.edu.pl/ for details. 

// Started on September, 2013. 

#include <sstream>
#include <vector>
#include "capd/auxil/ConfigFileReader.h"

const char fileContents[] = "\
   # Number of entries in a table                                        \n\
   number_of_entries = 10                                                \n\
   cols = 3  ; number of columns                                         \n\
   rows = 3  ; number of rows                                            \n\
   [data]                                                                \n\
   table = 1 2 3 4 5 6 7 8 9 10                                          \n\
   [raw_data]                                                            \n\
   2.3   12.1 12.1 \n\
   123.3  2.3  3.3 \n\
   33.4  22.3  2.2 \n\
   [next]          \n\
   message = \"The numbers are\" "; 
           
int main (int /*argc*/, char */*argv*/ [])
{
     
    // We construct ini reader   and set it up to stop parsing on raw_data section 
    capd::auxil::ConfigFileReader reader;
    reader.setStopTag("raw_data");
    
    // We preapare stream of data (a file, memory etc) 
    std::istringstream file(fileContents);
    
    // We parse it up to [raw_data] section
    file >> reader;

    // We read raw data section that contains array of given number of rows an columns
    int cols = reader.get<int>("cols",0);
    int rows = reader.get<int>("rows",0);
    for(int i = 0; i<cols*rows; ++i){
      double d;
      file >> d;
      std::cout << d << " ";
    }
    std::cout << std::endl;
    
    // We parse the rest of the file (assuming that it does not containd second raw_data section) 
    file>> reader;
    
    std::vector<int> table(reader.get<int>("number_of_entries",0));
    reader.getValues("data.table", table.begin(), table.end());
    
    std::string message;
    reader.getValue("next.message", message);
    
    std::cout << "\n" << message << "\n";
    for(int i=0; i<reader.get<int>("number_of_entries",0); ++i)
      std::cout << table[i] << " | ";
    std::cout << "\n---\n";
   	return 0;
} /* main */


