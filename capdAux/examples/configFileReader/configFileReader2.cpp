
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


const char fileContents[] = 
"   # This is an example                                                  \n"
"   command = copy                                                        \n"
"                                                                         \n"
"   # In the following we could also define section [source]              \n" 
"   # to shorten keys names                                               \n"
"   source.directory = /Documents and Settings/Jennifer/My Documents      \n"
"   source.size = 4   ; number of files                                   \n"
"   source.files = file1.cpp file2.cpp file3.cpp file4.cpp  ; file names cannot have white-spaces \n"
"   [destination]                                                         \n"
"   directory = C:\\Temp";
           
int main (int /*argc*/, char */*argv*/ [])
{
    // we should init capd logger to debug ConfigFileReader 
    // INIT_CAPD_LOGGER;
     
    // We construct ini reader   
    capd::auxil::ConfigFileReader reader;
    
    // We preapare stream of data (a file, memory etc) and parse it using reader
    std::istringstream file(fileContents);
    file >> reader;

    std::string command;
    reader.getValue("command", command);

    std::string sourceDir, destinationDir;
    reader.getValue("source.directory", sourceDir);
    reader.getValue("destination.directory", destinationDir);

    int numberOfFiles = reader.get<int>("source.size", 10);      // or
                                                                 //reader.getValue("source.size",numberOfFiles); 
    reader.getValue("source.number of files", numberOfFiles);    // it will not change the value because this key does not exist

    std::vector<std::string>  fileNames(numberOfFiles);          
    reader.getValues("source.files", fileNames.begin(), fileNames.end());

    std::cout << command << " from '" << sourceDir << "' to '" << destinationDir << "'\n files :";
    for(int i=0; i<numberOfFiles; ++i)
      std::cout << "\n - " << fileNames[i];
    std::cout << "\n--\n";
	return 0;
} /* main */

