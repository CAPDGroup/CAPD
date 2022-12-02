/////////////////////////////////////////////////////////////////////////////
//
/// @file ConfigFileReader.h
///
/// @author Tomasz Kapela   @date 2010-10-14
//
/////////////////////////////////////////////////////////////////////////////

#ifndef _CAPD_AUXIL_CONFIG_FILE_READER_H_
#define _CAPD_AUXIL_CONFIG_FILE_READER_H_

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include "capd/auxil/Logger.h"

namespace capd{
  namespace auxil{
/**---------------------------------------------------------------------------
 * The ConfigFileReader parses Windows .ini like file into map<string, string>
 *
 * The file is list of pairs, the key and the value, separated by an equal sign '='.
 * Pair is not permitted to span more than one line.
 *
 * Start of the section is denoted by [SectionName].
 * Section name is added to the following keys name e.g. SectionName.KeyName
 * Next section ends the previous one (use [] to set no section).
 *
 * Commentary :
 * - begins with a hash ('#') as the first non-space character on the line.
 * - begins with semi-colon (;) anywhere
 * and spans to the end of the line
 * 
 * By default it parses input file or stream to the end (EOF).
 * It can be set up to stop on given section name. 
 * It is useful when file contains not only pairs "key = value"
 * but also raw data that spans across several lines.  
 *
 * Example:
 *
   # This is an example
   command = copy

   # In the following we could also define section [source]
   # to shorten keys names
   source.directory = C:\Documents and Settings\Jennifer\My Documents
   source.size = 4   ; number of files
   source.files = file1.cpp file2.cpp file3.cpp file4.cpp  ; file names cannot have white-spaces
   [destination]
   directory = C:\Temp

 *
 *  Usage:
 *
    ConfigFileReader ini("iniFileName");

    std::string command;
    ini.getValue("command", command);

    std::string sourceDir, destinationDir;
    ini.getValue("source.directory", sourceDir);
    ini.getValue("destination.directory", destinationDir);

    int numberOfFiles = 10;                                   // or
    ini.getValue("source.size",numberOfFiles);                //  int numberOfFiles = ini.get<int>("source.size", 10);

    ini.getValue("source.number of files", numberOfFiles);    // it will not change the value because this key do not exist

    std::vector<std::string>  fileNames(numberOfFiles);       
    ini.getValues("source.files", fileNames.begin(), fileNames.end());

 * 
 *  More examples in  capdAux/examples/configFileReader
 * 
*/

class ConfigFileReader : public std::map<std::string, std::string> {
  std::string endSign;
public:
  ConfigFileReader() : endSign(""){
  }
  ConfigFileReader(const std::string & fileName) : endSign(""){
    std::ifstream file(fileName.c_str());
    if(file.good())
      parse(file);
    else
      std::cerr << "Cannon open INI file " << fileName << std::endl;
  }
  ConfigFileReader(std::istream & file) : endSign(""){
    if(file.good())
       parse(file);
    else
          std::cerr << "Cannon read from INI file "<< std::endl;
  }
  void setStopTag(const std::string & tag){
    endSign = tag;
  } 
  const std::string & getStopTag() const {
    return endSign;
  }
  
  // Here is a little convenience method...
  bool isKey(const std::string& s) const {
    return count(s) != 0;
  }

  template <typename T>
  T get(const std::string& key, T defaultValue = T()) {
    if(!isKey(key))
      return defaultValue;
    std::istringstream ss(this->operator [](key));
    T result;
    ss >> result;
    ss.get();  // needed for some types that stop reading after some marker e.g. {232,3223} and do not set eof flag
    if(!ss.eof()) return defaultValue;
    return result;
  }

  template <typename T>
  bool getValue(const std::string & key, T & d) {
    if(!isKey(key)) return false;
    T result;
    std::istringstream ss(this->operator [](key));
    ss >> result;
    ss.get();  // needed for some types that stop reading after some marker e.g. {232,3223} and do not set eof flag
    if(!ss.eof()){
      CAPD_WARN("\n[ Problem when reading parameter '" << key << "' from config file]\n" 
                << " value was : [" << this->operator [](key) << "]\n we read : [" << result << "]"
                << "\n We do not change variable value which is : [" << d << "]");
      return false;
    }
    d = result;
    CAPD_TRACE("parameter read : " << key << " -> ["<< d << "]");
    return true;
  }
  
  bool getValue(const std::string & key, std::string & d) {
      if(!isKey(key)) return false;
      d = (this->operator [](key));
      return true;
  }

  template <typename T>
  bool getValues(const std::string & key, T begin, T end) {
    if(!isKey(key))
      return false;
    std::istringstream ss(this->operator [](key));
    T actual = begin;
    while(actual != end){
      if(ss.eof()){
       CAPD_WARN("\n[ Problem when reading table '" << key << "' from config file]");
        return false;
      }
      ss >> *actual;
      ++actual;
    }
    if(!ss.eof())
      return false;
    return true;
  }
  template <typename T>
  bool getTable(const std::string & key, T table[], size_t size){
    return getValues(key, table, table+size);
  }

  void parse(std::istream& ins) {
    std::string s, key, value;
    CAPD_TRACE(" -- Parsing config file -- ");
    std::string section = "";
    // For each (key, value) pair in the file
    while(std::getline(ins, s)) {
      std::string::size_type begin = s.find_first_not_of(" \f\t\v");

      // Skip blank lines
      if(begin == std::string::npos) continue;

      // Skip commentary (we check if first letter is in given string "#;" )
      if(std::string("#;").find(s[begin]) != std::string::npos) continue;

      // Extract section definition
      if(s[begin]=='['){

        std::string::size_type end = s.find(']');
        // Skip if section tag is not closed
        if(end == std::string::npos) continue;

        // (No leading whitespace allowed)
        begin = s.find_first_not_of(" \f\n\r\t\v", begin+1);
        section = s.substr(begin, end-begin);

        // (No trailing whitespace allowed)
        section.erase(section.find_last_not_of(" \f\t\v") + 1);

        CAPD_TRACE("section [" << section << "]");

        // We check if we should end on that section name 
        if((endSign != "") and (endSign == section))
          break;

        continue;
      }

      // Extract the key value
      std::string::size_type end = s.find('=', begin);
      key = s.substr(begin, end - begin);

      // (No leading or trailing whitespace allowed)
      key.erase(key.find_last_not_of(" \f\t\v") + 1);

      // No blank keys allowed
      if(key.empty()) continue;

      if(section != "")
        key = section + "." + key;

      // Extract the value (no leading or trailing whitespace allowed)
      begin = s.find_first_not_of(" \f\n\r\t\v", end + 1);

      // we skip comments starting with ;
      if((end = s.find(';')) != std::string::npos)
        --end;
      end = s.find_last_not_of(" \f\n\r\t\v", end) + 1;

      value = s.substr(begin, end - begin);

      // Insert the properly extracted (key, value) pair into the map
      (*this)[key] = value;

      CAPD_TRACE("[" << key << "] ==> [" << value << "]");

    }
    CAPD_TRACE("-- end of parsing --");
  }
  // static function required for CAPD_* logger functions in class context.
  CAPD_CLASS_LOGGER
};

//---------------------------------------------------------------------------
// The extraction operator reads configuration::ConfigFileReader until EOF.
// Invalid data is ignored.
//
inline std::istream& operator >>(std::istream& ins, ConfigFileReader& d) {
  d.parse(ins);
  return ins;
}

//---------------------------------------------------------------------------
// The insertion operator writes all configuration::ConfigFileReader to stream.
//
inline std::ostream& operator <<(std::ostream& outs, const ConfigFileReader& d) {
  ConfigFileReader::const_iterator iter;
  for(iter = d.begin(); iter != d.end(); iter++)
    outs << iter->first << " = " << iter->second << std::endl;
  return outs;
}

}} // end of namespace capd::auxil

#endif //_CAPD_AUXIL_CONFIG_FILE_READER_H_
