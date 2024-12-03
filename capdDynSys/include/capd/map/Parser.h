/// @addtogroup map
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Parser.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#ifndef _CAPD_MAP_PARSER_H_
#define _CAPD_MAP_PARSER_H_

#include <string>
#include <vector>
#include <map>
#include "capd/autodiff/NodeType.h"

namespace capd{
namespace map{

int parseVariables(std::string expression, std::vector<std::string>& var);
int parseExpression(
    std::string& expression,
    const std::vector<std::string>& var,
    std::vector<capd::autodiff::Node>& dag,
    std::map<std::string,int>& knownNodes
);
int parseMap(unsigned numberOfVariables, std::string expression, const std::vector<std::string>& var, std::vector<capd::autodiff::Node>& dag, std::vector<int>& pos);

size_t parseBrackets(const std::string &e, size_t position = std::string::npos);
size_t searchForOperator(const std::string &e, unsigned char op, size_t position = std::string::npos);
size_t searchForFunction(const std::string &fun, const std::string &eq);
bool searchForFunction(const std::string &fun, const std::string &eq, std::string & params);
std::string & removeBrackets(std::string &eq);
void splitVariables(const std::string &, const std::string&,std::vector<std::string>& result);
bool isConstant(std::string&, double& value);
/// Converts given text to double (returns true on success)
bool stringToDouble(std::string const & text, double & result);
void removeWhiteSpaces(std::string & text);
void optimizeDAG(std::vector<capd::autodiff::Node>& dag, std::vector<int>& pos);
}} // namespace capd::map

#endif // _CAPD_MAP_PARSER_H_

/// @}
