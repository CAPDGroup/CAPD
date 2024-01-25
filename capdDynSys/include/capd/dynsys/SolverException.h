/////////////////////////////////////////////////////////////////////////////
/// @file SolverException.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_SOLVEREXCEPTION_H_ 
#define _CAPD_DYNSYS_SOLVEREXCEPTION_H_ 

#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

template<typename VectorType>
class SolverException : public std::runtime_error{
public:
  typedef typename VectorType::ScalarType ScalarType;

  VectorType V;
  ScalarType currentTime;
  ScalarType step;
  std::string message;

  SolverException(const std::string &info, const ScalarType& currentTime, const VectorType &_V, const ScalarType &S)
    : std::runtime_error(info), V(_V), currentTime(currentTime), step(S)
  {
    std::ostringstream d;
    d << std::runtime_error::what() << "\n";
    d << "the set: " << V << "\n";
    d << "current time: " << currentTime << "\n";
    d << "time step: " << step << "\n";
    message = d.str();
  }

  ~SolverException() throw() {}

  const char* what() const throw()
  {
    return message.c_str();
  }
};
/// @}

}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_SOLVEREXCEPTION_H_ 


