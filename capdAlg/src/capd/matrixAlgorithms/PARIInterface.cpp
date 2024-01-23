/////////////////////////////////////////////////////////////////////////////
/// @file PARIInterface.cpp
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2014-06-17
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD::RedHom Group.
//
// This file constitutes a part of the CAPD::RedHom library,
// distributed under the terms of the GNU General Public License.
// Consult  http://redhom.ii.edu.pl/ for details.

#include "PARIInterface.hpp"

using namespace capd::matrixAlgorithms;

#ifdef HAVE_PARI

const bool PARIInterface::_enabled = true;


PARIInterface::PARIInterface()
{
  pari_init(100000000, 0);
}

PARIInterface::~PARIInterface()
{
  pari_close();
}

#else

const bool PARIInterface::_enabled = false;

PARIInterface::PARIInterface()
{
}

PARIInterface::~PARIInterface()
{
}

#endif



    PARI_INTERFACE_SMITH_FORM_INSTANCE(short);
    PARI_INTERFACE_SMITH_FORM_INSTANCE(int);
    PARI_INTERFACE_SMITH_FORM_INSTANCE(long);
    typedef long long llong;
    PARI_INTERFACE_SMITH_FORM_INSTANCE(llong);
