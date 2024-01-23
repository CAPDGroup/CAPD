/////////////////////////////////////////////////////////////////////////////
/// @file Dll.h
///
/// @author Mateusz Juda <mateusz.juda@gmail.com>
///
/// @date 2016-11-24
/////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2016 by the CAPD Group.
//
// This file constitutes a part of the CAPD library (capdAux),
// distributed under the terms of the GNU General Public License.
// Consult http://capd.ii.uj.edu.pl and  http://redhom.ii.edu.pl/ for details.
/////////////////////////////////////////////////////////////////////////////

#ifndef CAPD_FILE_CAPDAUX_AUXIL_DLL_H
#define CAPD_FILE_CAPDAUX_AUXIL_DLL_H

// Code from https://gcc.gnu.org/wiki/Visibility
// Macros for symbols fvisibility. Read about -fvisibility compiler option.

#define DLL_PUBLIC __attribute__((visibility("default")))
#define DLL_LOCAL __attribute__((visibility("hidden")))

// #if defined _WIN32 || defined __CYGWIN__
// #ifdef BUILDING_DLL
// #ifdef __GNUC__
// #define DLL_PUBLIC __attribute__((dllexport))
// #else
// #define DLL_PUBLIC __declspec(dllexport)  // Note: actually gcc seems to also supports this syntax.
// #endif
// #else
// #ifdef __GNUC__
// #define DLL_PUBLIC __attribute__((dllimport))
// #else
// #define DLL_PUBLIC __declspec(dllimport)  // Note: actually gcc seems to also supports this syntax.
// #endif
// #endif
// #define DLL_LOCAL
// #else
// #if __GNUC__ >= 4
// #define DLL_PUBLIC __attribute__((visibility("default")))
// #define DLL_LOCAL __attribute__((visibility("hidden")))
// #else
// #define DLL_PUBLIC
// #define DLL_LOCAL
// #endif
// #endif

#endif  // CAPD_FILE_CAPDAUX_AUXIL_DLL_H
