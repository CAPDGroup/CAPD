

#ifndef _UNITTESTS_GRAPHICALSETTEST_CPP_
#define _UNITTESTS_GRAPHICALSETTEST_CPP_


#include "capd/capdlib.h"
#include "capd/dynset/C0GraphicalSet.h"
class Output{
public:
  void show(capd::C0Set & set){}
};
typedef capd::dynset::C0GraphicalSet<capd::C0Rect2Set, Output>  SetType;

#define FIXTURE_NAME Rect2SetTest

#include "AffineSetCommonTest.hpp"


#endif