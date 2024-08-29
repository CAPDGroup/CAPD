/*
 * DagIndexer.cpp
 *
 *  Created on: Aug 25, 2012
 *      Author: kapela
 */


#include "capd/autodiff/DagIndexer.hpp"
#include "capd/intervals/lib.h"

template class capd::autodiff::DagIndexer<double>;
template class capd::autodiff::DagIndexer<long double>;
template class capd::autodiff::DagIndexer<capd::DInterval>;
