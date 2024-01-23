#include "capd/autodiff/DagIndexer.hpp"

#include "capd/multiPrec/mplib.h"
#include "capd/intervals/mplib.h"

template class capd::autodiff::DagIndexer<capd::multiPrec::MpReal>;
template class capd::autodiff::DagIndexer<capd::MpInterval>;
