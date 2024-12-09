// @todo this file seems unused

#ifndef CAPD_HANDLE_UNEXPECTED_TERMINATE_H
#define CAPD_HANDLE_UNEXPECTED_TERMINATE_H

#include <iostream>
#include <cstdlib>
#include "capd/auxil/ofstreamcout.h"

extern ofstreamcout fcout;

void handle_unexpected() {
  fcout << "unexpected exception thrown" << std::endl;
  exit(1);
}

void handle_terminate() {
  fcout << "terminate() called" << std::endl;
  exit(1);
}

#endif // CAPD_HANDLE_UNEXPECTED_TERMINATE_H


