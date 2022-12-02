#ifndef _EXAMPLE_H
#define _EXAMPLE_H 1

#include "capd/auxil/Logger.h"

struct Example
{

  void log()
  {
    CAPD_INFO("Member function");

    CAPD_TRACE("msg: TRACE");
    CAPD_DEBUG("msg: DEBUG");
    CAPD_INFO("msg: INFO");
    CAPD_WARN("msg: WARN");
    CAPD_ERROR("msg: ERROR");
    CAPD_FATAL("msg: FATAL");
  }

  static void logStatic()
  {
    CAPD_INFO("Static function");
    CAPD_TRACE("msg: TRACE");
    CAPD_DEBUG("msg: DEBUG");
    CAPD_INFO("msg: INFO");
    CAPD_WARN("msg: WARN");
    CAPD_ERROR("msg: ERROR");
    CAPD_FATAL("msg: FATAL");
  }

private:

  // static function required for CAPD_* logger functions in class context.
  CAPD_CLASS_LOGGER;
};

void globalFunction()
{
  CAPD_INFO("Global function");
  CAPD_TRACE("msg: TRACE");
  CAPD_DEBUG("msg: DEBUG");
  CAPD_INFO("msg: INFO");
  CAPD_WARN("msg: WARN");
  CAPD_ERROR("msg: ERROR");
  CAPD_FATAL("msg: FATAL");
}

#endif /* _EXAMPLE_H */
