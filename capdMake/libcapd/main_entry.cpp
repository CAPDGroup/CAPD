#include <stdexcept>

#ifdef WIN32
#include <windows.h>

int WinMain (HINSTANCE hI, HINSTANCE hPI, LPSTR lpszCmd, int nCS)
{
  throw std::runtime_error("Not implemented");
}

#endif

int mainEntry (int /*argc*/, char */*argv*/ [])
{
  throw std::runtime_error("Not implemented");
}

int mainEntry (void)
{
  throw std::runtime_error("Not implemented");
}
