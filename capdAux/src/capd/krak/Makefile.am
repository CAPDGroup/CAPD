include ${capdMake}/make/common_makefile.mkf

noinst_LTLIBRARIES=libcapdkrak.la
libcapdkrak_la_SOURCES=atom.cpp  chain.cpp  color.cpp  coord3d.cpp  frame.cpp  frstring.cpp  item.cpp  job.cpp  kernel.cpp  noframe.cpp  primitiv.cpp  usermove.cpp
libcapdkrak_la_LDFLAGS=${LIBCAPD_LDFLAGS}
libcapdkrak_la_CXXFLAGS= ${AM_CXXFLAGS}  ${capdkrak_CFLAGS}  ${capdkrak_CFLAGS} ${X11_CFLAGS}
EXTRA_DIST=appl-mac.c dos16.c linux.c ms-dos.c ms-dosch.c ms-win31.c ms-win95.c sun.c sun-os.c sun_v.c unix_ch.c unix-ch.c wxWin.c



pkgconfigdir = @PKG_CONFIG_DIR@
pkgconfig_DATA =
