include ${capdMake}/make/common_programs_makefile.mkf
ACLOCAL_AMFLAGS= -I ../../capdMake/ac_aux/m4

if HAVE_capdkrak
   COND_SUBDIRS=krakdemo
endif

SUBDIRS=auxildemo configFileReader logger $(COND_SUBDIRS)
