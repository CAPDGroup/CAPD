AUTOMAKE_OPTIONS = subdir-objects

EXTRA_DIST = 
SUBDIRS =

noinst_LTLIBRARIES = libprim.la
libprim_la_LDFLAGS=-Wl,-Bstatic

libprim_la_SOURCES = \
	../ieee/primitive.cpp \
	../interval/stdfun/filib_consts.cpp \
	../fp_traits/fp_traits_double_setup.cpp \
	../fp_traits/fp_traits_float_setup.cpp \
	../fp_traits/fp_traits_base_double_setup.cpp \
	../fp_traits/fp_traits_base_float_setup.cpp \
	../fp_traits/fp_traits_sse_const.cpp \
	../fp_traits/fp_traits_x87_const.cpp
