AUTOMAKE_OPTIONS = gnu
lib_LTLIBRARIES = libsaneLib.la
libsaneLib_la_SOURCES = Corr_preprocess.cpp  map_making.cpp  NoCorr_preprocess.cpp  todprocess.cpp
include_HEADERS = Corr_preprocess.h  map_making.h  NoCorr_preprocess.h  todprocess.h

libsaneLib_la_LIBADD = ../nr/libnr.la ../iniparser/libiniparser.la ../saneIO/libsaneIO.la
INCLUDES = -I ../nr -I ../iniparser -I ../saneIO

libsaneLib_la_LDFLAGS = -version-info 0:0:0

ACLOCAL_AMFLAGS = -I m4