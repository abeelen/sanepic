lib_LTLIBRARIES = libsaneLib.la
libsaneLib_la_LDFLAGS = -version-info 0:3:0

libsaneLib_la_LIBADD = $(NR_LIB) $(INIPARSER_LIB) $(DEBUG)  \
	$(SANEIO_LIB)
libsaneLib_la_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) \
	$(SANEIO_INCLUDE) 


libsaneLib_la_SOURCES = \
cholesky.cpp \
Corr_preprocess.cpp \
fitpoly.cpp \
map_making.cpp \
NoCorr_preprocess.cpp \
todprocess.cpp \
cholesky.h \
Corr_preprocess.h \
fitpoly.h \
map_making.h \
NoCorr_preprocess.h \
todprocess.h
