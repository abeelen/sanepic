lib_LTLIBRARIES = libsaneLib.la
libsaneLib_la_LDFLAGS = -version-info 0:3:0

libsaneLib_la_LIBADD = $(NR_LIB) $(INIPARSER_LIB) $(DEBUG)  \
	$(SANEIO_LIB)
libsaneLib_la_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) \
	$(SANEIO_INCLUDE) $(OPENMP_CXXFLAGS)


libsaneLib_la_SOURCES = \
CorrPreprocess.cpp \
FitPoly.cpp \
MapMaking.cpp \
NoCorrPreprocess.cpp \
TodProcess.cpp \
CorrPreprocess.h \
FitPoly.h \
MapMaking.h \
NoCorrPreprocess.h \
TodProcess.h
