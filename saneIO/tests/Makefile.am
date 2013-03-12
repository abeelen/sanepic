# check_PROGRAMS = parserTest
# TESTS = $(check_PROGRAMS)
EXTRA_PROGRAMS = parserTest dirfileTest
#bin_PROGRAMS = parserTest dirfileTest

parserTest_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
        $(SANEIO_INCLUDE)

parserTest_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
        $(SANEIO_LIB)

dirfileTest_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
        $(SANEIO_INCLUDE)

dirfileTest_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
        $(SANEIO_LIB)

if PARA
   parserTest_CPPFLAGS  += -DUSE_MPI
   dirfileTest_CPPFLAGS += -DUSE_MPI
   ext=_PARA
endif

EXEEXT=$(ext)

parserTest_SOURCES = \
parserTest.cpp

dirfileTest_SOURCES = \
dirfileTest.cpp
