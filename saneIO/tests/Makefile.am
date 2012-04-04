# check_PROGRAMS = parserTest
# TESTS = $(check_PROGRAMS)
EXTRA_PROGRAMS = parserTest

parserTest_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
        $(SANEIO_INCLUDE)

parserTest_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
        $(SANEIO_LIB)

if PARAFRAME
   parserTest_CPPFLAGS += -DPARA_FRAME -DUSE_MPI
   ext=_PARA_FRAME
endif

if PARABOLO
   parserTest_CPPFLAGS += -DPARA_FRAME -DUSE_MPI
   ext=_PARA_FRAME
endif

parserTest_SOURCES = \
parserTest.cpp
