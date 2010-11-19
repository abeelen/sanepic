check_PROGRAMS = parserTest
TESTS = $(check_PROGRAMS)

parserTest_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
        $(SANEIO_INCLUDE)

parserTest_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
        $(SANEIO_LIB)

if PARAFRAME
   parserTest_CPPFLAGS += -DPARA_FRAME
   ext=_PARA_FRAME
endif

if PARABOLO
   parserTest_CPPFLAGS += -DPARA_BOLO
   ext=_PARA_BOLO
endif

parserTest_SOURCES = \
parserTest.cpp