bin_PROGRAMS = saneSplit

saneSplit_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE)
saneSplit_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
	$(SANEIO_LIB)

saneSplit_SOURCES = \
main.cpp \
SaneSplitParse.cpp \
SaneSplitParse.h \
SaneSplitTools.cpp \
SaneSplitTools.h
