bin_PROGRAMS = saneMerge

saneMerge_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE)
saneMerge_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
	$(SANEIO_LIB)

saneMerge_SOURCES = \
main.cpp \
parse_saneMerge.cpp \
parse_saneMerge.h \
tools.cpp \
tools.h
