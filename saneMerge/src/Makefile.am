bin_PROGRAMS = saneMerge

saneMerge_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE)
saneMerge_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
	$(SANEIO_LIB)

saneMerge_SOURCES = \
main.cpp \
SaneMergeParse.cpp \
SaneMergeParse.h \
SaneMergeTools.cpp \
SaneMergeTools.h
