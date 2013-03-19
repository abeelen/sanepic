bin_PROGRAMS = saneFix

saneFix_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE)
saneFix_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
	$(SANEIO_LIB)

if PARA
   saneFix_CPPFLAGS += -DUSE_MPI
   ext=_PARA
endif

EXEEXT=$(ext)

saneFix_SOURCES = \
main.cpp \
SaneFixTools.cpp \
SaneFixTools.h \
SaneFixParse.cpp \
SaneFixParse.h