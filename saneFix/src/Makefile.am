bin_PROGRAMS = saneFix

saneFix_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE)
saneFix_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
	$(SANEIO_LIB)

if PARAFRAME
   saneFix_CPPFLAGS += -DPARA_FRAME
   ext=_PARA_FRAME
endif

if PARABOLO
   saneFix_CPPFLAGS += -DPARA_FRAME
   ext=_PARA_FRAME
endif

EXEEXT=$(ext)


saneFix_SOURCES = \
main.cpp \
parse_saneFix.cpp \
parse_saneFix.h \
tools.cpp \
tools.h
