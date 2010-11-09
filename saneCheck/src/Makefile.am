bin_PROGRAMS = saneCheck

saneCheck_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE)
saneCheck_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
	$(SANEIO_LIB)

if PARAFRAME
   saneCheck_CPPFLAGS += -DPARA_FRAME
   ext=_PARA_FRAME
endif

if PARABOLO
   saneCheck_CPPFLAGS += -DPARA_FRAME
   ext=_PARA_FRAME
endif

EXEEXT=$(ext)


saneCheck_SOURCES = \
main.cpp \
parse_saneCheck.cpp \
parse_saneCheck.h \
tools.cpp \
tools.h 
