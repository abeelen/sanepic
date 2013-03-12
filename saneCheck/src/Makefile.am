bin_PROGRAMS = saneCheck

saneCheck_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE)
saneCheck_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
	$(SANEIO_LIB)

# if GPROF
#   saneCheck_CPPFLAGS += -pg
# endif	

if PARA
   saneCheck_CPPFLAGS += -DUSE_MPI
   ext=_PARA
endif

EXEEXT=$(ext)

saneCheck_SOURCES = \
main.cpp \
SaneCheckParse.cpp \
SaneCheckParse.h \
SaneCheckTools.cpp \
SaneCheckTools.h 
