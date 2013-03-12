bin_PROGRAMS = sanePos

sanePos_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE) $(SANELIB_INCLUDE)
sanePos_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
	$(SANEIO_LIB) $(SANELIB_LIB)

if PARA
   sanePos_CPPFLAGS += -DUSE_MPI
   ext=_PARA
endif

EXEEXT=$(ext)

sanePos_SOURCES = \
main.cpp \
SanePosMapMaking.cpp \
SanePosMapMaking.h \
SanePosPreProcess.cpp \
SanePosPreProcess.h \
Coord.cpp \
Coord.h