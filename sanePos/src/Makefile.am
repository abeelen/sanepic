bin_PROGRAMS = sanePos

sanePos_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE)
sanePos_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
	$(SANEIO_LIB)

if PARAFRAME
   sanePos_CPPFLAGS += -DPARA_FRAME
   ext=_PARA_FRAME
endif

if PARABOLO
   sanePos_CPPFLAGS += -DPARA_FRAME
   ext=_PARA_FRAME
endif

EXEEXT=$(ext)


sanePos_SOURCES = \
main.cpp \
sanePos_map_making.cpp \
sanePos_map_making.h \
sanePos_preprocess.cpp \
sanePos_preprocess.h
