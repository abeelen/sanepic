bin_PROGRAMS = sanePic

sanePic_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE) $(SANELIB_INCLUDE)
sanePic_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS)  \
	$(SANEIO_LIB) $(SANELIB_LIB)

if PARA
   sanePic_CPPFLAGS += -DUSE_MPI
   ext=_PARA
endif

EXEEXT=$(ext)

sanePic_SOURCES = \
main.cpp \
SanePicIO.cpp \
SanePicIO.h
