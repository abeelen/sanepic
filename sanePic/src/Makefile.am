bin_PROGRAMS = sanePic

sanePic_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE) $(SANELIB_INCLUDE)
sanePic_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS)  \
	$(SANEIO_LIB) $(SANELIB_LIB)

if PARAFRAME
   sanePic_CPPFLAGS += -DPARA_FRAME -DUSE_MPI
   ext=_PARA_FRAME
endif

if PARABOLO
   sanePic_CPPFLAGS += -DPARA_BOLO -DUSE_MPI
   ext=_PARA_BOLO
endif

EXEEXT=$(ext)

sanePic_SOURCES = \
main.cpp \
write_maps_to_disk.cpp \
write_maps_to_disk.h
