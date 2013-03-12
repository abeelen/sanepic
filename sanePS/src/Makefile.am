bin_PROGRAMS = sanePS

sanePS_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE) $(SANELIB_INCLUDE) $(GSL_CFLAGS)
sanePS_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
	$(SANEIO_LIB) $(SANELIB_LIB) $(GSL_LIBS)

if PARA
   sanePS_CPPFLAGS += -DUSE_MPI
   ext=_PARA
endif

EXEEXT=$(ext)


sanePS_SOURCES = \
sanePSIO.cpp \
sanePSIO.h \
EstimPSSteps.cpp \
EstimPSSteps.h \
main.cpp
