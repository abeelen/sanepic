bin_PROGRAMS = sanePS

sanePS_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE) $(SANELIB_INCLUDE)
sanePS_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
	$(SANEIO_LIB) $(SANELIB_LIB)

if PARAFRAME
   sanePS_CPPFLAGS += -DPARA_FRAME -DUSE_MPI
   ext=_PARA_FRAME
endif

if PARABOLO
   sanePS_CPPFLAGS += -DPARA_FRAME -DUSE_MPI
   ext=_PARA_FRAME
endif

EXEEXT=$(ext)


sanePS_SOURCES = \
estimPS.cpp \
estimPS.h \
estimPS_steps.cpp \
estimPS_steps.h \
main.cpp
