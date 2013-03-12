bin_PROGRAMS = saneInv

saneInv_CPPFLAGS =  $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE) $(SANELIB_INCLUDE) $(GSL_CFLAGS)
saneInv_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS) \
	$(SANEIO_LIB) $(SANELIB_LIB) $(GSL_LIBS)

if PARA
   saneInv_CPPFLAGS += -DUSE_MPI
   ext=_PARA
endif

EXEEXT=$(ext)

saneInv_SOURCES = \
SaneInvTools.cpp \
SaneInvTools.h \
main.cpp