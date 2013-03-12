bin_PROGRAMS = sanePre

sanePre_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE)
sanePre_LDADD    = $(NR_LIB) $(INIPARSER_LIB)  $(LIBS)  \
	$(SANEIO_LIB)

if PARA
   sanePre_CPPFLAGS += -DUSE_MPI
   ext =_PARA
endif

EXEEXT=$(ext)

sanePre_SOURCES = \
main.cpp
