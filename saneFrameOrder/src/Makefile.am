bin_PROGRAMS = saneFrameOrder

saneFrameOrder_LDADD = $(NR_LIB) $(LIBS) \
	$(INIPARSER_LIB) $(SANEIO_LIB)
saneFrameOrder_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE)

saneFrameOrder_CPPFLAGS += -DUSE_MPI

EXEEXT=_PARA

saneFrameOrder_SOURCES = \
main.cpp \
FrameOrder.cpp \
FrameOrder.h