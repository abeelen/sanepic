bin_PROGRAMS = saneFrameOrder

saneFrameOrder_LDADD = $(NR_LIB) $(LIBS) \
	$(INIPARSER_LIB) $(SANEIO_LIB)
saneFrameOrder_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) \
	$(SANEIO_INCLUDE)

saneFrameOrder_CPPFLAGS += -DPARA_FRAME

EXEEXT=_PARA_FRAME

saneFrameOrder_SOURCES = \
main.cpp \
parse_FBFO.cpp \
parse_FBFO.h
