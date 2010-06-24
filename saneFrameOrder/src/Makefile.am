AUTOMAKE_OPTIONS = gnu
bin_PROGRAMS = saneFrameOrder


saneFrameOrder_CPPFLAGS= -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK
saneFrameOrder_SOURCES = main.cpp parse_FBFO.cpp parse_FBFO.h

LDADD = ../saneIO/libsaneIO.la ../nr/libnr.la ../iniparser/libiniparser.la
INCLUDES = -I../saneIO -I../nr -I../iniparser
