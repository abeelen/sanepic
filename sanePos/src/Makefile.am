AUTOMAKE_OPTIONS = gnu

if MPC
bin_PROGRAMS = sanePos sanePos_MPI
else
bin_PROGRAMS = sanePos
endif


if MPC
sanePos_MPI_SOURCES = main.cpp  sanePos_map_making.cpp  sanePos_preprocess.cpp sanePos_map_making.h  sanePos_preprocess.h
sanePos_MPI_CPPFLAGS=-DUSE_MPI -DMPICH_IGNORE_CXX_SEEK
endif

sanePos_SOURCES = main.cpp  sanePos_map_making.cpp  sanePos_preprocess.cpp sanePos_map_making.h  sanePos_preprocess.h

LDADD = ../saneIO/libsaneIO.la ../nr/libnr.la ../iniparser/libiniparser.la
INCLUDES = -I../saneIO -I../nr -I../iniparser


