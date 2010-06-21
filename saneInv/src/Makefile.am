AUTOMAKE_OPTIONS = gnu

if MPC
bin_PROGRAMS = saneInv saneInv_MPI
else
bin_PROGRAMS = saneInv
endif

if MPC
saneInv_MPI_CPPFLAGS=-DUSE_MPI -DMPICH_IGNORE_CXX_SEEK
saneInv_MPI_SOURCES = invMatrix.cpp  main.cpp  parseInv.cpp invMatrix.h parseInv.h
endif

saneInv_SOURCES = invMatrix.cpp  main.cpp  parseInv.cpp invMatrix.h parseInv.h

LDADD = ../saneIO/libsaneIO.la ../nr/libnr.la ../iniparser/libiniparser.la ../saneLib/libsaneLib.la
INCLUDES = -I../saneIO -I../nr -I../iniparser -I../saneLib
