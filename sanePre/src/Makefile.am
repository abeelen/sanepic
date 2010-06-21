AUTOMAKE_OPTIONS = gnu

if MPC
   bin_PROGRAMS = sanePre sanePre_MPI sanePre_PARA_BOLO
else
   bin_PROGRAMS = sanePre
endif


if MPC
   sanePre_MPI_CPPFLAGS= -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK
   sanePre_PARA_BOLO_CPPFLAGS= -DPARA_BOLO -DMPICH_IGNORE_CXX_SEEK
   sanePre_MPI_SOURCES = main.cpp 
   sanePre_PARA_BOLO_SOURCES = main.cpp 
   sanePre_SOURCES = main.cpp 
else
   sanePre_SOURCES = main.cpp 
endif


LDADD = ../saneLib/libsaneLib.la ../saneIO/libsaneIO.la ../nr/libnr.la ../iniparser/libiniparser.la
INCLUDES = -I../saneLib -I../saneIO -I../nr -I../iniparser
