AUTOMAKE_OPTIONS = gnu
bin_PROGRAMS = sanePre

if MPC
   sanePre_CPPFLAGS = -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK
   ext=_PARA_FRAME
endif

if PARAC 
   sanePre_CPPFLAGS = -DPARA_BOLO -DMPICH_IGNORE_CXX_SEEK
   ext=_PARA_BOLO
endif

EXEEXT=$(ext)

sanePre_SOURCES = main.cpp

LDADD = ../saneLib/libsaneLib.la ../saneIO/libsaneIO.la ../nr/libnr.la ../iniparser/libiniparser.la
INCLUDES = -I../saneLib -I../saneIO -I../nr -I../iniparser
