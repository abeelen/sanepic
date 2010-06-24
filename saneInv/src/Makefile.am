AUTOMAKE_OPTIONS = gnu
bin_PROGRAMS = saneInv


if MPC
   saneInv_CPPFLAGS= -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK
   ext=_PARA_FRAME
endif

if MPC2
   saneInv_CPPFLAGS= -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK
   ext=_PARA_FRAME
endif

EXEEXT=$(ext)

saneInv_SOURCES = invMatrix.cpp  main.cpp  parseInv.cpp invMatrix.h parseInv.h

LDADD = ../saneIO/libsaneIO.la ../nr/libnr.la ../iniparser/libiniparser.la ../saneLib/libsaneLib.la
INCLUDES = -I../saneIO -I../nr -I../iniparser -I../saneLib
