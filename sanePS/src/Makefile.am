AUTOMAKE_OPTIONS = gnu
bin_PROGRAMS = sanePS



if MPC
   sanePS_CPPFLAGS= -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK
   ext=_PARA_FRAME
endif

if MPC2
   sanePS_CPPFLAGS= -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK
   ext=_PARA_FRAME
endif

EXEEXT=$(ext)

sanePS_SOURCES = main.cpp  estimPS.cpp estimPS.h estimPS_steps.cpp estimPS_steps.h

LDADD = ../saneIO/libsaneIO.la ../nr/libnr.la ../iniparser/libiniparser.la ../saneLib/libsaneLib.la
INCLUDES = -I../saneIO -I../nr -I../iniparser -I../saneLib


