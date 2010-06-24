AUTOMAKE_OPTIONS = gnu
bin_PROGRAMS = sanePos



if MPC
   sanePos_CPPFLAGS= -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK
   ext=_PARA_FRAME
endif

if MPC2
   sanePos_CPPFLAGS= -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK
   ext=_PARA_FRAME
endif

EXEEXT=$(ext)

sanePos_SOURCES = main.cpp  sanePos_map_making.cpp  sanePos_preprocess.cpp sanePos_map_making.h  sanePos_preprocess.h

LDADD = ../saneIO/libsaneIO.la ../nr/libnr.la ../iniparser/libiniparser.la
INCLUDES = -I../saneIO -I../nr -I../iniparser


