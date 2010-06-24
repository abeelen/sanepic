AUTOMAKE_OPTIONS = gnu
bin_PROGRAMS = sanePic

if MPC
   sanePic_CPPFLAGS = -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK
   ext=_PARA_FRAME
endif

if PARAC 
   sanePic_CPPFLAGS = -DPARA_BOLO -DMPICH_IGNORE_CXX_SEEK
   ext=_PARA_BOLO
endif

EXEEXT=$(ext)

sanePic_SOURCES = main.cpp sanepic_preprocess.cpp  write_maps_to_disk.cpp sanepic_preprocess.h  write_maps_to_disk.h

LDADD = ../saneLib/libsaneLib.la ../saneIO/libsaneIO.la ../nr/libnr.la ../iniparser/libiniparser.la
INCLUDES = -I../saneLib -I../saneIO -I../nr -I../iniparser
