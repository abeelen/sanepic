AUTOMAKE_OPTIONS = gnu

if MPC
bin_PROGRAMS = sanePic sanePic_MPI sanePic_PARA_BOLO
else
bin_PROGRAMS = sanePic
endif

if MPC
 sanePic_MPI_CPPFLAGS=-DUSE_MPI -DMPICH_IGNORE_CXX_SEEK
sanePic_PARA_BOLO_CPPFLAGS=-DPARA_BOLO -DMPICH_IGNORE_CXX_SEEK
sanePic_MPI_SOURCES = main.cpp sanepic_preprocess.cpp  write_maps_to_disk.cpp sanepic_preprocess.h  write_maps_to_disk.h
sanePic_PARA_BOLO_SOURCES = main.cpp sanepic_preprocess.cpp  write_maps_to_disk.cpp sanepic_preprocess.h  write_maps_to_disk.h
endif

sanePic_SOURCES = main.cpp sanepic_preprocess.cpp  write_maps_to_disk.cpp sanepic_preprocess.h  write_maps_to_disk.h


LDADD = ../saneLib/libsaneLib.la ../saneIO/libsaneIO.la ../nr/libnr.la ../iniparser/libiniparser.la
INCLUDES = -I../saneLib -I../saneIO -I../nr -I../iniparser
