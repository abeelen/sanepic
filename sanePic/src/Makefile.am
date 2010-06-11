AUTOMAKE_OPTIONS = gnu
bin_PROGRAMS = sanePic

sanePic_SOURCES = main.cpp  parseSanepic.cpp  sanepic_preprocess.cpp  write_maps_to_disk.cpp
sanePic_LDADD = ../saneLib/libsaneLib.la ../saneIO/libsaneIO.la ../nr/libnr.la ../iniparser/libiniparser.la

INCLUDES = -I../saneLib -I../saneIO -I../nr -I../iniparser
