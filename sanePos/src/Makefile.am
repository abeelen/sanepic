AUTOMAKE_OPTIONS = gnu
bin_PROGRAMS = sanePos

sanePos_SOURCES = main.cpp  parsePos.cpp  sanePos_map_making.cpp  sanePos_preprocess.cpp
sanePos_LDADD = ../saneIO/libsaneIO.la ../nr/libnr.la ../iniparser/libiniparser.la

INCLUDES = -I../saneIO -I../nr -I../iniparser
