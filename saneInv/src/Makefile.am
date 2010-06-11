AUTOMAKE_OPTIONS = gnu
bin_PROGRAMS = saneInv

saneInv_SOURCES = invMatrix.cpp  main.cpp  parseInv.cpp
saneInv_LDADD = ../saneIO/libsaneIO.la ../nr/libnr.la ../iniparser/libiniparser.la
INCLUDES = -I../saneIO -I../nr -I../iniparser
