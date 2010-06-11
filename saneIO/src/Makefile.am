AUTOMAKE_OPTIONS = gnu
lib_LTLIBRARIES = libsaneIO.la
libsaneIO_la_SOURCES = covMatrix_IO.cpp imageIO.cpp dataIO.cpp inline_IO2.cpp inputFileIO.cpp mpi_architecture_builder.cpp parser_functions.cpp psdIO.cpp
include_HEADERS = covMatrix_IO.h dataIO.h imageIO.h inline_IO2.h inputFileIO.h mpi_architecture_builder.h parser_functions.h psdIO.h struct_definition.h

libsaneIO_la_LIBADD = ../nr/libnr.la ../iniparser/libiniparser.la
INCLUDES = -I ../nr -I ../iniparser

libsaneIO_la_LDFLAGS = -version-info 0:0:0

ACLOCAL_AMFLAGS = -I m4
