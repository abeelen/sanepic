lib_LTLIBRARIES = libsaneIO.la
libsaneIO_la_LDFLAGS = -version-info 0:3:0

libsaneIO_la_LIBADD = $(NR_LIB) $(INIPARSER_LIB)
libsaneIO_la_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG)

if PARAFRAME
  libsaneIO_la_CPPFLAGS += -DUSE_MPI
endif

if PARABOLO
  libsaneIO_la_CPPFLAGS += -DUSE_MPI
endif

libsaneIO_la_SOURCES = \
covMatrix_IO.cpp \
crc.cpp \
dataIO.cpp \
imageIO.cpp \
inputFileIO.cpp \
mpi_architecture_builder.cpp \
parser_functions.cpp \
temporary_IO.cpp \
utilities.cpp \
covMatrix_IO.h \
crc.h \
dataIO.h \
error_code.h \
imageIO.h \
inputFileIO.h \
mpi_architecture_builder.h \
parser_functions.h \
struct_definition.h \
temporary_IO.h \
utilities.h