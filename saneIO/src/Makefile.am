lib_LTLIBRARIES = libsaneIO.la
libsaneIO_la_LDFLAGS = -version-info 0:3:0

libsaneIO_la_LIBADD = $(NR_LIB) $(INIPARSER_LIB) $(GSL_LIBS)
libsaneIO_la_CPPFLAGS = $(NR_INCLUDE) $(INIPARSER_INCLUDE) $(DEBUG) $(GSL_CFLAGS)

if PARA
  libsaneIO_la_CPPFLAGS += -DUSE_MPI
endif

libsaneIO_la_SOURCES = \
CovMatrixIO.cpp \
Crc.cpp \
DataIO.cpp \
ImageIO.cpp \
InputFileIO.cpp \
MPIConfiguration.cpp \
ParserFunctions.cpp \
TemporaryIO.cpp \
Utilities.cpp \
CovMatrixIO.h \
Crc.h \
DataIO.h \
ErrorCode.h \
ImageIO.h \
InputFileIO.h \
MPIConfiguration.h \
ParserFunctions.h \
StructDefinition.h \
TemporaryIO.h \
Utilities.h