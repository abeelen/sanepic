#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([saneFrameOrder], [0.3], [sanepic@ias.u-psud.fr])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
LT_PREREQ([2.2])
LT_INIT([dlopen])

AC_LANG(C++)

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL

AX_MPI
CXX="$MPICXX"
LIBS="-DUSE_MPI -DMPICH_IGNORE_CXX_SEEK $MPILIBS $LIBS"

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS(algorithm,,AC_MSG_WARN(STL classes missing ?))
AC_CHECK_HEADERS(cstdio,,AC_MSG_WARN(STL classes missing ?))
AC_CHECK_HEADERS(cstdlib,,AC_MSG_WARN(STL classes missing ?))
AC_CHECK_HEADERS(iostream,,AC_MSG_WARN(STL classes missing ?))
AC_CHECK_HEADERS(iomanip,,AC_MSG_WARN(STL classes missing ?))
AC_CHECK_HEADERS(fstream,,AC_MSG_WARN(STL classes missing ?))
AC_CHECK_HEADERS(vector,,AC_MSG_WARN(STL classes missing ?))
AC_CHECK_HEADERS(string,,AC_MSG_WARN(STL classes missing ?))

AC_CHECK_HEADERS(time.h,,AC_MSG_WARN(STL classes missing ?))
AC_CHECK_HEADERS(fcntl.h,,AC_MSG_WARN(STL classes missing ?))

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for internal library functions.
if test -d "$srcdir/nr"; then
   NRINC='$(top_srcdir)/nr/src'
   NRLIB='$(top_builddir)/nr/src'
else
   if test -d "$srcdir/../nr"; then
      NRINC='$(top_srcdir)/../nr/src'
      NRLIB='$(top_builddir)/../nr/src'
   fi
fi
if test -n "$NRLIB"; then
   AC_SUBST([NR_LIB],["$NRLIB/libnr.la"])
   AC_SUBST([NR_INCLUDE],["-I$NRINC"])
fi

if test -d "$srcdir/iniparser"; then
   INIPARSERINC='$(top_srcdir)/iniparser/src'
   INIPARSERLIB='$(top_builddir)/iniparser/src'
else
   if test -d "$srcdir/../iniparser"; then
      INIPARSERINC='$(top_srcdir)/../iniparser/src'
      INIPARSERLIB='$(top_builddir)/../iniparser/src'
   fi
fi
if test -n "$INIPARSERLIB"; then
   AC_SUBST([INIPARSER_LIB],["$INIPARSERLIB/libiniparser.la"])
   AC_SUBST([INIPARSER_INCLUDE],["-I$INIPARSERINC"])
fi

if test -d "$srcdir/saneIO"; then
   SANEIOINC='$(top_srcdir)/saneIO/src'
   SANEIOLIB='$(top_builddir)/saneIO/src'
else
   if test -d "$srcdir/../saneIO"; then
      SANEIOINC='$(top_srcdir)/../saneIO/src'
      SANEIOLIB='$(top_builddir)/../saneIO/src'
   fi
fi
if test -n "$SANEIOLIB"; then
   AC_SUBST([SANEIO_LIB],["$SANEIOLIB/libsaneIO.la"])
   AC_SUBST([SANEIO_INCLUDE],["-I$SANEIOINC"])
fi
 
# Configure DEBUG source code, if requested.
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug],
	[enable the verbose/debug output (default is no)])],
	[debug=yes],[debug=no])

if test "x$debug" = xyes; then
   AC_SUBST([DEBUG],["-DDEBUG"])
fi

AC_CONFIG_FILES([Makefile
		 src/Makefile])
AC_OUTPUT

echo "
 ($PACKAGE_NAME) version $PACKAGE_VERSION
 Prefix............: $prefix
 Debug.............: $debug
 C++ Compiler......: $CXX $CXXFLAGS $CPPFLAGS
 Linker............: $LD $LDFLAGS $LIBS
"
