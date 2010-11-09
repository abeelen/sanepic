#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([sanePS], [0.3], [sanepic@ias.u-psud.fr])
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


# Configure options
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug],
	[enable the verbose/debug output (default is no)])],
	[debug=yes],[debug=no])
if test "x$debug" = xyes; then
AC_SUBST([DEBUG], ["-DDEBUG"])
fi

AC_ARG_ENABLE([paraframe],
	[AS_HELP_STRING([--enable-paraframe],
	[enable frame parallel processing output (default is no)])],
	[paraframe=yes],[paraframe=no])

AC_ARG_ENABLE([parabolo],
	[AS_HELP_STRING([--enable-parabolo],
	[enable channel parallel processing output (default is no)])],
	[parabolo=yes],[parabolo=no])

if test "x$paraframe" = xyes || test "x$parabolo" = xyes; then
   AX_MPI
   CXX="$MPICXX"
   CPPFLAGS="$CPPFLAGS -DMPICH_IGNORE_CXX_SEEK"
fi

AM_CONDITIONAL([PARAFRAME],[test "x$paraframe" = xyes])
AM_CONDITIONAL([PARABOLO],[test "x$parabolo" = xyes])

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
if test -n "$SANEIOLIB"; then
   AC_SUBST([SANEIO_LIB],["$SANEIOLIB/libsaneIO.la"])
   AC_SUBST([SANEIO_INCLUDE],["-I$SANEIOINC"])
fi

if test -d "$srcdir/saneLib"; then
   SANELIBINC='$(top_srcdir)/saneLib/src'
   SANELIBLIB='$(top_builddir)/saneLib/src'
else
   if test -d "$srcdir/../saneLib"; then
      SANELIBINC='$(top_srcdir)/../saneLib/src'
      SANELIBLIB='$(top_builddir)/../saneLib/src'
   fi
fi
if test -n "$SANELIBLIB"; then
   AC_SUBST([SANELIB_LIB],["$SANELIBLIB/libsaneLib.la"])
   AC_SUBST([SANELIB_INCLUDE],["-I$SANELIBINC"])
fi

if test -d "$srcdir/saneLib"; then
   SANELIBINC='$(top_srcdir)/saneLib/src'
   SANELIBLIB='$(top_builddir)/saneLib/src'
else
   if test -d "$srcdir/../saneLib"; then
      SANELIBINC='$(top_srcdir)/../saneLib/src'
      SANELIBLIB='$(top_builddir)/../saneLib/src'
   fi
fi
if test -n "$SANELIBLIB"; then
   AC_SUBST([SANELIB_LIB],["$SANELIBLIB/libsaneLib.la"])
   AC_SUBST([SANELIB_INCLUDE],["-I$SANELIBINC"])
fi

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS(cstring,,AC_MSG_WARN(STL classes missing ?))
AC_CHECK_HEADERS(fcntl.h,,AC_MSG_WARN(STL classes missing ?))
AC_CHECK_HEADERS(list,,AC_MSG_WARN(STL classes missing ?))
AC_CHECK_HEADERS(stdio.h,,AC_MSG_WARN(STL classes missing ?))
AC_CHECK_HEADERS(time.h,,AC_MSG_WARN(STL classes missing ?))
AC_CHECK_HEADERS(wcslib/cel.h,,AC_MSG_WARN(wcslib missing ?))
AC_CHECK_HEADERS(wcslib/prj.h,,AC_MSG_WARN(wcslib missing ?))
AC_CHECK_HEADERS(wcslib/wcshdr.h,,AC_MSG_WARN(wcslib missing ?))
AC_CHECK_HEADERS(wcslib/wcs.h,,AC_MSG_WARN(wcslib missing ?))
AC_CHECK_HEADERS(wcslib/wcsmath.h,,AC_MSG_WARN(wcslib missing ?))
AC_CHECK_HEADERS(wcslib/wcstrig.h,,AC_MSG_WARN(wcslib missing ?))



# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

AC_CONFIG_FILES([Makefile
		 src/Makefile])
AC_OUTPUT

echo "
 ($PACKAGE_NAME) version $PACKAGE_VERSION
 Prefix............: $prefix
 Debug.............: $debug.
 Para Frame........: $paraframe;
 Para Bolo.........: $parabolo;
 C++ Compiler......: $CXX $CXXFLAGS $CPPFLAGS
 Linker............: $LD $LDFLAGS $LIBS
"
