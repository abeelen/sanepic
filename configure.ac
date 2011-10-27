#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
define([svnversion], esyscmd([sh -c "svnversion|tr -d '\n'"]))dnl
AC_DEFINE(SVN_REVISION, "svnversion", [SVN Revision])
AC_INIT([sanepic],[0.1.svnversion],[sanepic@ias.u-psud.fr])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([-Wall -Werror])

AC_LANG(C++)

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([sanePic/src/main.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
LT_INIT
#AC_OPENMP

# Configure options
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug],
	[enable the verbose/debug output (default is no)])],
	[debug=yes],[debug=no])
if test "x$debug" = xyes; then
AC_SUBST([DEBUG], ["-DDEBUG"])
fi

AC_CONFIG_FILES([
Makefile
nr/Makefile
nr/src/Makefile
iniparser/Makefile
iniparser/src/Makefile
saneIO/Makefile
saneIO/src/Makefile
saneIO/tests/Makefile
saneLib/Makefile
saneLib/src/Makefile
])

AC_ARG_ENABLE([paraframe],
	[AS_HELP_STRING([--enable-paraframe],
	[enable frame parallel processing output (default is no)])],
	[paraframe=yes],[paraframe=no])

AC_ARG_ENABLE([parabolo],
	[AS_HELP_STRING([--enable-parabolo],
	[enable channel parallel processing output (default is no)])],
	[parabolo=yes],[parabolo=no])

AM_CONDITIONAL([PARAFRAME],[test "x$paraframe" = xyes])
AM_CONDITIONAL([PARABOLO],[test "x$parabolo" = xyes])

if test "x$paraframe" = xyes || test "x$parabolo" = xyes; then
   AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
   AC_CHECK_PROGS(MPICXX, mpicxx, $CXX)
   CXX="$MPICXX"
   AC_SUBST(MPICXX)
dnl   AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])
   AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
   CPPFLAGS="$CPPFLAGS -DMPICH_IGNORE_CXX_SEEK"
   LIBS="$LIBS $MPILIBS"
fi

if test -d "sanePS" && test -d "saneCheck" && test -d "saneFix" && test -d "saneMerge" && test -d "saneSplit"; then
   dnl AC_CONFIG_SUBDIRS([sanePS saneCheck saneFix saneMerge saneSplit])
   AC_SUBST([OPT_SUBDIRS],["sanePS saneCheck saneFix saneMerge saneSplit"])
   AC_CONFIG_FILES([saneCheck/Makefile
		    saneCheck/src/Makefile
		    saneFix/Makefile
		    saneFix/src/Makefile
		    saneMerge/Makefile
		    saneMerge/src/Makefile
		    sanePS/Makefile
		    sanePS/src/Makefile
		    saneSplit/Makefile
		    saneSplit/src/Makefile])

fi


# Checks for internal library
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
   AC_SUBST([SANELIB_LIB],["$SANELIBLIB/libsaneLib.la $OPENMP_CXXFLAGS"])
   AC_SUBST([SANELIB_INCLUDE],["-I$SANELIBINC $OPENMP_CXXFLAGS"])
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
   AC_SUBST([SANELIB_LIB],["$SANELIBLIB/libsaneLib.la  $OPENMP_CXXFLAGS"])
   AC_SUBST([SANELIB_INCLUDE],["-I$SANELIBINC  $OPENMP_CXXFLAGS"])
fi



# Checks for lib

AC_CHECK_LIB(m,main,,AC_MSG_ERROR(math needed for GSL))
AC_CHECK_LIB(m,main,,AC_MSG_ERROR(math needed for GSL))
AC_CHECK_LIB(fftw3, fftw_plan_dft_c2r_1d,,AC_MSG_ERROR(fftw3 needed))
AC_CHECK_LIB(gslcblas,main,,AC_MSG_ERROR(GSLCBLAS needed for GSL ))
AC_CHECK_LIB(gsl,main,,AC_MSG_ERROR(GSL needed))
AC_CHECK_LIB(cfitsio, ffvers,,AC_MSG_ERROR(cfitsio needed))
AC_CHECK_LIB(wcs, wcsini,,AC_MSG_ERROR(wcslib needed))
AC_CHECK_LIB(getdata, gd_open,,AC_MSG_ERROR(getdata needed))


# Checks for header files.
AC_CHECK_HEADERS(algorithm,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(cmath,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(cstdio,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(cstdlib,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(cstring,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(fcntl.h,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(fftw3.h,,AC_MSG_WARN(fftw3 header missing ?))
AC_CHECK_HEADERS(fitsio.h,,AC_MSG_WARN(CFITSIO header missing ?))
AC_CHECK_HEADERS(fstream,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(gsl/gsl_math.h,,AC_MSG_WARN(GSL header missing ?))
AC_CHECK_HEADERS(iomanip,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(iostream,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(list,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(ostream,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(sstream,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(stdio.h,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(string,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(sysexits.h,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(sysexits.h,,AC_MSG_WARN(wcslib header missing ?))
AC_CHECK_HEADERS(time.h,,AC_MSG_WARN(time.h header missing ?))
AC_CHECK_HEADERS(unistd.h,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(vector,,AC_MSG_WARN(STL classes header missing ?))
AC_CHECK_HEADERS(wcslib/cel.h,,AC_MSG_WARN(wcslib header missing ?))
AC_CHECK_HEADERS(wcslib/prj.h,,AC_MSG_WARN(wcslib header missing ?))
AC_CHECK_HEADERS(wcslib/wcs.h,,AC_MSG_WARN(wcslib header missing ?))
AC_CHECK_HEADERS(wcslib/wcshdr.h,,AC_MSG_WARN(wcslib header missing ?))
AC_CHECK_HEADERS(wcslib/wcsmath.h,,AC_MSG_WARN(wcslib header missing ?))
AC_CHECK_HEADERS(wcslib/wcstrig.h,,AC_MSG_WARN(wcslib header missing ?))
AC_CHECK_HEADERS(getdata.h,,AC_MSG_WARN(getdata header missing ?))


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

AC_CONFIG_FILES([
saneFrameOrder/Makefile
saneFrameOrder/src/Makefile
saneInv/Makefile
saneInv/src/Makefile
sanePic/Makefile
sanePic/src/Makefile
sanePos/Makefile
sanePos/src/Makefile
sanePre/Makefile
sanePre/src/Makefile
])

AC_CONFIG_MACRO_DIR([m4])
AC_OUTPUT

echo "
 ($PACKAGE_NAME) version $PACKAGE_VERSION
 Prefix............: $prefix
 Debug.............: $debug.
 Para Frame........: $paraframe;
 Para Bolo.........: $parabolo;
 C++ Compiler......: $CXX $CXXFLAGS $CPPFLAGS $OPENMP_CXXFLAGS
 Linker............: $LD $LDFLAGS $LIBS
"
