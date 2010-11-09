#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([sanepic], [0.3], [sanepic@ias.u-psud.fr])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
LT_PREREQ([2.2])
LT_INIT([dlopen])

# Configure options
AC_ARG_ENABLE([debug],
	[AS_HELP_STRING([--enable-debug],
	[enable the verbose/debug output (default is no)])],
	[debug=yes],[debug=no])

AC_ARG_ENABLE([paraframe],
	[AS_HELP_STRING([--enable-paraframe],
	[enable frame parallel processing output (default is no)])],
	[paraframe=yes],[paraframe=no])

AC_ARG_ENABLE([parabolo],
	[AS_HELP_STRING([--enable-parabolo],
	[enable channel parallel processing output (default is no)])],
	[parabolo=yes],[parabolo=no])

if test "x$paraframe" = xyes || test "x$parabolo" = xyes; then
   AC_CONFIG_SUBDIRS([saneFrameOrder])
fi

if test -d "sanePS" && test -d "saneCheck" && test -d "saneFix" && test -d "saneMerge" && test -d "saneSplit"; then
   AC_CONFIG_SUBDIRS([sanePS saneCheck saneFix saneMerge saneSplit])
   AC_SUBST([OPT_SUBDIRS],["sanePS saneCheck saneFix saneMerge saneSplit"])
fi

AM_CONDITIONAL([PARAFRAME],[test "x$paraframe" = xyes])
AM_CONDITIONAL([PARABOLO],[test "x$parabolo" = xyes])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SUBDIRS([nr iniparser saneIO saneLib saneInv sanePos sanePre sanePic])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

echo "
 ($PACKAGE_NAME) version $PACKAGE_VERSION
 Prefix............: $prefix
 Debug.............: $debug
 Para Frame........: $paraframe;
 Para Bolo.........: $parabolo;
"
