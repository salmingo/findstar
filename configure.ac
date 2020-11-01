dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(findstar, 1.0)


AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE()

AC_PROG_CXX
AC_PROG_CC

AC_ARG_ENABLE(debug,
AS_HELP_STRING([--enable-debug],
               [enable debugging, default: no]),
[case "${enableval}" in
             yes) debug=true ;;
             no)  debug=false ;;
             *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
esac],
[debug=false])

AM_CONDITIONAL(DEBUG, test x"$debug" = x"true")

AC_CHECK_HEADER_STDBOOL
AC_CHECK_FUNCS([memset pow sqrt])
AC_FUNC_MALLOC

AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT

