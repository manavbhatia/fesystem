dnl $Id: configure.ac,v 1.8.4.1 2007-03-04 03:37:11 manav Exp $

 
dnl Process this file with autoconf to produce a configure script.


AC_INIT(FESystem, 1.0, manav@u.washington.edu)
AC_CONFIG_SRCDIR([src/fem_main.C])


AC_MSG_RESULT()
AC_MSG_RESULT(-----------------   FESystem  ----------------- )
AC_MSG_RESULT(----------------- configuration --------------- )
AC_MSG_RESULT()

AM_INIT_AUTOMAKE

dnl find out the information for variable 'host'
AC_CANONICAL_HOST


AM_CONFIG_HEADER(src/FESystem/FESystemConfig.h)

AC_PROG_CC
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_PROG_FC
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS

AC_ARG_ENABLE(shared, 
		AC_HELP_STRING([--enable-shared],
				[build shared libraries for all packages]),
		enableshared=$enableval,
		enableshared=yes)
AC_SUBST(enableshared)

dnl sunstitute the variables that for this package make system
AC_SUBST(host)
AC_SUBST(abs_top_srcdir)

AC_SUBST(CXX)
AC_SUBST(CC)
AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(FCFLAGS)
AC_SUBST(FCLIBS)

AC_SUBST(PETSC_ARCH)
AC_SUBST(PETSC_DIR)
AC_SUBST(LIBMESH_ARCH)
AC_SUBST(SLEPC_DIR)
AC_SUBST(HDF5_DIR)

dnl configure the packaged contined in this distribution
AC_CONFIG_SUBDIRS([contrib/deal.II contrib/libmesh])


dnl now build the output
AC_MSG_RESULT()
AC_MSG_RESULT(----------------- creating output --------------- )
AC_MSG_RESULT()
AC_OUTPUT([	Makefile build/Makefile 
		contrib/Makefile
		contrib/ARPACK/Makefile])

AC_MSG_RESULT()
AC_MSG_RESULT(---------------------   FESystem  --------------------- )
AC_MSG_RESULT(----------------- configuration complete --------------- )
AC_MSG_RESULT()

dnl configure.ac ends here
