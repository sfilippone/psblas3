dnl $Id$

dnl
dnl Michele Martone's proposal for a transition towards the GNU build system for the PSBLAS library
dnl

dnl Process this file with autoconf to produce a configure script.
dnl
dnl usage : aclocal -I config/ && autoconf && ./configure && make 

dnl then  : VAR=VAL ./configure

dnl In some configurations (AIX) the next line is needed:
dnl MPIFC=mpxlf95 ./configure

dnl then  : ./configure VAR=VAL 
dnl then  : ./configure --help=short
dnl then  : ./configure --help

dnl the PSBLAS modules get this task difficult to accomplish!
dnl SEE  : --module-path --include-path

dnl NOTE : There is no cross compilation support.



###############################################################################
###############################################################################
# 
#  This script is used by the PSBLAS to determine the compilers, linkers, and
# libraries to build its libraries executable code.
#  Its behaviour is driven on the compiler it finds or it is dictated to work
# with.
#
###############################################################################
###############################################################################

# NOTE: the literal for version (the second argument to AC_INIT should be a literal!)
AC_INIT([PSBLAS],3.4, salvatore.filippone@uniroma2.it)

# VERSION is the file containing the PSBLAS version code
# FIXME
psblas_cv_version="3.4"

# A sample source file
AC_CONFIG_SRCDIR([base/modules/psb_base_mod.f90])
dnl AC_CONFIG_AUX_DIR(config)

# Our custom M4 macros are in the 'config' directory
AC_CONFIG_MACRO_DIR([config])
AC_MSG_NOTICE([
--------------------------------------------------------------------------------
	Welcome to the $PACKAGE_NAME $psblas_cv_version configure Script.

	This creates Make.inc, but if you read carefully the
	documentation, you can make your own by hand for your needs.

	Be sure to specify the library paths of your interest. Examples:
        ./configure --with-libs=-L/some/directory/LIB    <- will append to LIBS
        FC=mpif90 F77=$FC CC=mpicc ./configure           <- will force FC,CC,F77

	See ./configure --help=short fore more info.
--------------------------------------------------------------------------------
		])
###############################################################################
# Compilers detection: FC,F77,CC should be set, if found.
###############################################################################
#
# Installation. 
#
#
AC_PROG_INSTALL

AC_MSG_CHECKING([where to install])
case $prefix in  
   \/* )   eval "INSTALL_DIR=$prefix";;
   * ) eval "INSTALL_DIR=/usr/local/psblas";;
esac
case $libdir in 
   \/* )   eval "INSTALL_LIBDIR=$libdir";;
   * ) eval "INSTALL_LIBDIR=$INSTALL_DIR/lib";;
esac
case $includedir in 
   \/* )   eval "INSTALL_INCLUDEDIR=$includedir";;
   * ) eval "INSTALL_INCLUDEDIR=$INSTALL_DIR/include";;
esac
case $docsdir in 
   \/* )   eval "INSTALL_DOCSDIR=$docsdir";;
   * ) eval "INSTALL_DOCSDIR=$INSTALL_DIR/docs";;
esac
case $samplesdir in 
   \/* )   eval "INSTALL_SAMPLESDIR=$samplesdir";;
   * ) eval "INSTALL_SAMPLESDIR=$INSTALL_DIR/samples";;
esac
AC_MSG_RESULT([$INSTALL_DIR $INSTALL_INCLUDEDIR $INSTALL_LIBDIR $INSTALL_DOCSDIR $INSTALL_SAMPLESDIR])

# Note that the following line won't save from troubles.
AC_PROG_FC([ftn xlf2003_r xlf2003 xlf95_r xlf95 xlf90 xlf pgf95 pgf90 ifort ifc  nagfor gfortran])
AC_PROG_F77([ftn xlf pgf77 ifort ifc nagfor gfortran])
AC_PROG_CC([xlc pgcc icc gcc cc ])
dnl AC_PROG_CXX

dnl AC_PROG_F90 doesn't exist, at the time of writing this !
dnl AC_PROG_F90

# Sanity checks, although redundant (useful when debugging this configure.ac)!
if test "X$FC" == "X" ; then
	AC_MSG_ERROR([Problem : No Fortran compiler specified nor found!])
fi

if test "X$CC" == "X" ; then
	AC_MSG_ERROR([Problem : No C compiler specified nor found!])
fi

###############################################################################
#			Suitable MPI compilers detection
###############################################################################
# Note: Someday we will contemplate a fake MPI - configured version of PSBLAS
###############################################################################
# First check whether the user required our serial (fake) mpi.
PAC_ARG_SERIAL_MPI

#Note : we miss the name of the Intel C compiler
if test x"$pac_cv_serial_mpi" == x"yes" ; then
   FAKEMPI="fakempi.o";
   MPIFC="$FC";
   MPIF77="$F77";
   MPICC="$CC";
else 
AC_LANG([C])
if test "X$MPICC" = "X" ; then
    # This is our MPICC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPICC],[mpxlc mpcc mpicc cc])
fi
ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for C]])])


AC_LANG([Fortran])

if test "X$MPIFC" = "X" ; then
    # This is our MPIFC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPIFC],[mpxlf2003_r mpxlf2003 mpxlf95_r mpxlf90 mpf95 mpf90 mpifort mpif95 mpif90 ftn ])
fi

ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for Fortran]])])

AC_LANG(Fortran 77)
if test "X$MPIF77" = "X" ; then
    # This is our MPIFC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPIF77],[mpxlf mpf77 mpif77 ftn])
fi
ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for Fortran 77]])])
FC="$MPIFC" ;
F77="$MPIF77";
CC="$MPICC";
fi

# We leave a default language for the next checks.
dnl AC_LANG([Fortran 77])
AC_LANG([C])

dnl Now on, MPIFC should be set, as MPIF77 and MPICC

###############################################################################
# Sanity checks, although redundant (useful when debugging this configure.ac)!
###############################################################################

if test "X$MPIFC" == "X" ; then
	AC_MSG_ERROR([Problem : No MPI Fortran compiler specified nor found!])
fi

if test "X$MPICC" == "X" ; then
	AC_MSG_ERROR([Problem : No MPI C compiler specified nor found!])
fi

###############################################################################
# 			FLAGS and LIBS user customization
###############################################################################

dnl NOTE : no spaces before the comma, and no brackets before the second argument!
PAC_ARG_WITH_FLAGS(ccopt,CCOPT)
PAC_ARG_WITH_FLAGS(fcopt,FCOPT)
#PAC_ARG_WITH_FLAGS(f90copt,F90COPT)
#PAC_ARG_WITH_FLAGS(ldflags,LDFLAGS)
PAC_ARG_WITH_LIBS
PAC_ARG_WITH_FLAGS(clibs,CLIBS)
PAC_ARG_WITH_FLAGS(flibs,FLIBS)

dnl candidates for removal:
PAC_ARG_WITH_FLAGS(library-path,LIBRARYPATH)
PAC_ARG_WITH_FLAGS(include-path,INCLUDEPATH)
PAC_ARG_WITH_FLAGS(module-path,MODULE_PATH)

# we just gave the user the chance to append values to these variables

###############################################################################

dnl  Library oriented Autotools facilities (we don't care about this for now)
dnl AC_PROG_LIBTOOL
dnl AM_MAINTAINER_MODE
AC_PROG_RANLIB
dnl system's ranlib will be found, too

AM_INIT_AUTOMAKE
dnl Specify required version of autoconf.
AC_PREREQ(2.59)


###############################################################################
#	Compiler identification (sadly, it is necessary)
###############################################################################
psblas_cv_fc=""

dnl Do we use gfortran & co ? Compiler identification.
dnl NOTE : in /autoconf/autoconf/fortran.m4 there are plenty of better tests!

PAC_CHECK_HAVE_GFORTRAN(
	[psblas_cv_fc="gcc"],
)

PAC_CHECK_HAVE_CRAYFTN(
	[psblas_cv_fc="cray"],
)

if test x"$psblas_cv_fc" == "x" ; then
 if eval "$MPIFC -qversion 2>&1 | grep XL 2>/dev/null" ; then
	psblas_cv_fc="xlf"
	# Some configurations of the XLF want "-WF," prepended to -D.. flags.
	# TODO : discover the exact conditions when the usage of -WF is needed.
	psblas_cv_define_prepend="-WF,"
        if eval "$MPIFC -qversion 2>&1 | grep -e\"Version: 10\.\" 2>/dev/null"; then
        	FDEFINES="$psblas_cv_define_prepend-DXLF_10 $FDEFINES"		  
        fi

 # Note : there could be problems with old xlf compiler versions ( <10.1 )
 # since (as far as it is known to us) -WF, is not used in earlier versions.
 # More problems could be undocumented yet.
 elif eval "$MPIFC -V  2>&1 | grep Sun 2>/dev/null" ; then
   # Sun compiler detection
  
   psblas_cv_fc="sun"
 elif eval "$MPIFC -V  2>&1 | grep Portland 2>/dev/null" ; then
   # Portland group compiler detection
   
   psblas_cv_fc="pg"
 elif eval "$MPIFC -V 2>&1 | grep Intel.*Fortran.*Compiler 2>/dev/null" ; then
   # Intel compiler identification
    
   psblas_cv_fc="ifc"
 elif eval "$MPIFC -v  2>&1 | grep NAG 2>/dev/null" ; then
   psblas_cv_fc="nag"
   FC="$MPIFC"
   F77="$MPIFC"
 else
   psblas_cv_fc=""
  # unsupported MPI Fortran compiler
     AC_MSG_NOTICE([[Unknown Fortran compiler, proceeding with fingers crossed !]])
 fi
fi
if test "X$psblas_cv_fc" == "Xgcc" ; then
PAC_HAVE_MODERN_GFORTRAN(
	[],
        [AC_MSG_ERROR([Sorry, we require GNU Fortran 4.6 or later.])]
)
fi
# TODO : SEE _AC_PROG_FC_V
# TODO : AC_MSG_ERROR(see "$ac_cv_prog_FC_fc")
# AC_MSG_NOTICE( "ac_cv_prog_FC_fc : $ac_cv_prog_FC_fc")
###############################################################################

###############################################################################
#	Linking, symbol mangling, and misc tests
###############################################################################

# Note : This is functional to Make.inc rules and structure (see below).
AC_LANG([C])
AC_CHECK_SIZEOF(void *)
# Define for platforms with 64 bit (void *     ) pointers
if test X"$ac_cv_sizeof_void_p" == X"8" ; then
  CDEFINES="-DPtr64Bits $CDEFINES"
fi
AC_LANG([Fortran])
if test "X$psblas_cv_fc" == X"pg" ; then
   save_FC=$FC
   FC=$F77
fi
__AC_FC_NAME_MANGLING
if test "X$psblas_cv_fc" == X"pg" ; then
   FC=$save_FC
fi
AC_LANG([C])
dnl AC_MSG_NOTICE([Fortran name mangling: $ac_cv_fc_mangling])
[pac_fc_case=${ac_cv_fc_mangling%%,*}]
[pac_fc_under=${ac_cv_fc_mangling#*,}]
[pac_fc_sec_under=${pac_fc_under#*,}]
[pac_fc_sec_under=${pac_fc_sec_under# }]
[pac_fc_under=${pac_fc_under%%,*}]
[pac_fc_under=${pac_fc_under# }]
AC_MSG_CHECKING([defines for C/Fortran name interfaces])
if test "x$pac_fc_case" == "xlower case"; then
   if test "x$pac_fc_under" == "xunderscore"; then
      if test "x$pac_fc_sec_under" == "xno extra underscore"; then 
          pac_f_c_names="-DLowerUnderscore"
      elif test "x$pac_fc_sec_under" == "xextra underscore"; then
          pac_f_c_names="-DLowerDoubleUnderscore"
      else
          pac_f_c_names="-DUNKNOWN"
dnl  AC_MSG_NOTICE([Fortran name mangling extra underscore unknown case])	
      fi
   elif test "x$pac_fc_under" == "xno underscore"; then
      pac_f_c_names="-DLowerCase"
   else		       
      pac_f_c_names="-DUNKNOWN"
dnl  AC_MSG_NOTICE([Fortran name mangling underscore unknown case])	
   fi
elif test "x$pac_fc_case" == "xupper case"; then
   if test "x$pac_fc_under" == "xunderscore"; then
      if test "x$pac_fc_sec_under" == "xno extra underscore"; then 
          pac_f_c_names="-DUpperUnderscore"
      elif test "x$pac_fc_sec_under" == "xextra underscore"; then
          pac_f_c_names="-DUpperDoubleUnderscore"
      else
          pac_f_c_names="-DUNKNOWN"
dnl  AC_MSG_NOTICE([Fortran name mangling extra underscore unknown case])	
      fi
   elif test "x$pac_fc_under" == "xno underscore"; then
      pac_f_c_names="-DUpperCase"
   else		       
      pac_f_c_names="-DUNKNOWN"
dnl  AC_MSG_NOTICE([Fortran name mangling underscore unknown case])	
   fi
dnl  AC_MSG_NOTICE([Fortran name mangling UPPERCASE not handled])	
else
   pac_f_c_names="-DUNKNOWN"   
dnl   AC_MSG_NOTICE([Fortran name mangling unknown case])	
fi
CDEFINES="$pac_f_c_names $CDEFINES"

AC_MSG_RESULT([ $pac_f_c_names ])

###############################################################################
# 			Make.inc generation logic
###############################################################################
F90COPT="$FCOPT"

if test "X$CCOPT" == "X" ; then
 if test "X$psblas_cv_fc" == "Xgcc" ; then
 # note that no space should be placed around the equality symbol in assignements
 # Note : 'native' is valid _only_ on GCC/x86  (32/64 bits)
   CCOPT="-O3  $CCOPT"

 elif test "X$psblas_cv_fc" == X"xlf" ; then
   # XL compiler : consider using -qarch=auto
     CCOPT="-O3 -qarch=auto $CCOPT"
 elif test "X$psblas_cv_fc" == X"ifc" ; then
   # other compilers ..
     CCOPT="-O3 $CCOPT"
 elif test "X$psblas_cv_fc" == X"pg" ; then
   # other compilers ..
     CCOPT="-fast $CCOPT"
  # NOTE : PG & Sun use -fast instead -O3
 elif test "X$psblas_cv_fc" == X"sun" ; then
   # other compilers ..
     CCOPT="-fast $CCOPT"
 elif test "X$psblas_cv_fc" == X"cray" ; then
     CCOPT="-O3 $CCOPT"
     MPICC="cc"
 elif test "X$psblas_cv_fc" == X"nag" ; then
   # using GCC in conjunction with NAG.
     CCOPT="-O2"
 else 
     CCOPT="-O2 $CCOPT"
 fi
fi
#CFLAGS="${CCOPT}"

if test "X$FCOPT" == "X" ; then
 if test "X$psblas_cv_fc" == "Xgcc" ; then
 # note that no space should be placed around the equality symbol in assignations
 # Note : 'native' is valid _only_ on GCC/x86  (32/64 bits)
   FCOPT="-O3  $FCOPT"
 elif test "X$psblas_cv_fc" == X"xlf" ; then
   # XL compiler : consider using -qarch=auto
     FCOPT="-O3 -qarch=auto  -qfixed -qsuffix=f=f:cpp=F -qlanglvl=extended $FCOPT"
     FCFLAGS="-qhalt=e $FCFLAGS"
 elif test "X$psblas_cv_fc" == X"ifc" ; then
   # other compilers ..
     FCOPT="-O3 $FCOPT"
 elif test "X$psblas_cv_fc" == X"pg" ; then
   # other compilers ..
     FCOPT="-fast $FCOPT"
  # NOTE : PG & Sun use -fast instead -O3
 elif test "X$psblas_cv_fc" == X"sun" ; then
   # other compilers ..
     FCOPT="-fast $FCOPT"
 elif test "X$psblas_cv_fc" == X"cray" ; then
     FCOPT="-O3 -em $FCOPT"
 elif test "X$psblas_cv_fc" == X"nag" ; then
   # NAG  compiler ..
     FCOPT="-O2 "
  # NOTE : PG & Sun use -fast instead -O3
 else 
     FCOPT="-O2 $FCOPT"
 fi
fi
if test "X$psblas_cv_fc" == X"nag" ; then
   # Add needed options 
   FCOPT="$FCOPT  -dcfuns -f2003 -wmismatch=mpi_scatterv,mpi_alltoallv,mpi_gatherv,mpi_allgatherv"
fi
#FFLAGS="${FCOPT}"

if test "X$F90COPT" == "X" ; then
 if test "X$psblas_cv_fc" == "Xgcc" ; then
 # note that no space should be placed around the equality symbol in assignations
 # Note : 'native' is valid _only_ on GCC/x86  (32/64 bits)
   F90COPT="-O3  $F90COPT"
 elif test "X$psblas_cv_fc" == X"xlf" ; then
   # XL compiler : consider using -qarch=auto
     F90COPT="-O3 -qarch=auto -qsuffix=f=f90:cpp=F90 -qlanglvl=extended $F90COPT"
 elif test "X$psblas_cv_fc" == X"ifc" ; then
   # other compilers ..
     F90COPT="-O3 $F90COPT"
 elif test "X$psblas_cv_fc" == X"pg" ; then
   # other compilers ..
     F90COPT="-fast $F90COPT"
 elif test "X$psblas_cv_fc" == X"sun" ; then
     F90COPT="-fast $F90COPT"
 elif test "X$psblas_cv_fc" == X"cray" ; then
     MPIFC="ftn"
     F90COPT="-O3 -em $F90COPT"
 elif test "X$psblas_cv_fc" == X"nag" ; then
   # NAG compiler
     F90COPT="-O2"
 else 
   # other compilers ..
     F90COPT="-O2 $F90COPT"
 fi
else
 echo "Found FCFLAGS $F90COPT"
  #F90COPT="${FCFLAGS}"
fi
if test "X$psblas_cv_fc" == X"nag" ; then
   # Add needed options 
     F90COPT="$F90COPT -dcfuns -f2003 -wmismatch=mpi_scatterv,mpi_alltoallv,mpi_gatherv,mpi_allgatherv"
     EXTRA_OPT="-mismatch_all"
fi

#FCFLAGS="${F90COPT}"

# COPT,FCOPT, F90COPT are aliases for FFLAGS,CFLAGS,FCFLAGS .

##############################################################################
# Compilers variables selection 
##############################################################################
F90=${FC}
MPF90=${MPIFC}
FC=${FC}
MPF77=${MPIFC}
CC=${CC}
MPCC=${MPICC}



##############################################################################
# Choice of our compilers, needed by Make.inc
##############################################################################
if test "X$FLINK" == "X" ; then
   FLINK=${MPF77}
fi

if test "X$F90LINK" == "X" ; then
   F90LINK=${MPF90}
fi

##############################################################################
#
#      Include flags and Module file names for dependencies.
#
##############################################################################
BASEMODNAME=psb_base_mod
PRECMODNAME=psb_prec_mod
METHDMODNAME=psb_krylov_mod
UTILMODNAME=psb_util_mod

if test "X$psblas_cv_fc" == X"cray" 
then
   MODEXT=".mod"
   FMFLAG="-I"
   FIFLAG="-I"
   BASEMODNAME=PSB_BASE_MOD
   PRECMODNAME=PSB_PREC_MOD
   METHDMODNAME=PSB_KRYLOV_MOD
   UTILMODNAME=PSB_UTIL_MOD

else 
   AX_F90_MODULE_EXTENSION
   AX_F90_MODULE_FLAG
   MODEXT=".$ax_cv_f90_modext"
   FMFLAG="${ax_cv_f90_modflag%%[ ]*}"
   FIFLAG=-I
fi
###############################################################################

###############################################################################
# NOTE :
# Missing stuff : 
# 		In the case the detected fortran compiler is ifort, icc or gcc
#		should be valid options.
#		The same for pg (Portland Group compilers).
###############################################################################


###############################################################################
# Custom test : do we have a module or include for MPI Fortran interface?
if test x"$pac_cv_serial_mpi" == x"yes" ; then
   FDEFINES="$psblas_cv_define_prepend-DSERIAL_MPI $psblas_cv_define_prepend-DMPI_MOD $FDEFINES";
else 
 PAC_FORTRAN_CHECK_HAVE_MPI_MOD_F08()
 if test x"$pac_cv_mpi_f08" == x"yes" ; then
dnl    FDEFINES="$psblas_cv_define_prepend-DMPI_MOD_F08 $FDEFINES";
    FDEFINES="$psblas_cv_define_prepend-DMPI_MOD $FDEFINES";
  else					     	     
    PAC_FORTRAN_CHECK_HAVE_MPI_MOD(
 	  [FDEFINES="$psblas_cv_define_prepend-DMPI_MOD $FDEFINES"],
	  [FDEFINES="$psblas_cv_define_prepend-DMPI_H $FDEFINES"])	  
 fi
fi

PAC_ARG_LONG_INTEGERS
if test x"$pac_cv_long_integers" == x"yes" ; then
   FDEFINES="$psblas_cv_define_prepend-DLONG_INTEGERS $FDEFINES";
fi

#
# Tests for support of various Fortran features; some of them are critical,
# some optional
#

#
# Critical features
#
PAC_FORTRAN_TEST_TR15581( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for TR15581. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.6.])]
)

PAC_FORTRAN_TEST_EXTENDS( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for EXTENDS. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.6.])]
)

PAC_FORTRAN_TEST_CLASS_TBP( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for CLASS and type bound procedures. 
 Please get a Fortran compiler that supports them, e.g. GNU Fortran 4.6.])]
)

PAC_FORTRAN_TEST_SOURCE( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for SOURCE= allocation. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.6.])]
)

PAC_FORTRAN_HAVE_MOVE_ALLOC(
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_MOVE_ALLOC $FDEFINES"],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for MOVE_ALLOC. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.6.])]
)

PAC_FORTRAN_HAVE_SUBMODULE(
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for SUBMODULES. 
 Please get a Fortran compiler that supports it.])]
)

PAC_FORTRAN_TEST_ISO_C_BIND( 
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_ISO_C_BINDING $FDEFINES"],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for ISO_C_BINDING. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.6.])]
)

#
# Optional features
#

PAC_FORTRAN_TEST_VOLATILE( 
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_VOLATILE $FDEFINES"],
)

PAC_FORTRAN_TEST_GENERICS( 
	[],
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_BUGGY_GENERICS $FDEFINES"]
)

PAC_FORTRAN_TEST_FLUSH( 
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_FLUSH_STMT $FDEFINES"],
)

PAC_FORTRAN_TEST_ISO_FORTRAN_ENV(
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_ISO_FORTRAN_ENV $FDEFINES"],
)

PAC_FORTRAN_TEST_FINAL( 
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_FINAL $FDEFINES"],
)

PAC_FORTRAN_TEST_SAME_TYPE( 
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_SAME_TYPE_AS $FDEFINES"],
)

PAC_FORTRAN_TEST_EXTENDS_TYPE( 
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_EXTENDS_TYPE_OF $FDEFINES"],
)

PAC_FORTRAN_TEST_MOLD( 
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_MOLD $FDEFINES"],
)

###############################################################################
#	Additional pathname stuff (yes, it is redundant and confusing...)
###############################################################################

# -I
if test x"$INCLUDEPATH" != "x" ; then
	FINCLUDES="$FINCLUDES $INCLUDEPATH"
	CINCLUDES="$CINCLUDES $INCLUDEPATH"
fi

# -L
if test x"$LIBRARYPATH" != "x" ; then
	FINCLUDES="$FINCLUDES $LIBRARYPATH"
fi

# -I
if test x"$MODULE_PATH" != "x" ; then
	FINCLUDES="$FINCLUDES $MODULE_PATH"
fi

###############################################################################
#	Parachute rules for ar and ranlib ... (could cause problems)
###############################################################################

if test "X$AR" == "X" ; then
  AR="ar"
fi

if test "X$RANLIB" == "X" ; then
  RANLIB="ranlib"
fi

# This should be portable
AR="${AR} -cur"


###############################################################################
#		BLAS  library presence checks
###############################################################################

# Note : The libmkl.a (Intel Math Kernel Library) library could be used, too.
# It is sufficient to specify it as -lmkl in the CLIBS or FLIBS or LIBS 
# and specify its path adjusting -L/path in CFLAGS.

# Right now it is a matter of user's taste when linking custom applications.
# But PSBLAS examples could take advantage of these libraries, too.

PAC_BLAS([], [AC_MSG_ERROR([[Cannot find BLAS library, specify a path using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])
PAC_LAPACK(
[FDEFINES="$psblas_cv_define_prepend-DHAVE_LAPACK $FDEFINES"],
)
AC_LANG([C])


###############################################################################
#		BLACS library presence checks
###############################################################################
#AC_LANG([C])
#if test x"$pac_cv_serial_mpi" == x"no" ; then
#save_FC="$FC";
#save_CC="$CC";
#FC="$MPIFC";
#CC="$MPICC";
#PAC_CHECK_BLACS
#FC="$save_FC";
#CC="$save_CC";
#fi

PAC_MAKE_IS_GNUMAKE

###############################################################################
#		METIS, SuperLU, SuperLU_Dist UMFPack libraries presence checks
###############################################################################

# Note : also SuperLUStat,superlu_gridexit,.. 
# Note : SuperLU_Dist produces superlu.a by default..
# Note : THESE WERE NOT TESTED
#AC_CHECK_LIB(superlu,[superlu_malloc_dist],psblas_cv_have_superludist=yes,psblas_cv_have_superludist=no)
#AC_CHECK_LIB(superlu,[superlu_malloc],psblas_cv_have_superlu=yes,psblas_cv_have_superlu=no)
# Note : also umfdi_local_search, ...
#AC_CHECK_LIB(umf,umfpack_di_solve,psblas_cv_have_umfpack=yes,psblas_cv_have_umfpack=no,[amd])

AC_ARG_WITH(rsb, AC_HELP_STRING([--with-rsb], [Specify Recursive Sparse BLAS library linkage info (that is, the output of librsb-config  --static --ldflags, or a directory where the usual bin/include/lib subdirs with a regular RSB installation resides, or nothing to make the configure script invoke librsb-config)]), [if test x"$withval" = xno; then
want_rsb_libs= ; else if test x"$withval" = xyes ; then want_rsb_libs=yes ; else want_rsb_libs="$withval" ; fi ; fi], [want_rsb_libs=""])
if test x"$want_rsb_libs" != x ; then
	if test x"$want_rsb_libs" = xyes ; then
		want_rsb_libs="`librsb-config --static --ldflags`"
	else
		if test -d "$want_rsb_libs" ; then
			want_rsb_libs="`$want_rsb_libs/bin/librsb-config --static --ldflags`"
		else
			true;
			# we assume want_rsb_libs are linkage parameters
		fi
	fi
	FDEFINES="$FDEFINES $psblas_cv_define_prepend-DHAVE_LIBRSB"
fi
RSB_LIBS="$want_rsb_libs"
LIBS="$RSB_LIBS ${LIBS}"
dnl AC_CHECK_HEADERS([rsb.h], [ LIBS="${LIBS} $want_rsb_libs"], [])

PAC_CHECK_METIS
if test "x$psblas_cv_have_metis" == "xyes" ; then 
   FDEFINES="$psblas_cv_define_prepend-DHAVE_METIS $FDEFINES"
   CDEFINES="-DHAVE_METIS_ $psblas_cv_metis_includes $CDEFINES"
fi

PAC_CHECK_AMD
if test "x$psblas_cv_have_amd" == "xyes" ; then 
   FDEFINES="$psblas_cv_define_prepend-DHAVE_AMD $FDEFINES"
   CDEFINES="-DHAVE_AMD_ $psblas_cv_amd_includes $CDEFINES"
fi



# Note  : We would like to detect PSBLAS, but this is complicated by the
#         module symbols mangling rules, which are compiler specific !
#
#	  Moreover, the PSBLAS doesn't have an installer, currently.

###############################################################################
#		Library target directory and archive files.
###############################################################################

LIBDIR=lib
BASELIBNAME=libpsb_base.a
PRECLIBNAME=libpsb_prec.a
METHDLIBNAME=libpsb_krylov.a
UTILLIBNAME=libpsb_util.a

###############################################################################
# Variable substitutions : the Make.inc.in will have these @VARIABLES@
# substituted.

AC_SUBST(FDEFINES)
AC_SUBST(CDEFINES)
AC_SUBST(BASEMODNAME)
AC_SUBST(PRECMODNAME)
AC_SUBST(METHDMODNAME)
AC_SUBST(UTILMODNAME)
AC_SUBST(BASELIBNAME)
AC_SUBST(F90)
AC_SUBST(F90COPT)
AC_SUBST(MPF90)
AC_SUBST(MPF77)
AC_SUBST(MPCC)
AC_SUBST(FCOPT)
AC_SUBST(CCOPT)
AC_SUBST(EXTRA_OPT)
AC_SUBST(FAKEMPI)
AC_SUBST(FIFLAG)
AC_SUBST(FMFLAG)
AC_SUBST(MODEXT)
AC_SUBST(MPIF77)
AC_SUBST(FLINK)
AC_SUBST(F90LINK)
AC_SUBST(LIBS)
AC_SUBST(AR)
AC_SUBST(RANLIB)
AC_SUBST(MPIFC)
AC_SUBST(MPIFCC)

AC_SUBST(INSTALL)
AC_SUBST(INSTALL_DATA)
AC_SUBST(INSTALL_DIR)
AC_SUBST(INSTALL_LIBDIR)
AC_SUBST(INSTALL_INCLUDEDIR)
AC_SUBST(INSTALL_DOCSDIR)
AC_SUBST(INSTALL_SAMPLESDIR)

AC_SUBST(BLAS_LIBS)
AC_SUBST(AMD_LIBS)
AC_SUBST(METIS_LIBS)
AC_SUBST(LAPACK_LIBS)

AC_SUBST(CINCLUDES)
AC_SUBST(FINCLUDES)


PSBLASRULES='
PSBLDLIBS=$(LAPACK) $(BLAS) $(METIS_LIB) $(AMD_LIB) $(LIBS)
CDEFINES=$(PSBCDEFINES)
FDEFINES=$(PSBFDEFINES)
	

# These should be portable rules, arent they?
.c.o:
	$(CC) $(CCOPT) $(CINCLUDES) $(CDEFINES) -c $<
.f.o:
	$(FC) $(FCOPT) $(FINCLUDES)  -c $<
.f90.o:
	$(F90) $(F90COPT) $(FINCLUDES) -c $<
.F.o:
	$(FC)  $(FCOPT)   $(FINCLUDES) $(FDEFINES) -c $<
.F90.o:
	$(F90) $(F90COPT) $(FINCLUDES) $(FDEFINES) -c $<'


AC_SUBST(PSBLASRULES)
AC_SUBST(LIBDIR)
AC_SUBST(RSB_LIBS)
AC_SUBST(BASELIBNAME)
AC_SUBST(PRECLIBNAME)
AC_SUBST(METHDLIBNAME)
AC_SUBST(UTILLIBNAME)

###############################################################################
# the following files will be created by Automake

AC_CONFIG_FILES([Make.inc])
AC_OUTPUT()
#AC_OUTPUT(Make.inc Makefile)
###############################################################################

dnl Please note that brackets around variable identifiers are absolutely needed for compatibility..
AC_MSG_NOTICE([
	${PACKAGE_NAME} ${psblas_cv_version} has been configured as follows:

   	MPF90                 : ${MPF90}
	MPF77                 : ${MPF77}
	MPCC                  : ${MPICC}
dnl	F90LINK               : ${F90LINK}
	FLINK                 : ${FLINK}

	FDEFINES              : ${FDEFINES}
	CDEFINES              : ${CDEFINES}
dnl	CFLAGS                : ${CFLAGS}
dnl	FFLAGS                : ${FFLAGS}
dnl	FCFLAGS               : ${FCFLAGS}
	MODEXT                : ${MODEXT}
	FMFLAG                : ${FMFLAG}
   	F90COPT               : ${F90COPT}
   	FCOPT                 : ${FCOPT}
   	CCOPT                 : ${CCOPT}

dnl	ESSL/PESSL            : ${psblas_cv_have_essl} / ${psblas_cv_have_pessl}

	BLAS                  : ${BLAS_LIBS}
dnl	BLACS                 : ${BLACS_LIBS}

	METIS detected        : ${psblas_cv_have_metis}
	AMD detected          : ${psblas_cv_have_amd}
dnl	SuperLU detected      : ${psblas_cv_have_superlu}
dnl	SuperLU_Dist detected : ${psblas_cv_have_superludist}
dnl	UMFPack detected      : ${psblas_cv_have_umfpack}
dnl	F90                   : ${F90}
dnl	FC                    : ${FC}
dnl	CC                    : ${CC}
dnl	F77                   : ${F77}
dnl	FMFLAG                : ${FMFLAG}
dnl	FIFLAG                : ${FIFLAG}

dnl 	MPICXX                : ${MPICXX}
dnl	MPIF77                : ${MPIF77}
dnl	MPIFC                 : ${MPIFC}
dnl	MPIF90                : ${MPIF90}
dnl	MPCC                  : ${MPCC}
dnl	AR                    : ${AR}
dnl	RANLIB                : ${RANLIB}
        LIBS                  : ${LIBS}
dnl FIXME : CLIBS seems an useless variable..
dnl FIXME : FLIBS seems an useless variable..
dnl FIXME :  LIBS seems an useless variable..
dnl	CLIBS                 : ${CLIBS}
dnl	FLIBS                 : ${CLIBS}
dnl	 LIBS                 : ${LIBS}

dnl Note : we should use LDLIBS sooner or later!
	LDLIBS                : ${LDLIBS}

	LIBRARYPATH           : ${LIBRARYPATH}
	INCLUDEPATH           : ${INCLUDEPATH}
	MODULE_PATH            : ${MODULE_PATH}

	If you are satisfied, run 'make' to build ${PACKAGE_NAME} and its documentation; otherwise
	type  ./configure --help=short for a complete list of configure options specific to ${PACKAGE_NAME}.
dnl	To install the program and its documentation, run 'make install' if you are root,
dnl	or run 'su -c "make install"' if you are not root.
])

###############################################################################

