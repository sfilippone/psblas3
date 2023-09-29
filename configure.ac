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
AC_INIT([PSBLAS],3.8.1, [https://github.com/sfilippone/psblas3/issues])

# VERSION is the file containing the PSBLAS version code
# FIXME
psblas_cv_version="3.8.1"

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
        ./configure --with-libs=-L/some/directory/LIB  <- will append to LIBS
	       [ --with-spgpu=/path/to/spgpu] 
        FC=mpif90 CC=mpicc ./configure           <- will force FC,CC

	See ./configure --help=short fore more info.
--------------------------------------------------------------------------------
		])
###############################################################################
# Compilers detection: FC,CC should be set, if found.
###############################################################################
#
# Installation. 
#
#
AC_PROG_INSTALL
INSTALL="${INSTALL} -p"
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
INSTALL_MODULESDIR=$INSTALL_DIR/modules
AC_MSG_RESULT([$INSTALL_DIR $INSTALL_INCLUDEDIR $INSTALL_MODULESDIR $INSTALL_LIBDIR $INSTALL_DOCSDIR $INSTALL_SAMPLESDIR])

dnl
dnl We set our own FC flags, ignore those from AC_PROG_FC but not those from the
dnl environment variable. Same for C
dnl 
save_FCFLAGS="$FCFLAGS";
AC_PROG_FC([ftn xlf2003_r xlf2003 xlf95_r xlf95 xlf90 xlf pgf95 pgf90 ifort ifc  nagfor gfortran])
FCFLAGS="$save_FCFLAGS";
save_CFLAGS="$CFLAGS";
AC_PROG_CC([xlc pgcc icc gcc cc ])
if test "x$ac_cv_prog_cc_stdc" == "xno" ; then
   AC_MSG_ERROR([Problem : Need a C99 compiler ! ])
else
   C99OPT="$ac_cv_prog_cc_stdc";
fi	
CFLAGS="$save_CFLAGS";
AC_PROG_CXX([CC xlc++  icpc g++])

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
   MPICC="$CC";
   MPICXX="$CXX";	
else 
AC_LANG([C])
if test "X$MPICC" = "X" ; then
    # This is our MPICC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPICC],[mpxlc mpiicc mpcc mpicc cc])
fi
ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for C]])])

AC_LANG([Fortran])
AC_LANG([C++])
if test "X$MPICXX" = "X" ; then
    # This is our MPICC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPICXX],[mpxlc++ mpiicpc  mpicxx])
fi
ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for C++]])])
AC_LANG([Fortran])

if test "X$MPIFC" = "X" ; then
    # This is our MPIFC compiler preference: it will override ACX_MPI's first try.
    AC_CHECK_PROGS([MPIFC],[mpxlf2003_r mpxlf2003 mpxlf95_r mpxlf90 mpiifort mpf95 mpf90 mpifort mpif95 mpif90 ftn ])
fi

ACX_MPI([], [AC_MSG_ERROR([[Cannot find any suitable MPI implementation for Fortran]])])

FC="$MPIFC" ;
CC="$MPICC";
CXX="$MPICXX";
fi

AC_LANG([C])

dnl Now on, MPIFC should be set, and MPICC

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
PAC_ARG_WITH_FLAGS(cxxopt,CXXOPT)
PAC_ARG_WITH_FLAGS(fcopt,FCOPT)
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

AC_PROG_RANLIB

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
 else
   psblas_cv_fc=""
  # unsupported MPI Fortran compiler
     AC_MSG_NOTICE([[Unknown Fortran compiler, proceeding with fingers crossed !]])
 fi
fi
if test "X$psblas_cv_fc" == "Xgcc" ; then
PAC_HAVE_MODERN_GFORTRAN(
	[],
        [AC_MSG_ERROR([Bailing out.])]
)
fi


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
# Honor CFLAGS if they were specified explicitly, but --with-ccopt take precedence
if test "X$CCOPT" == "X" ; then
  CCOPT="$CFLAGS";
fi	
if test "X$CCOPT" == "X" ; then
 if test "X$psblas_cv_fc" == "Xgcc" ; then
 # note that no space should be placed around the equality symbol in assignements
 # Note : 'native' is valid _only_ on GCC/x86  (32/64 bits)
   CCOPT="-g -O3  $CCOPT"

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
     CCOPT="-g -O2 $CCOPT"
 fi
fi
#CFLAGS="${CCOPT}"
if test "X$CXXOPT" == "X" ; then
  CXXOPT="$CXXFLAGS";
fi	
if test "X$CXXOPT" == "X" ; then
 if test "X$psblas_cv_fc" == "Xgcc" ; then
 # note that no space should be placed around the equality symbol in assignements
 # Note : 'native' is valid _only_ on GCC/x86  (32/64 bits)
   CXXOPT="-g -O3  $CXXOPT"

 elif test "X$psblas_cv_fc" == X"xlf" ; then
   # XL compiler : consider using -qarch=auto
     CXXOPT="-O3 -qarch=auto $CXXOPT"
 elif test "X$psblas_cv_fc" == X"ifc" ; then
   # other compilers ..
     CXXOPT="-O3 $CXXOPT"
 elif test "X$psblas_cv_fc" == X"pg" ; then
   # other compilers ..
     CXXCOPT="-fast $CXXOPT"
  # NOTE : PG & Sun use -fast instead -O3
 elif test "X$psblas_cv_fc" == X"sun" ; then
   # other compilers ..
     CXXOPT="-fast $CXXOPT"
 elif test "X$psblas_cv_fc" == X"cray" ; then
     CXXOPT="-O3 $CXXOPT"
     MPICXX="CC"
 else 
     CXXOPT="-g -O3 $CXXOPT"
 fi
fi

# Honor FCFLAGS if they were specified explicitly, but --with-fcopt take precedence
if test "X$FCOPT" == "X" ; then
  FCOPT="$FCFLAGS";
fi
   

if test "X$FCOPT" == "X" ; then
 if test "X$psblas_cv_fc" == "Xgcc" ; then
 # note that no space should be placed around the equality symbol in assignations
 # Note : 'native' is valid _only_ on GCC/x86  (32/64 bits)
   FCOPT="-g -O3 -frecursive $FCOPT"  
 elif test "X$psblas_cv_fc" == X"xlf" ; then
   # XL compiler : consider using -qarch=auto
     FCOPT="-O3 -qarch=auto  -qlanglvl=extended -qxlf2003=polymorphic:autorealloc $FCOPT"
     FCFLAGS="-qhalt=e -qlanglvl=extended -qxlf2003=polymorphic:autorealloc $FCFLAGS"
 elif test "X$psblas_cv_fc" == X"ifc" ; then
   # other compilers ..
     FCOPT="-O3 -recursive $FCOPT"
 elif test "X$psblas_cv_fc" == X"pg" ; then
   # other compilers ..
     FCOPT="-fast $FCOPT"
  # NOTE : PG & Sun use -fast instead -O3
 elif test "X$psblas_cv_fc" == X"sun" ; then
   # other compilers ..
     FCOPT="-fast $FCOPT"
 elif test "X$psblas_cv_fc" == X"cray" ; then
     FCOPT="-O3 -em -J. $FCOPT"
 elif test "X$psblas_cv_fc" == X"nag" ; then
   # NAG  compiler ..
     FCOPT="-O2 "
  # NOTE : PG & Sun use -fast instead -O3
 else 
     FCOPT="-g -O2 $FCOPT"
 fi
fi
if test "X$psblas_cv_fc" == X"nag" ; then
   # Add needed options 
     FCOPT="$FCOPT -dcfuns -f2003 -wmismatch=mpi_scatterv,mpi_alltoallv,mpi_gatherv,mpi_allgatherv"
     EXTRA_OPT="-mismatch_all"
fi
if test "X$psblas_cv_fc" == "Xgcc" ; then
   FCOPT="-frecursive $FCOPT"  
   PAC_HAVE_GFORTRAN_10(
	[FCOPT="-fallow-argument-mismatch $FCOPT"],
        [])  
fi   


# COPT,FCOPT are aliases for CFLAGS,FCFLAGS .
##############################################################################
# Compilers variables selection 
##############################################################################
FC=${FC}
CC=${CC}
CXX=${CXX}
CCOPT="$CCOPT $C99OPT"


##############################################################################
# Choice of our compilers, needed by Make.inc
##############################################################################
if test "X$FLINK" == "X" ; then
   FLINK=${MPF90}
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

PAC_ARG_WITH_IPK
PAC_ARG_WITH_LPK
# Defaults for IPK/LPK
if test x"$pac_cv_ipk_size" == x"" ; then
  pac_cv_ipk_size=4
fi
if test x"$pac_cv_lpk_size" == x"" ; then
  pac_cv_lpk_size=8
fi
# Enforce sensible combination
if (( $pac_cv_lpk_size < $pac_cv_ipk_size )); then
   AC_MSG_NOTICE([[Invalid combination of size specs IPK ${pac_cv_ipk_size} LPK ${pac_cv_lpk_size}. ]]);
   AC_MSG_NOTICE([[Forcing equal values]])	
   pac_cv_lpk_size=$pac_cv_ipk_size;	
fi				
FDEFINES="$psblas_cv_define_prepend-DIPK${pac_cv_ipk_size} $FDEFINES";
FDEFINES="$psblas_cv_define_prepend-DLPK${pac_cv_lpk_size} $FDEFINES";	
CDEFINES="-DIPK${pac_cv_ipk_size} -DLPK${pac_cv_lpk_size} $CDEFINES"

FLINK="$MPIFC"
PAC_ARG_OPENMP()
if test x"$pac_cv_openmp" == x"yes" ; then
   FDEFINES="$psblas_cv_define_prepend-DOPENMP $FDEFINES";
   CDEFINES="-DOPENMP $CDEFINES";	
   FCOPT="$FCOPT $pac_cv_openmp_fcopt";
   CCOPT="$CCOPT $pac_cv_openmp_ccopt";
   CXXOPT="$CXXOPT $pac_cv_openmp_cxxopt";	
   FLINK="$FLINK $pac_cv_openmp_fcopt";
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
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_TEST_EXTENDS( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for EXTENDS. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_TEST_CLASS_TBP( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for CLASS and type bound procedures. 
 Please get a Fortran compiler that supports them, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_TEST_SOURCE( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for SOURCE= allocation. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_HAVE_MOVE_ALLOC(
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for MOVE_ALLOC. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_TEST_ISO_C_BIND( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for ISO_C_BINDING. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_TEST_SAME_TYPE( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for SAME_TYPE_AS. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_TEST_EXTENDS_TYPE( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for EXTENDS_TYPE_OF. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_TEST_MOLD( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for MOLD= allocation. 
 Please get a Fortran compiler that supports it, e.g. GNU Fortran 4.8.])]
)

PAC_FORTRAN_TEST_VOLATILE(
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for	VOLATILE])]	
)

PAC_FORTRAN_TEST_ISO_FORTRAN_ENV(
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for	ISO_FORTRAN_ENV])]	
)

PAC_FORTRAN_TEST_FINAL( 
	[],
        [AC_MSG_ERROR([Sorry, cannot build PSBLAS without support for FINAL])]	
)


#
# Optional features
#


PAC_FORTRAN_TEST_GENERICS( 
	[],
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_BUGGY_GENERICS $FDEFINES"]
)

PAC_FORTRAN_TEST_FLUSH( 
	[FDEFINES="$psblas_cv_define_prepend-DHAVE_FLUSH_STMT $FDEFINES"],
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

AC_ARG_WITH(rsb, AS_HELP_STRING([--with-rsb], [Specify Recursive Sparse BLAS library linkage info (that is, the output of librsb-config  --static --ldflags, or a directory where the usual bin/include/lib subdirs with a regular RSB installation resides, or nothing to make the configure script invoke librsb-config)]), [if test x"$withval" = xno; then
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

AC_MSG_CHECKING([Compatibility between metis and LPK])
if  test  "x$pac_cv_lpk_size" == "x4" ; then
  if test  "x$pac_cv_metis_idx" == "x64" ; then
  dnl mismatch between metis size and PSBLAS LPK
       psblas_cv_have_metis="no";
 dnl 
   fi
 fi			   
 if test "x$pac_cv_lpk_size" == "x8" ; then
   if test "x$pac_cv_metis_idx" == "x32" ; then 
   dnl mismatch between metis size and PSBLAS LPK
      psblas_cv_have_metis="no";
   fi
 fi			   
AC_MSG_RESULT([$psblas_cv_have_metis])

if test "x$pac_cv_metis_idx" == "xunknown" ; then
  dnl mismatch between metis size and PSBLAS LPK
     AC_MSG_NOTICE([Unknown  METIS bitsize.])	
     $psblas_cv_have_metis = "no";
fi		   
if test "x$pac_cv_metis_real" == "xunknown" ; then
  dnl mismatch between metis size and PSBLAS LPK
     AC_MSG_NOTICE([Unknown  METIS REAL bitsize.])	
     $psblas_cv_have_metis = "no";
fi		   
if test "x$psblas_cv_have_metis" == "xyes" ; then 
   FDEFINES="$psblas_cv_define_prepend-DHAVE_METIS $psblas_cv_define_prepend-DMETIS_$pac_cv_metis_idx $psblas_cv_define_prepend-DMETIS_REAL_$pac_cv_metis_real $FDEFINES"
   CDEFINES="-DHAVE_METIS_ $psblas_cv_metis_includes $CDEFINES -DMETIS_$pac_cv_metis_idx -DMETIS_REAL_$pac_cv_metis_real"
   METISINCFILE=$psblas_cv_metisincfile
fi

PAC_CHECK_AMD
if test "x$psblas_cv_have_amd" == "xyes" ; then 
   FDEFINES="$psblas_cv_define_prepend-DHAVE_AMD $FDEFINES"
   CDEFINES="-DHAVE_AMD_ $psblas_cv_amd_includes $CDEFINES"
fi



PAC_CHECK_CUDA()

if test "x$pac_cv_have_cuda" == "xyes"; then

PAC_CHECK_CUDA_VERSION()
CUDA_VERSION="$pac_cv_cuda_version";
CUDA_SHORT_VERSION=$(expr $pac_cv_cuda_version / 1000);
PAC_CHECK_SPGPU()
if test "x$pac_cv_have_spgpu" == "xyes" ; then 
   GPUD=gpud;
   GPULD=gpuld;	
   EXTRALDLIBS="-lstdc++";
fi
AC_MSG_NOTICE([At this point GPUTARGET is $GPUD $GPULD])

PAC_ARG_WITH_CUDACC()
if test "x$pac_cv_cudacc" == "x"; then
   pac_cv_cudacc="30,35,37,50,60";
fi
CUDEFINES="";
for cc in `echo $pac_cv_cudacc|sed 's/,/ /gi'`
do
  CUDEFINES="$CUDEFINES -gencode arch=compute_$cc,code=sm_$cc";
done
if test "x$pac_cv_cuda_version" != "xunknown"; then
   CUDEFINES="$CUDEFINES -DCUDA_SHORT_VERSION=${CUDA_SHORT_VERSION} -DCUDA_VERSION=${CUDA_VERSION}"
   FDEFINES="$FDEFINES -DCUDA_SHORT_VERSION=${CUDA_SHORT_VERSION} -DCUDA_VERSION=${CUDA_VERSION}"
   CDEFINES="$CDEFINES -DCUDA_SHORT_VERSION=${CUDA_SHORT_VERSION} -DCUDA_VERSION=${CUDA_VERSION}"    
fi

fi

if test "x$pac_cv_ipk_size" != "x4"; then
  AC_MSG_NOTICE([For CUDA I need psb_ipk_ to be 4 bytes but it is $pac_cv_ipk_size,  disabling CUDA/SPGPU])
  GPUD="";
  GPULD="";
  CUDEFINES="";
  CUDA_INCLUDES="";
  CUDA_LIBS="";
fi




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
AC_SUBST(MPIFC)
AC_SUBST(MPICC)
AC_SUBST(FCOPT)
AC_SUBST(CCOPT)
AC_SUBST(CXXOPT)
AC_SUBST(EXTRA_OPT)
AC_SUBST(FAKEMPI)
AC_SUBST(FIFLAG)
AC_SUBST(FMFLAG)
AC_SUBST(MODEXT)
AC_SUBST(FLINK)
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
AC_SUBST(INSTALL_MODULESDIR)

AC_SUBST(BLAS_LIBS)
AC_SUBST(AMD_LIBS)
AC_SUBST(METIS_LIBS)
AC_SUBST(LAPACK_LIBS)

AC_SUBST(CINCLUDES)
AC_SUBST(FINCLUDES)

PSBLASRULES='
PSBLDLIBS=$(LAPACK) $(BLAS) $(METIS_LIB) $(AMD_LIB) $(LIBS)
CXXDEFINES=$(PSBCXXDEFINES)
CDEFINES=$(PSBCDEFINES)
FDEFINES=$(PSBFDEFINES)
	

# These should be portable rules, arent they?
.c.o:
	$(CC) $(CCOPT) $(CINCLUDES) $(CDEFINES) -c $< -o $@
.f90.o:
	$(FC) $(FCOPT) $(FINCLUDES) -c $< -o $@
.F90.o:
	$(FC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $< -o $@
.cpp.o:
	$(CXX) $(CXXOPT) $(CXXINCLUDES) $(CXXDEFINES) -c $< -o $@'


AC_SUBST(PSBLASRULES)
AC_SUBST(LIBDIR)
AC_SUBST(RSB_LIBS)
AC_SUBST(BASELIBNAME)
AC_SUBST(PRECLIBNAME)
AC_SUBST(METHDLIBNAME)
AC_SUBST(UTILLIBNAME)
AC_SUBST(METISINCFILE)
AC_SUBST(SPGPU_FLAGS)
AC_SUBST(SPGPU_LIBS)
AC_SUBST(SPGPU_DIR)
AC_SUBST(SPGPU_DEFINES)
AC_SUBST(SPGPU_INCLUDES)
AC_SUBST(SPGPU_INCDIR)
AC_SUBST(EXTRALDLIBS)
AC_SUBST(CUDA_DIR)
AC_SUBST(CUDA_DEFINES)
AC_SUBST(CUDA_INCLUDES)
AC_SUBST(CUDA_LIBS)
AC_SUBST(CUDA_VERSION)
AC_SUBST(CUDA_SHORT_VERSION)
AC_SUBST(CUDA_NVCC)
AC_SUBST(CUDEFINES)
AC_SUBST(GPUD)
AC_SUBST(GPULD)
###############################################################################
# the following files will be created by Automake

AC_CONFIG_FILES([Make.inc])
AC_CONFIG_FILES([util/psb_metis_int.h])	
AC_OUTPUT()
#AC_OUTPUT(Make.inc Makefile)
###############################################################################

dnl Please note that brackets around variable identifiers are absolutely needed for compatibility..
AC_MSG_NOTICE([
	${PACKAGE_NAME} ${psblas_cv_version} has been configured as follows:

   	MPIFC		      : ${MPIFC}
	MPICC                 : ${MPICC}
	FLINK                 : ${FLINK}

	FDEFINES              : ${FDEFINES}
	CDEFINES              : ${CDEFINES}
	MODEXT                : ${MODEXT}
	FMFLAG                : ${FMFLAG}
   	FCOPT                 : ${FCOPT}
   	CCOPT                 : ${CCOPT}

	BLAS                  : ${BLAS_LIBS}

	METIS usable          : ${psblas_cv_have_metis}
	METIS bitsize         : ${pac_cv_metis_idx}
	METIS realsize        : ${pac_cv_metis_real}
	AMD detected          : ${psblas_cv_have_amd}
        LIBS                  : ${LIBS}
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

