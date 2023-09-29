dnl
dnl $Id$
dnl
dnl 20080206
dnl M4 macros for the PSBLAS library and useful for packages using PSBLAS.
dnl

dnl @synopsis PAC_CHECK_LIBS
dnl
dnl Tries to detect the presence of a specific function among various libraries, using AC_CHECK_LIB
dnl repeatedly on the specified libraries.
dnl 
dnl Example use:
dnl
dnl PAC_CHECK_LIBS([atlas blas],
dnl		[dgemm],
dnl		[have_dgemm=yes],
dnl		[have_dgemm=no])
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
dnl 20080211 modified slighty from original.
AC_DEFUN([PAC_CHECK_LIBS],
[
 pac_check_libs_ok=no
 [for pac_check_libs_f in $2 
 do ]
 [for pac_check_libs_l in $1 
 do ]
    if test x"$pac_check_libs_ok" == xno ; then
     AC_CHECK_LIB([$pac_check_libs_l],[$pac_check_libs_f], [pac_check_libs_ok=yes; pac_check_libs_LIBS="-l$pac_check_libs_l"],[],[$5])
    fi
  done
  done
 # Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
 [ if test x"$pac_check_libs_ok" = xyes ; then
	$3
 else
        pac_check_libs_ok=no
        $4
 fi
 ]
])dnl 

dnl # How do we determine netlib blacs ? Here it is :
dnl AC_CHECK_LIB(blacsCinit_MPI-LINUX-0,BI_Asend,
dnl AC_CHECK_LIB(blacs_MPI-LINUX-0,BI_Asend,
dnl                         [acx_cv_blacs_ok=yes; BLACS_LIBS="-lblacsCinit_MPI-LINUX-0"])
dnl AC_F77_FUNC([BI_Asend])
dnl AC_CHECK_LIB(blacsCinit_MPI-LINUX-0,[$BI_Asend],
dnl                        [acx_cv_blacs_ok=yes; BLACS_LIBS="-lblacsCinit_MPI-LINUX-0"])
dnl AC_FC_FUNC([BI_Iam])
dnl AC_CHECK_LIB(blacsCinit_MPI-LINUX-0,[$BI_Asend],
dnl                         [acx_cv_blacs_ok=yes; BLACS_LIBS="-lblacsCinit_MPI-LINUX-0"])


dnl @synopsis PAC_FORTRAN_HAVE_MOVE_ALLOC( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program with move_alloc (a Fortran 2003 function).
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN([PAC_FORTRAN_HAVE_MOVE_ALLOC],
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran MOVE_ALLOC intrinsic])
 AC_LANG_PUSH([Fortran])
 ac_ext='f90';
 AC_COMPILE_IFELSE([ program test_move_alloc
		       integer, allocatable :: a(:), b(:)
		       allocate(a(3))
		       call move_alloc(a, b)
		       print *, allocated(a), allocated(b)
		       print *, b
		     end program test_move_alloc],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])



dnl @synopsis PAC_CHECK_HAVE_CRAYFTN( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will check if MPIFC is $FC.
dnl The check will proceed by compiling a small Fortran program
dnl containing the _CRAYFTN macro, which should be defined in the
dnl gfortran compiled programs.
dnl
dnl On pass, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_HAVE_CRAYFTN,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([for Cray Fortran])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program main
#ifdef _CRAYFTN 
              print *, "Cray FTN!"
#else
        this program will fail
#endif
           end],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])



dnl @synopsis PAC_CHECK_HAVE_GFORTRAN( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will check if MPIFC is $FC.
dnl The check will proceed by compiling a small Fortran program
dnl containing the __GNUC__ macro, which should be defined in the
dnl gfortran compiled programs.
dnl
dnl On pass, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_HAVE_GFORTRAN,
[AC_MSG_CHECKING([for GNU Fortran])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program main
#ifdef __GNUC__ 
              print *, "GCC!"
#else
        this program will fail
#endif
           end],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_HAVE_MODERN_GFORTRAN( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will check if the GNU fortran version is suitable for PSBLAS.
dnl If yes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl Note : Will use MPIFC; if unset, will use '$FC'.
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN(PAC_HAVE_MODERN_GFORTRAN,
 [AC_MSG_CHECKING([for recent GNU Fortran])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program main
#if ( __GNUC__ >= 4 && __GNUC_MINOR__ >= 9 ) || ( __GNUC__ > 4 )
              print *, "ok"
#else
        this program will fail
#endif
           end],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])
		     AC_MSG_NOTICE([Sorry, we require GNU Fortran version 4.9 or later.])
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_HAVE_GFORTRAN_10( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will check if the GNU fortran version is suitable for PSBLAS.
dnl If yes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl Note : Will use MPIFC; if unset, will use '$FC'.
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN(PAC_HAVE_GFORTRAN_10,
 [AC_MSG_CHECKING([for version 10 or later of GNU Fortran])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program main
#if ( __GNUC__ >= 10  ) 
              print *, "ok"
#else
        this program will fail
#endif
           end],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_CHECK_HAVE_MPI_MOD( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will determine if the fortran compiler MPIFC needs to include mpi.h or needs
dnl to use the mpi module.
dnl
dnl If yes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl Modified Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_CHECK_HAVE_MPI_MOD,
 [AC_MSG_CHECKING([for Fortran MPI mod])
  AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program test
             use mpi
           end program test],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])


dnl @synopsis PAC_FORTRAN_CHECK_HAVE_MPI_MOD_F08( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will determine if the fortran compiler MPIFC provides mpi_f08
dnl
dnl If yes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl 
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl Modified Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_CHECK_HAVE_MPI_MOD_F08,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([MPI Fortran 2008 interface])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program test
             use mpi_f08
           end program test],
		   [  AC_MSG_RESULT([yes])
		      pac_cv_mpi_f08="yes";
		      ifelse([$1], , :, [ $1])],
		   [  AC_MSG_RESULT([no])
	              pac_cv_mpi_f08="no";
		      echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		      cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		      ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])



dnl @synopsis PAC_ARG_WITH_FLAGS(lcase_name, UCASE_NAME)
dnl
dnl Test for --with-lcase_name="compiler/loader flags".  if defined, prepend 
dnl flags to standard UCASE_NAME definition.
dnl
dnl Use this macro to facilitate additional special flags that should be
dnl passed on to the preprocessor/compilers/loader.
dnl
dnl NOTE : Renamed after TAC_ARG_WITH_FLAGS as in the Trilinos-8.0.4 package.
dnl 
dnl NOTE : This macro works in a way the user should invoke
dnl         --with-flags=...
dnl	   only once, otherwise the first one will take effect.
dnl
dnl Example use:
dnl 
dnl PAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
dnl 
dnl tests for --with-cxxflags and pre-pends to CXXFLAGS
dnl 
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl @notes  Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_WITH_FLAGS],
[
AC_MSG_CHECKING([whether additional [$2] flags should be added (should be invoked only once)])
dnl AC_MSG_CHECKING([whether additional [$2] flags should be added])
AC_ARG_WITH($1,
AS_HELP_STRING([--with-$1], 
[additional [$2] flags to be added: will prepend to [$2]]),
[
$2="${withval} ${$2}"
AC_MSG_RESULT([$2 = ${$2}])
],
AC_MSG_RESULT(no)
)
])


dnl @synopsis PAC_ARG_WITH_LIBS
dnl
dnl Test for --with-libs="name(s)".
dnl 
dnl Prepends the specified name(s) to the list of libraries to link 
dnl with.  
dnl
dnl note: Renamed after PAC_ARG_WITH_LIBS as in the Trilinos package.
dnl
dnl Example use:
dnl
dnl PAC_ARG_WITH_LIBS
dnl 
dnl tests for --with-libs and pre-pends to LIBS
dnl
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([PAC_ARG_WITH_LIBS],
[
AC_MSG_CHECKING([whether additional libraries are needed])
AC_ARG_WITH(libs,
AS_HELP_STRING([--with-libs], 
[List additional link flags  here.  For example, --with-libs=-lspecial_system_lib
or --with-libs=-L/path/to/libs]),
[
LIBS="${withval} ${LIBS}"
AC_MSG_RESULT([LIBS = ${LIBS}])
],
AC_MSG_RESULT(no)
)
]
)


dnl @synopsis PAC_ARG_SERIAL_MPI
dnl
dnl Test for --enable-serial
dnl 
dnl 
dnl
dnl Example use:
dnl
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_SERIAL_MPI],
[AC_MSG_CHECKING([whether we want serial  mpi stubs])
AC_ARG_ENABLE(serial,
AS_HELP_STRING([--enable-serial], 
[Specify whether to enable a fake mpi library to run in serial mode. ]),
[
pac_cv_serial_mpi="yes";
]
dnl ,
dnl [pac_cv_serial_mpi="no";]
)
if test x"$pac_cv_serial_mpi" == x"yes" ; then
   AC_MSG_RESULT([yes.])
else
 pac_cv_serial_mpi="no";
 AC_MSG_RESULT([no.])
fi
]
)

dnl @synopsis PAC_ARG_OPENMP
dnl
dnl Test for --enable-openmp
dnl 
dnl 
dnl
dnl Example use:
dnl
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_OPENMP],
[AC_MSG_CHECKING([whether we want openmp ])
AC_ARG_ENABLE(openmp,
AS_HELP_STRING([--enable-openmp], 
[Specify whether to enable openmp. ]),
[
pac_cv_openmp="yes";
]
dnl ,
dnl [pac_cv_openmp="no";]
	     )
if test x"$pac_cv_openmp" == x"yes" ; then
   AC_MSG_RESULT([yes.])
   AC_LANG_PUSH([Fortran])
   AC_OPENMP() 
   pac_cv_openmp_fcopt="$OPENMP_FCFLAGS";
   AC_LANG_POP()
   AC_LANG_PUSH([C])
   AC_OPENMP() 
   pac_cv_openmp_ccopt="$OPENMP_CFLAGS";
   AC_LANG_POP()
   AC_LANG_PUSH([C++])
   AC_OPENMP() 
   pac_cv_openmp_cxxopt="$OPENMP_CXXFLAGS";
   AC_LANG_POP()
else
 pac_cv_openmp="no";
 AC_MSG_RESULT([no.])
fi
]
)

dnl @synopsis PAC_ARG_LONG_INTEGERS
dnl
dnl Test for --enable-long-integers
dnl 
dnl 
dnl
dnl Example use:
dnl
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_LONG_INTEGERS],
[
AC_MSG_CHECKING([whether we want long (8 bytes) integers])
AC_ARG_ENABLE(long-integers,
AS_HELP_STRING([--enable-long-integers], 
[Specify usage of 64 bits integers. ]),
[
pac_cv_long_integers="yes";
]
dnl ,
dnl [pac_cv_long_integers="no";]
)
if test x"$pac_cv_long_integers" == x"yes" ; then
   AC_MSG_RESULT([yes.])
else
 pac_cv_long_integers="no";
 AC_MSG_RESULT([no.])
fi
]
)

dnl @synopsis PAC_ARG_WITH_IPK
dnl
dnl Test for --with-ipk 
dnl 
dnl 
dnl
dnl Example use: --with-ipk=4
dnl
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_WITH_IPK],
[
AC_MSG_CHECKING([what size in bytes we want for local indices and data])
AC_ARG_WITH(ipk,
	    AS_HELP_STRING([--with-ipk=<bytes>], 
			   [Specify the size in bytes for local indices and data, default 4 bytes. ]),
	    [pac_cv_ipk_size=$withval;],
	    [pac_cv_ipk_size=4;]
	   )
if test x"$pac_cv_ipk_size" == x"4"  || test  x"$pac_cv_ipk_size" == x"8" ; then
   AC_MSG_RESULT([Size: $pac_cv_ipk_size.])
else
  AC_MSG_RESULT([Unsupported value for IPK: $pac_cv_ipk_size, defaulting to 4.])
  pac_cv_ipk_size=4;
fi
]
)

dnl @synopsis PAC_ARG_WITH_LPK
dnl
dnl Test for --with-lpk 
dnl 
dnl 
dnl
dnl Example use: --with-lpk=8
dnl
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_WITH_LPK],
[
 AC_MSG_CHECKING([what size in bytes we want for global indices and data])
 AC_ARG_WITH(lpk,
	     AS_HELP_STRING([--with-lpk=<bytes>], 
			    [Specify the size in bytes for global indices and data, default 8 bytes. ]),
	     [pac_cv_lpk_size=$withval;],
	     [pac_cv_lpk_size=8;]
	    )
if test x"$pac_cv_lpk_size" == x"4" || test x"$pac_cv_lpk_size" == x"8"; then
  AC_MSG_RESULT([Size: $pac_cv_lpk_size.])
else
  AC_MSG_RESULT([Unsupported value for LPK: $pac_cv_lpk_size, defaulting to 8.])
  pac_cv_lpk_size=8;
fi
]
)


dnl @synopsis PAC_FORTRAN_HAVE_PSBLAS( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program using the PSBLAS library
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_HAVE_PSBLAS,
ac_exeext=''
ac_ext='f90'
ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FCFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
dnl Warning : square brackets are EVIL!
[cat > conftest.$ac_ext <<EOF
           program test
	       use psb_base_mod
           end program test
EOF
if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
  ifelse([$1], , :, [rm -rf conftest*
  $1])
else
  echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
  cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
ifelse([$2], , , [  rm -rf conftest*
  $2
])dnl
fi
rm -f conftest*])

dnl @synopsis PAC_FORTRAN_TEST_TR15581( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the TR15581 Fortran extension support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_TR15581,
[AC_MSG_CHECKING([support for Fortran allocatables TR15581])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
module conftest
  type outer
    integer,  allocatable :: v(:)
  end type outer

  interface foo
    module procedure foov, food
  end interface
contains

  subroutine foov(a,b)

    implicit none
    integer, allocatable, intent(inout) :: a(:)
    integer, allocatable, intent(out) :: b(:)


    allocate(b(size(a)))

  end subroutine foov
  subroutine food(a,b)

    implicit none
    type(outer), intent(inout) :: a
    type(outer), intent(out) :: b


    allocate(b%v(size(a%v)))

  end subroutine food

end module conftest



program testtr15581
  use conftest
  type(outer) :: da, db
  integer, allocatable :: a(:), b(:)

  allocate(a(10),da%v(10))
  a = (/ (i,i=1,10) /)
  da%v = (/ (i,i=1,10) /)
  call foo(a,b)
  call foo(da,db)
  write(*,*) b
  write(*,*) db%v

end program testtr15581],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_TEST_VOLATILE( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the VOLATILE Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_VOLATILE,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran VOLATILE])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program conftest
  integer, volatile :: i, j
end program conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])


dnl @synopsis PAC_FORTRAN_TEST_GENERICS( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile a program checking the GENERIC Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_GENERICS,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([test GENERIC interfaces])
AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
module conftest

  interface foo 
    subroutine i_sub_foo(v)
      integer, intent(inout) :: v(:)
    end subroutine i_sub_foo
  end interface foo

  interface bar
    procedure i_sub_foo
  end interface bar

end module conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])


dnl @synopsis PAC_FORTRAN_TEST_EXTENDS( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the EXTENDS Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_EXTENDS,
ac_exeext=''
ac_ext='f90'
ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FCFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran EXTENDS])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program conftest
  type foo
    integer :: i
  end type foo
  type, extends(foo) :: bar
    integer j
  end type bar 
  type(bar) :: barvar
end program conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_TEST_CLASS_TBP( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the TBP Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_CLASS_TBP,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran CLASS TBP])
AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='F90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
module conftest_mod
  type foo
    integer :: i 
  contains
    procedure, pass(a) :: doit
    procedure, pass(a) :: getit
  end type foo

  private doit,getit
contains
  subroutine  doit(a) 
    class(foo) :: a
    
    a%i = 1
    write(*,*) 'FOO%DOIT base version'
  end subroutine doit
  function getit(a) result(res)
    class(foo) :: a
    integer :: res

    res = a%i
  end function getit

end module conftest_mod
program conftest
  use conftest_mod
  type(foo) :: foovar
end program conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_TEST_FINAL( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the FINAL Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_FINAL,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran FINAL])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
module conftest_mod
  type foo
    integer :: i 
  contains
    final  :: destroy_foo
  end type foo

  private destroy_foo
contains
  subroutine destroy_foo(a)
    type(foo) :: a
     ! Just a test
  end subroutine destroy_foo
end module conftest_mod
program conftest
  use conftest_mod
  type(foo) :: foovar
end program conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_TEST_SAME_TYPE( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the SAME_TYPE_AS Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_SAME_TYPE,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran SAME_TYPE_AS])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program stt
  type foo
    integer :: i
  end type foo
  type, extends(foo) :: new_foo
    integer :: j
  end type new_foo
  type(foo) :: foov
  type(new_foo) :: nfv1, nfv2

    
  write(*,*) 'foov == nfv1? ', same_type_as(foov,nfv1)
  write(*,*) 'nfv2 == nfv1? ', same_type_as(nfv2,nfv1)
end program stt],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_TEST_EXTENDS_TYPE( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the EXTENDS_TYPE_OF Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_EXTENDS_TYPE,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran EXTENDS_TYPE_OF])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program xtt
  type foo
    integer :: i
  end type foo
  type, extends(foo) :: new_foo
    integer :: j
  end type new_foo
  type(foo) :: foov
  type(new_foo) :: nfv1, nfv2

  write(*,*) 'nfv1 extends foov? ', extends_type_of(nfv1,foov)
end program xtt],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_CHECK_BLACS
dnl
dnl Will try to find the BLACS
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_BLACS,
[AC_ARG_WITH(blacs, AS_HELP_STRING([--with-blacs=LIB], [Specify BLACSLIBNAME or -lBLACSLIBNAME or the absolute library filename.]),
        [psblas_cv_blacs=$withval],
        [psblas_cv_blacs=''])

case $psblas_cv_blacs in
	yes | "") ;;
	-* | */* | *.a | *.so | *.so.* | *.o) 
	     BLACS_LIBS="$psblas_cv_blacs" ;;
	*) BLACS_LIBS="-l$psblas_cv_blacs" ;;
esac

#
# Test user-defined BLACS
#
if test x"$psblas_cv_blacs" != "x" ; then
      save_LIBS="$LIBS";
      AC_LANG([Fortran])
      LIBS="$BLACS_LIBS $LIBS"
      AC_MSG_CHECKING([for dgesd2d in $BLACS_LIBS])
      AC_TRY_LINK_FUNC(dgesd2d, [psblas_cv_blacs_ok=yes], [psblas_cv_blacs_ok=no;BLACS_LIBS=""])
      AC_MSG_RESULT($psblas_cv_blacs_ok)

     if test x"$psblas_cv_blacs_ok" == x"yes";  then 
     AC_MSG_CHECKING([for blacs_pinfo in $BLACS_LIBS])
     AC_TRY_LINK_FUNC(blacs_pinfo, [psblas_cv_blacs_ok=yes], [psblas_cv_blacs_ok=no;BLACS_LIBS=""])
     AC_MSG_RESULT($psblas_cv_blacs_ok)
     fi 
     LIBS="$save_LIBS";
fi
AC_LANG([C])	

######################################
# System BLACS with PESSL default names. 
######################################
if test x"$BLACS_LIBS" == "x" ; then
   AC_LANG([Fortran])
   PAC_CHECK_LIBS([blacssmp blacsp2 blacs], 
	[dgesd2d],
	[psblas_cv_blacs_ok=yes; LIBS="$LIBS $pac_check_libs_LIBS "  ]
	[BLACS_LIBS="$pac_check_libs_LIBS" ]
	AC_MSG_NOTICE([BLACS libraries detected.]),[]
    )
    if test x"$BLACS_LIBS" != "x"; then 
          save_LIBS="$LIBS";
          LIBS="$BLACS_LIBS $LIBS"
          AC_MSG_CHECKING([for blacs_pinfo in $BLACS_LIBS])
          AC_LANG([Fortran])
	  AC_TRY_LINK_FUNC(blacs_pinfo, [psblas_cv_blacs_ok=yes], [psblas_cv_blacs_ok=no;BLACS_LIBS=""])
          AC_MSG_RESULT($psblas_cv_blacs_ok)
          LIBS="$save_LIBS";	
    fi 
fi
######################################
# Maybe we're looking at PESSL BLACS?#
######################################
if  test x"$BLACS_LIBS" != "x" ; then
    save_LIBS="$LIBS";
    LIBS="$BLACS_LIBS $LIBS"
    AC_MSG_CHECKING([for PESSL BLACS])
    AC_LANG([Fortran])
    AC_TRY_LINK_FUNC(esvemonp, [psblas_cv_pessl_blacs=yes], [psblas_cv_pessl_blacs=no])
    AC_MSG_RESULT($psblas_cv_pessl_blacs)
    LIBS="$save_LIBS";
fi    
if test "x$psblas_cv_pessl_blacs" == "xyes";  then
   FDEFINES="$psblas_cv_define_prepend-DHAVE_ESSL_BLACS $FDEFINES"
fi 
    

##############################################################################
#	Netlib BLACS library with default names
##############################################################################

if test x"$BLACS_LIBS" == "x" ; then
   save_LIBS="$LIBS";
   AC_LANG([Fortran])
   PAC_CHECK_LIBS([ blacs_MPI-LINUX-0 blacs_MPI-SP5-0 blacs_MPI-SP4-0 blacs_MPI-SP3-0 blacs_MPI-SP2-0 blacsCinit_MPI-ALPHA-0 blacsCinit_MPI-IRIX64-0 blacsCinit_MPI-RS6K-0 blacsCinit_MPI-SPP-0 blacsCinit_MPI-SUN4-0 blacsCinit_MPI-SUN4SOL2-0 blacsCinit_MPI-T3D-0 blacsCinit_MPI-T3E-0 
	], 
	[dgesd2d],
	[psblas_cv_blacs_ok=yes; LIBS="$LIBS $pac_check_libs_LIBS " 
	psblas_have_netlib_blacs=yes;  ]
	[BLACS_LIBS="$pac_check_libs_LIBS" ]
	AC_MSG_NOTICE([BLACS libraries detected.]),[]
    )
    
    if test x"$BLACS_LIBS" != "x" ; then	
      AC_LANG([Fortran])	   
      PAC_CHECK_LIBS([ blacsF77init_MPI-LINUX-0 blacsF77init_MPI-SP5-0 blacsF77init_MPI-SP4-0 blacsF77init_MPI-SP3-0 blacsF77init_MPI-SP2-0 blacsF77init_MPI-ALPHA-0 blacsF77init_MPI-IRIX64-0 blacsF77init_MPI-RS6K-0 blacsF77init_MPI-SPP-0 blacsF77init_MPI-SUN4-0 blacsF77init_MPI-SUN4SOL2-0 blacsF77init_MPI-T3D-0 blacsF77init_MPI-T3E-0 
 	], 
	[blacs_pinfo],
	[psblas_cv_blacs_ok=yes; LIBS="$pac_check_libs_LIBS $LIBS" ]
	[BLACS_LIBS="$pac_check_libs_LIBS $BLACS_LIBS" ]
	AC_MSG_NOTICE([Netlib BLACS Fortran initialization libraries detected.]),[]
       )
    fi

    if test x"$BLACS_LIBS" != "x" ; then	
    
      AC_LANG([C])
      PAC_CHECK_LIBS([ blacsCinit_MPI-LINUX-0 blacsCinit_MPI-SP5-0 blacsCinit_MPI-SP4-0 blacsCinit_MPI-SP3-0 blacsCinit_MPI-SP2-0 blacsCinit_MPI-ALPHA-0 blacsCinit_MPI-IRIX64-0 blacsCinit_MPI-RS6K-0 blacsCinit_MPI-SPP-0 blacsCinit_MPI-SUN4-0 blacsCinit_MPI-SUN4SOL2-0 blacsCinit_MPI-T3D-0 blacsCinit_MPI-T3E-0 
	], 
	[Cblacs_pinfo],
	[psblas_cv_blacs_ok=yes; LIBS="$pac_check_libs_LIBS $LIBS" ]
	[BLACS_LIBS="$BLACS_LIBS $pac_check_libs_LIBS" ]
	AC_MSG_NOTICE([Netlib BLACS C initialization libraries detected.]),[]
       )
    fi
    LIBS="$save_LIBS";	
fi

if test x"$BLACS_LIBS" == "x" ; then
	AC_MSG_ERROR([
	No BLACS library detected! $PACKAGE_NAME will be unusable.
	Please make sure a BLACS implementation is accessible (ex.: --with-blacs="-lblacsname -L/blacs/dir" )
	])
else 
      save_LIBS="$LIBS";
      LIBS="$BLACS_LIBS $LIBS"
      AC_MSG_CHECKING([for ksendid in $BLACS_LIBS])
      AC_LANG([Fortran])
      AC_TRY_LINK_FUNC(ksendid, [psblas_cv_have_sendid=yes],[psblas_cv_have_sendid=no])
      AC_MSG_RESULT($psblas_cv_have_sendid)
      LIBS="$save_LIBS"
      AC_LANG([C])
      if test "x$psblas_cv_have_sendid" == "xyes";  then
        FDEFINES="$psblas_cv_define_prepend-DHAVE_KSENDID $FDEFINES"
      fi 
fi
])dnl


dnl @synopsis PAC_MAKE_IS_GNUMAKE
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
define(PAC_MAKE_IS_GNUMAKE,[
AC_MSG_CHECKING(for gnumake)
MAKE=${MAKE:-make}

if $MAKE --version 2>&1 | grep -e"GNU Make" >/dev/null; then 
    AC_MSG_RESULT(yes)
    psblas_make_gnumake='yes'
else
    AC_MSG_RESULT(no)
    psblas_make_gnumake='no'
fi
])dnl


dnl @synopsis PAC_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl modified from ACX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the BLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/). On
dnl success, it sets the BLAS_LIBS output variable to hold the
dnl requisite library linkages.
dnl
dnl To link with BLAS, you should link with:
dnl
dnl 	$BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order. FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS), and
dnl is sometimes necessary in order to link with F77 libraries. Users
dnl will also need to use AC_F77_DUMMY_MAIN (see the autoconf manual),
dnl for the same reason.
dnl
dnl Many libraries are searched for, from ATLAS to CXML to ESSL. The
dnl user may also use --with-blas=<lib> in order to use some specific
dnl BLAS library <lib>. In order to link successfully, however, be
dnl aware that you will probably need to use the same Fortran compiler
dnl (which can be set via the F77 env. var.) as was used to compile the
dnl BLAS library.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a BLAS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands to
dnl run it if it is not found. If ACTION-IF-FOUND is not specified, the
dnl default action will define HAVE_BLAS.
dnl
dnl This macro requires autoconf 2.50 or later.
dnl
dnl @category InstalledPackages
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl @version 2001-12-13
dnl @license GPLWithACException
dnl
dnl modified by salvatore.filippone@uniroma2.it
dnl
dnl shifted check for ESSL as it was generating erroneous results on
dnl AIX SP5. 
dnl Modified with new name to handle Fortran compilers (such as NAG) 
dnl for which the linking MUST be done with the compiler (i.e.: 
dnl trying to link the Fortran version of the BLAS with the C compiler 
dnl would fail even when linking in the compiler's library)
dnl
dnl Modified by salvatore.filippone@cranfield.ac.uk to include 
dnl new tests for MKL from the 2008 version (see license below)
dnl 
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

      
AC_DEFUN([PAC_BLAS], [
AC_PREREQ(2.50)
dnl AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
pac_blas_ok=no 

AC_ARG_WITH(blas,
	[AS_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) pac_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac
AC_ARG_WITH(blasdir,
	[AS_HELP_STRING([--with-blasdir=<dir>], [search for BLAS library in <dir>])])
case $with_blasdir in
  "") ;;
      *) if test -d $with_blasdir; then
	    BLAS_LIBDIR="-L$with_blasdir";
	fi ;;
esac
# Get fortran linker names of BLAS functions to check for.
#AC_FC_FUNC(sgemm)
#AC_FC_FUNC(dgemm)

pac_blas_save_LIBS="$LIBS"
#LIBS="$LIBS $FLIBS"
AC_LANG([Fortran])

# First, check BLAS_LIBS environment variable
if test $pac_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $BLAS_LIBDIR $LIBS"
	AC_MSG_CHECKING([for sgemm in $BLAS_LIBS])
	AC_TRY_LINK_FUNC(sgemm, [pac_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($pac_blas_ok)
	LIBS="$save_LIBS"
fi
fi

LIBS="$BLAS_LIBDIR $save_LIBS "
# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $pac_blas_ok = no; then
	AC_LANG([C])
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_LANG([Fortran])
		 AC_CHECK_LIB(f77blas, sgemm,
		[AC_LANG([C])
		 AC_CHECK_LIB(cblas, cblas_dgemm,
			[pac_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])

fi
if test $pac_blas_ok = no; then
	AC_LANG([C])
	AC_CHECK_LIB(satlas, ATL_xerbla,
		[AC_LANG([Fortran])
		 AC_CHECK_LIB(satlas, sgemm,
		[AC_LANG([C])
		 AC_CHECK_LIB(satlas, cblas_dgemm,
			[pac_blas_ok=yes
			 BLAS_LIBS="-lsatlas"],
			[], [-lsatlas])],
			[], [-lsatlas])])

fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $pac_blas_ok = no; then
        AC_LANG([Fortran])
	AC_CHECK_LIB(blas, sgemm,
		[AC_CHECK_LIB(dgemm, dgemm,
		[AC_CHECK_LIB(sgemm, sgemm,
			[pac_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi


# BLAS in OpenBLAS? 
if test $pac_blas_ok = no; then
  AC_LANG([Fortran])
  AC_CHECK_LIB(openblas, sgemm, [pac_blas_ok=yes;BLAS_LIBS="-lopenblas"])
fi
				# BLAS in Intel MKL library?
sgemm="sgemm";
if test $pac_blas_ok = no; then
	# MKL for gfortran
	if test x"$ac_cv_fc_compiler_gnu" = xyes; then
		# 64 bit
		if test $host_cpu = x86_64; then
			AC_CHECK_LIB(mkl_gf_lp64, $sgemm,
			[pac_blas_ok=yes;BLAS_LIBS="-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread"],,
			[-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread])
		# 32 bit
		elif test $host_cpu = i686; then
			AC_CHECK_LIB(mkl_gf, $sgemm,
				[pac_blas_ok=yes;BLAS_LIBS="-lmkl_gf -lmkl_sequential -lmkl_core -lpthread"],,
				[-lmkl_gf -lmkl_sequential -lmkl_core -lpthread])
		fi
	# MKL for other compilers (Intel, PGI, ...?)
	else
		# 64-bit
		if test $host_cpu = x86_64; then
			AC_CHECK_LIB(mkl_intel_lp64, $sgemm,
				[pac_blas_ok=yes;BLAS_LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread"],,
				[-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread])
		# 32-bit
		elif test $host_cpu = i686; then
			AC_CHECK_LIB(mkl_intel, $sgemm,
				[pac_blas_ok=yes;BLAS_LIBS="-lmkl_intel -lmkl_sequential -lmkl_core -lpthread"],,
				[-lmkl_intel -lmkl_sequential -lmkl_core -lpthread])
		fi
	fi
fi
# Old versions of MKL
if test $pac_blas_ok = no; then
	AC_CHECK_LIB(mkl, $sgemm, [pac_blas_ok=yes;BLAS_LIBS="-lmkl -lguide -lpthread"],,[-lguide -lpthread])
fi

# BLAS in Apple vecLib library?
if test $pac_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="-framework vecLib $LIBS"
	AC_MSG_CHECKING([for $sgemm in -framework vecLib])
	AC_TRY_LINK_FUNC($sgemm, [pac_blas_ok=yes;BLAS_LIBS="-framework vecLib"])
	AC_MSG_RESULT($pac_blas_ok)
	LIBS="$save_LIBS"
fi
# BLAS in Alpha CXML library? 
if test $pac_blas_ok = no; then
	AC_CHECK_LIB(cxml, sgemm, [pac_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $pac_blas_ok = no; then
	AC_CHECK_LIB(dxml, sgemm, [pac_blas_ok=yes;BLAS_LIBS="-ldxml"])

fi

# BLAS in Sun Performance library?
if test $pac_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, sgemm,
        			[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 pac_blas_ok=yes],[],[-lsunmath])])

	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $pac_blas_ok = no; then
	AC_CHECK_LIB(scs, sgemm, [pac_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $pac_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [pac_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $pac_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, sgemm,
			[pac_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
	fi
# BLAS in generic BLAS library? 
if test $pac_blas_ok = no; then
  AC_LANG([Fortran])
  AC_CHECK_LIB(blas, sgemm, , [pac_blas_ok=yes;BLAS_LIBS="-lblas"])
fi
	
# BLAS linked to by default?  (happens on some supercomputers)
if test $pac_blas_ok = no; then
	AC_TRY_LINK_FUNC(sgemm, [pac_blas_ok=yes], [BLAS_LIBS=""])
dnl	AC_CHECK_FUNC(sgemm, [pac_blas_ok=yes])
fi

# Generic BLAS library?
if test $pac_blas_ok = no; then
  AC_LANG([Fortran])
  AC_CHECK_LIB(blas, sgemm, [pac_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

dnl AC_SUBST(BLAS_LIBS)

LIBS="$pac_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$pac_blas_ok" = xyes; then
        if  test "x$BLAS_LIBDIR" != "x" ; then
            BLAS_LIBS="$BLAS_LIBS $BLAS_LIBDIR";
        fi
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        pac_blas_ok=no
        $2
fi
])dnl PAC_BLAS


dnl @synopsis PAC_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl @synopsis ACX_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the LAPACK
dnl linear-algebra interface (see http://www.netlib.org/lapack/). On
dnl success, it sets the LAPACK_LIBS output variable to hold the
dnl requisite library linkages.
dnl
dnl To link with LAPACK, you should link with:
dnl
dnl     $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order. BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically. FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS), and
dnl is sometimes necessary in order to link with F77 libraries. Users
dnl will also need to use AC_F77_DUMMY_MAIN (see the autoconf manual),
dnl for the same reason.
dnl
dnl The user may also use --with-lapack=<lib> in order to use some
dnl specific LAPACK library <lib>. In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as was
dnl used to compile the LAPACK and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a LAPACK
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands to
dnl run it if it is not found. If ACTION-IF-FOUND is not specified, the
dnl default action will define HAVE_LAPACK.
dnl
dnl @category InstalledPackages
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl @version 2002-03-12
dnl @license GPLWithACException
dnl modified by salvatore.filippone@uniroma2.it
dnl shifted check for ESSL as it was generating erroneous results on
dnl AIX SP5. 
dnl Modified with new name to handle Fortran compilers (such as NAG) 
dnl for which the linking MUST be done with the compiler (i.e.: 
dnl trying to link the Fortran version of the BLAS with the C compiler 
dnl would fail even when linking in the compiler's library)

AC_DEFUN([PAC_LAPACK], [
AC_REQUIRE([PAC_BLAS])
pac_lapack_ok=no

AC_ARG_WITH(lapack,
        [AS_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
case $with_lapack in
        yes | "") ;;
        no) pac_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
        *) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
#AC_FC_FUNC(cheev)

# We cannot use LAPACK if BLAS is not found
if test "x$pac_blas_ok" != xyes; then
        pac_lapack_ok=noblas
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for cheev in $LAPACK_LIBS])
	AC_LANG([Fortran])
	dnl Warning : square brackets are EVIL!
	cat > conftest.$ac_ext <<EOF
        program test_cheev 
          call cheev
        end 
EOF
	if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
	  pac_lapack_ok=yes
	  AC_MSG_RESULT([yes])	
	else
	  AC_MSG_RESULT([no])	
	  echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
	  cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
	fi 
	rm -f conftest*
        LIBS="$save_LIBS"
        if test pac_lapack_ok = no; then
                LAPACK_LIBS=""
        fi
        AC_LANG([C])
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $pac_lapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
        AC_MSG_CHECKING([for cheev in default libs])
	AC_LANG([Fortran])
	dnl Warning : square brackets are EVIL!
	cat > conftest.$ac_ext <<EOF
        program test_cheev 
          call cheev
        end 
EOF
	if AC_TRY_EVAL(ac_link) && test -s conftest${ac_exeext}; then
	  pac_lapack_ok=yes
	  AC_MSG_RESULT([yes])	
	else
	  AC_MSG_RESULT([no])	
	  echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
	  cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
	fi 
	rm -f conftest*
        LIBS="$save_LIBS"
        AC_LANG([C])
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
        if test $pac_lapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
		AC_LANG([Fortran])
		AC_CHECK_LIB($lapack, cheev,
                    [pac_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$FLIBS])
		AC_LANG([C])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$pac_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        pac_lapack_ok=no
        $2
fi
])dnl PAC_LAPACK

dnl @synopsis PAC_FORTRAN_TEST_FLUSH( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the FLUSH Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_FLUSH,
dnl Warning : square brackets are EVIL!
[AC_MSG_CHECKING([support for Fortran FLUSH statement])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program conftest
   integer :: iunit=10
   open(10)
   write(10,*) 'Test '
   flush(10)
   close(10)
end program conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_TEST_ISO_FORTRAN_ENV( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will determine if the fortran compiler MPIFC supports ISO_FORTRAN_ENV
dnl
dnl If yes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl 
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_FORTRAN_TEST_ISO_FORTRAN_ENV,
[AC_MSG_CHECKING([support for ISO_FORTRAN_ENV])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
           program test
             use iso_fortran_env
           end program test],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_TEST_ISO_C_BIND( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the ISO C Binding  Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Michele Martone <michele.martone@uniroma2.it>
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_ISO_C_BIND,
[AC_MSG_CHECKING([support for Fortran ISO_C_BINDING module])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program conftest
  use iso_c_binding
end program conftest],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_FORTRAN_TEST_MOLD( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the MOLD=  Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_MOLD,
[AC_MSG_CHECKING([support for Fortran MOLD= allocation])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program xtt
  type foo
    integer :: i
  end type foo
  type, extends(foo) :: new_foo
    integer :: j
  end type new_foo
  class(foo), allocatable  :: fooab
  type(new_foo) :: nfv 
  integer :: info

  allocate(fooab, mold=nfv, stat=info)

end program xtt],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])


dnl @synopsis PAC_FORTRAN_TEST_SOURCE( [ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl
dnl Will try to compile and link a program checking the SOURCE=  Fortran support.
dnl
dnl Will use MPIFC, otherwise '$FC'.
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
AC_DEFUN(PAC_FORTRAN_TEST_SOURCE,
[AC_MSG_CHECKING([support for Fortran SOURCE= allocation])
 AC_LANG_PUSH([Fortran])
 ac_exeext=''
 ac_ext='f90'
 dnl ac_link='${MPIFC-$FC} -o conftest${ac_exeext} $FFLAGS $LDFLAGS conftest.$ac_ext $LIBS 1>&5'
 ac_fc=${MPIFC-$FC};
 AC_COMPILE_IFELSE([
program xtt
  type foo
    integer :: i
  end type foo
  type, extends(foo) :: new_foo
    integer :: j
  end type new_foo
  class(foo), allocatable  :: fooab
  type(new_foo) :: nfv 
  integer :: info

  allocate(fooab, source=nfv, stat=info)

end program xtt],
		  [  AC_MSG_RESULT([yes])
		     ifelse([$1], , :, [ $1])],
		  [  AC_MSG_RESULT([no])	
		     echo "configure: failed program was:" >&AS_MESSAGE_LOG_FD
		     cat conftest.$ac_ext >&AS_MESSAGE_LOG_FD  
		     ifelse([$2], , , [ $2])])
AC_LANG_POP([Fortran])
])

dnl @synopsis PAC_CHECK_AMD
dnl
dnl Will try to find the AMD library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_AMD,
[AC_ARG_WITH(amd, AS_HELP_STRING([--with-amd=LIBNAME], [Specify the library name for AMD library. 
Default: "-lamd"]),
        [psblas_cv_amd=$withval],
        [psblas_cv_amd='-lamd'])
AC_ARG_WITH(amddir, AS_HELP_STRING([--with-amddir=DIR], [Specify the directory for AMD library and includes.]),
        [psblas_cv_amddir=$withval],
        [psblas_cv_amddir=''])
AC_ARG_WITH(amdincdir, AS_HELP_STRING([--with-amdincdir=DIR], [Specify the directory for AMD includes.]),
        [psblas_cv_amdincdir=$withval],
        [psblas_cv_amdincdir=''])
AC_ARG_WITH(amdlibdir, AS_HELP_STRING([--with-amdlibdir=DIR], [Specify the directory for AMD library.]),
        [psblas_cv_amdlibdir=$withval],
        [psblas_cv_amdlibdir=''])

AC_LANG([C])
SAVE_LIBS="$LIBS"
SAVE_CPPFLAGS="$CPPFLAGS"
if test "x$psblas_cv_amddir" != "x"; then 
   AMD_LIBDIR="-L$psblas_cv_amddir"
   LIBS="-L$psblas_cv_amddir $LIBS"
   AMD_INCLUDES="-I$psblas_cv_amddir"
   CPPFLAGS="$AMD_INCLUDES $CPPFLAGS"
fi
if test "x$psblas_cv_amdincdir" != "x"; then 
   AMD_INCLUDES="-I$psblas_cv_amdincdir"
   CPPFLAGS="$AMD_INCLUDES $CPPFLAGS"
fi
if test "x$psblas_cv_amdlibdir" != "x"; then 
   LIBS="-L$psblas_cv_amdlibdir $LIBS"
   AMD_LIBDIR="-L$psblas_cv_amdlibdir"
fi

AC_MSG_NOTICE([amd dir $psblas_cv_amddir])
AC_CHECK_HEADER([amd.h],
 [pac_amd_header_ok=yes],
 [pac_amd_header_ok=no; AMD_INCLUDES=""])
if test "x$pac_amd_header_ok" == "xno" ; then
dnl Maybe Include or include subdirs? 
  unset ac_cv_header_amd_h
  AMD_INCLUDES="-I$psblas_cv_amddir/include -I$psblas_cv_amddir/Include "
  CPPFLAGS="$AMD_INCLUDES $SAVE_CPPFLAGS"

 AC_MSG_CHECKING([for amd_h in $AMD_INCLUDES])
 AC_CHECK_HEADER([amd.h],
    [pac_amd_header_ok=yes],
    [pac_amd_header_ok=no; AMD_INCLUDES=""])
fi
if test "x$pac_amd_header_ok" == "xno" ; then
dnl Maybe new structure with AMD UFconfig AMD? 
   unset ac_cv_header_amd_h
   AMD_INCLUDES="-I$psblas_cv_amddir/UFconfig -I$psblas_cv_amddir/AMD/Include -I$psblas_cv_amddir/AMD/Include"
   CPPFLAGS="$AMD_INCLUDES $SAVE_CPPFLAGS"
   AC_CHECK_HEADER([amd.h],
     [pac_amd_header_ok=yes],
     [pac_amd_header_ok=no; AMD_INCLUDES=""])
fi


if test "x$pac_amd_header_ok" == "xyes" ; then 
      psblas_cv_amd_includes="$AMD_INCLUDES"
      if test "x$AMD_LIBDIR" == "x" ; then 
	 AMD_LIBS="$psblas_cv_amd"
      else
	AMD_LIBS="$psblas_cv_amd $AMD_LIBDIR"
      fi
      LIBS="$AMD_LIBS -lm $LIBS";
      AC_MSG_CHECKING([for amd_order in $AMD_LIBS])
      AC_TRY_LINK_FUNC(amd_order, 
       [psblas_cv_have_amd=yes;pac_amd_lib_ok=yes; ],
       [psblas_cv_have_amd=no;pac_amd_lib_ok=no; AMD_LIBS=""])
      AC_MSG_RESULT($pac_amd_lib_ok)
     if test "x$pac_amd_lib_ok" == "xno" ; then 
        dnl Maybe Lib or lib? 
        AMD_LIBDIR="-L$psblas_cv_amddir/Lib -L$psblas_cv_amddir/lib"
        AMD_LIBS="$psblas_cv_amd $AMD_LIBDIR"
        LIBS="$AMD_LIBS -lm $SAVE_LIBS"
        
      AC_MSG_CHECKING([for amd_order in $AMD_LIBS])
      AC_TRY_LINK_FUNC(amd_order, 
       [psblas_cv_have_amd=yes;pac_amd_lib_ok=yes; ],
       [psblas_cv_have_amd=no;pac_amd_lib_ok=no; AMD_LIBS=""])
      AC_MSG_RESULT($pac_amd_lib_ok)
     fi
     if test "x$pac_amd_lib_ok" == "xno" ; then 
        dnl Maybe AMD/Lib? 
        AMD_LIBDIR="-L$psblas_cv_amddir/AMD/Lib -L$psblas_cv_amddir/AMD/Lib"
        AMD_LIBS="$psblas_cv_amd $AMD_LIBDIR"
        LIBS="$AMD_LIBS -lm $SAVE_LIBS"
      AC_MSG_CHECKING([for amd_order in $AMD_LIBS])
      AC_TRY_LINK_FUNC(amd_order, 
       [psblas_cv_have_amd=yes;pac_amd_lib_ok=yes; ],
       [psblas_cv_have_amd=no;pac_amd_lib_ok=no; AMD_LIBS=""])
      AC_MSG_RESULT($pac_amd_lib_ok)
     fi
fi
LIBS="$SAVE_LIBS";
CPPFLAGS="$SAVE_CPPFLAGS";
])dnl 

dnl @synopsis PAC_CHECK_METIS
dnl
dnl Will try to find the METIS library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_METIS,
[AC_ARG_WITH(metis, AS_HELP_STRING([--with-metis=LIBNAME], [Specify the library name for METIS library. 
Default: "-lmetis"]),
        [psblas_cv_metis=$withval],
        [psblas_cv_metis='-lmetis'])
AC_ARG_WITH(metisincfile, AS_HELP_STRING([--with-metisincfile=DIR], [Specify the name  for METIS include file.]),
        [psblas_cv_metisincfile=$withval],
        [psblas_cv_metisincfile='metis.h'])
AC_ARG_WITH(metisdir, AS_HELP_STRING([--with-metisdir=DIR], [Specify the directory for METIS library and includes.]),
        [psblas_cv_metisdir=$withval],
        [psblas_cv_metisdir=''])
AC_ARG_WITH(metisincdir, AS_HELP_STRING([--with-metisincdir=DIR], [Specify the directory for METIS includes.]),
        [psblas_cv_metisincdir=$withval],
        [psblas_cv_metisincdir=''])
AC_ARG_WITH(metislibdir, AS_HELP_STRING([--with-metislibdir=DIR], [Specify the directory for METIS library.]),
        [psblas_cv_metislibdir=$withval],
        [psblas_cv_metislibdir=''])

AC_LANG([C])
SAVE_LIBS="$LIBS"
SAVE_CPPFLAGS="$CPPFLAGS"
if test "x$psblas_cv_metisdir" != "x"; then 
   METIS_LIBDIR="-L$psblas_cv_metisdir"
   LIBS="-L$psblas_cv_metisdir $LIBS"
   METIS_INCLUDES="-I$psblas_cv_metisdir"
   CPPFLAGS="$METIS_INCLUDES $CPPFLAGS"
fi
if test "x$psblas_cv_metisincdir" != "x"; then 
   METIS_INCLUDES="-I$psblas_cv_metisincdir"
   CPPFLAGS="$METIS_INCLUDES $CPPFLAGS"
fi
if test "x$psblas_cv_metislibdir" != "x"; then 
   LIBS="-L$psblas_cv_metislibdir $LIBS"
   METIS_LIBDIR="-L$psblas_cv_metislibdir"
fi

AC_MSG_NOTICE([metis dir $psblas_cv_metisdir])
AC_CHECK_HEADERS([limits.h "$psblas_cv_metisincfile"],
 [pac_metis_header_ok=yes],
 [pac_metis_header_ok=no; METIS_INCLUDES=""])
if test "x$pac_metis_header_ok" == "xno" ; then
dnl Maybe Include or include subdirs? 
  unset ac_cv_header_metis_h
  METIS_INCLUDES="-I$psblas_cv_metisdir/include -I$psblas_cv_metisdir/Include "
  CPPFLAGS="$METIS_INCLUDES $SAVE_CPPFLAGS"

 AC_MSG_CHECKING([for metis_h in $METIS_INCLUDES])
 AC_CHECK_HEADERS([limits.h  "$psblas_cv_metisincfile"],
    [pac_metis_header_ok=yes],
    [pac_metis_header_ok=no; METIS_INCLUDES=""])
fi
if test "x$pac_metis_header_ok" == "xno" ; then
dnl Maybe new structure with METIS UFconfig METIS? 
   unset ac_cv_header_metis_h
   METIS_INCLUDES="-I$psblas_cv_metisdir/UFconfig -I$psblas_cv_metisdir/METIS/Include -I$psblas_cv_metisdir/METIS/Include"
   CPPFLAGS="$METIS_INCLUDES $SAVE_CPPFLAGS"
   AC_CHECK_HEADERS([limits.h  "$psblas_cv_metisincfile"],
     [pac_metis_header_ok=yes],
     [pac_metis_header_ok=no; METIS_INCLUDES=""])
fi

if test "x$pac_metis_header_ok" == "xyes" ; then
   AC_LANG_PUSH([C])
   AC_MSG_CHECKING([for METIS integer size])
   AC_LINK_IFELSE([AC_LANG_SOURCE(
	#include <stdio.h>
	#include "$psblas_cv_metisincfile"
        void main(){
		    printf("%d\n",IDXTYPEWIDTH);
		    }
	       )],
	       [pac_cv_metis_idx=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [pac_cv_metis_idx="unknown"])
      AC_MSG_RESULT($pac_cv_metis_idx)

   AC_LANG_POP()
fi

if test "x$pac_metis_header_ok" == "xyes" ; then
   AC_LANG_PUSH([C])
   AC_MSG_CHECKING([for METIS real size])
   AC_LINK_IFELSE([AC_LANG_SOURCE(
	#include <stdio.h>
	#include "$psblas_cv_metisincfile"
        void main(){
		    printf("%d\n",REALTYPEWIDTH);
		    }
	       )],
	       [pac_cv_metis_real=`./conftest${ac_exeext} | sed 's/^ *//'`],
	       [pac_cv_metis_real="unknown"])
      AC_MSG_RESULT($pac_cv_metis_real)

   AC_LANG_POP()
fi

if test "x$pac_metis_header_ok" = "xyes" ; then 
      psblas_cv_metis_includes="$METIS_INCLUDES"
      if  test "x$METIS_LIBDIR" == "x" ; then
	  METIS_LIBS="$psblas_cv_metis"
      else
	METIS_LIBS="$psblas_cv_metis $METIS_LIBDIR"
      fi  
      LIBS="$METIS_LIBS -lm $LIBS";
      AC_MSG_CHECKING([for METIS_PartGraphKway in $METIS_LIBS])
      AC_TRY_LINK_FUNC(METIS_PartGraphKway, 
       [psblas_cv_have_metis=yes;pac_metis_lib_ok=yes; ],
       [psblas_cv_have_metis=no;pac_metis_lib_ok=no; METIS_LIBS=""])
      AC_MSG_RESULT($pac_metis_lib_ok)
     if test "x$pac_metis_lib_ok" = "xno" ; then 
        dnl Maybe Lib or lib? 
        METIS_LIBDIR="-L$psblas_cv_metisdir/Lib -L$psblas_cv_metisdir/lib"
        METIS_LIBS="$psblas_cv_metis $METIS_LIBDIR"
        LIBS="$METIS_LIBS -lm $SAVE_LIBS"
        
      AC_MSG_CHECKING([for METIS_PartGraphKway in $METIS_LIBS])
      AC_TRY_LINK_FUNC(METIS_PartGraphKway, 
       [psblas_cv_have_metis=yes;pac_metis_lib_ok=yes; ],
       [psblas_cv_have_metis=no;pac_metis_lib_ok=no; METIS_LIBS=""])
      AC_MSG_RESULT($pac_metis_lib_ok)
     fi

     if test "x$pac_metis_lib_ok" = "xno" ; then 
        dnl Maybe METIS/Lib? 
        METIS_LIBDIR="-L$psblas_cv_metisdir/METIS/Lib -L$psblas_cv_metisdir/METIS/Lib"
        METIS_LIBS="$psblas_cv_metis $METIS_LIBDIR"
        LIBS="$METIS_LIBS -lm $SAVE_LIBS"
      AC_MSG_CHECKING([for METIS_PartGraphKway in $METIS_LIBS])
      AC_TRY_LINK_FUNC(METIS_PartGraphKway, 
       [psblas_cv_have_metis=yes;pac_metis_lib_ok="yes"; ],
       [psblas_cv_have_metis=no;pac_metis_lib_ok="no"; METIS_LIBS=""])
      AC_MSG_RESULT($pac_metis_lib_ok)
      fi
 fi
dnl AC_MSG_NOTICE([ metis lib ok $pac_metis_lib_ok])

 if test "x$pac_metis_lib_ok" = "xyes" ; then 
      AC_MSG_CHECKING([for METIS_SetDefaultOptions in $LIBS])
      AC_TRY_LINK_FUNC(METIS_SetDefaultOptions, 
        [psblas_cv_have_metis=yes;pac_metis_lib_ok=yes; ],
        [psblas_cv_have_metis=no;pac_metis_lib_ok="no. Unusable METIS version, sorry."; METIS_LIBS=""
      ])
      AC_MSG_RESULT($pac_metis_lib_ok)

fi

LIBS="$SAVE_LIBS";
CPPFLAGS="$SAVE_CPPFLAGS";
])dnl 


dnl @synopsis PAC_CHECK_SPGPU
dnl
dnl Will try to find the spgpu library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_SPGPU,
	 [SAVE_LIBS="$LIBS"
	  SAVE_CPPFLAGS="$CPPFLAGS"
	  if test "x$pac_cv_have_cuda" == "x"; then  
             PAC_CHECK_CUDA()
          fi
dnl AC_MSG_NOTICE([From CUDA: $pac_cv_have_cuda ])
	  if test "x$pac_cv_have_cuda" == "xyes"; then  
	  AC_ARG_WITH(spgpu, AC_HELP_STRING([--with-spgpu=DIR], [Specify the directory for SPGPU library and includes.]),
		      [pac_cv_spgpudir=$withval],
		      [pac_cv_spgpudir=''])
	  
	  AC_LANG([C])
	  if test "x$pac_cv_spgpudir" != "x"; then 
	  LIBS="-L$pac_cv_spgpudir/lib $LIBS"
	  GPU_INCLUDES="-I$pac_cv_spgpudir/include"
	  CPPFLAGS="$GPU_INCLUDES $CUDA_INCLUDES $CPPFLAGS"
	  GPU_LIBDIR="-L$pac_cv_spgpudir/lib"
	  fi
	  AC_MSG_CHECKING([spgpu dir $pac_cv_spgpudir])
	  AC_CHECK_HEADER([core.h],
			  [pac_gpu_header_ok=yes],
			  [pac_gpu_header_ok=no; GPU_INCLUDES=""])
	  
	  if test "x$pac_gpu_header_ok" == "xyes" ; then 
	  GPU_LIBS="-lspgpu $GPU_LIBDIR"
	  LIBS="$GPU_LIBS $CUDA_LIBS -lm $LIBS";
	  AC_MSG_CHECKING([for spgpuCreate in $GPU_LIBS])
	  AC_TRY_LINK_FUNC(spgpuCreate, 
			   [pac_cv_have_spgpu=yes;pac_gpu_lib_ok=yes; ],
			   [pac_cv_have_spgpu=no;pac_gpu_lib_ok=no; GPU_LIBS=""])
	  AC_MSG_RESULT($pac_gpu_lib_ok)
	  if test "x$pac_cv_have_spgpu" == "xyes" ; then 
	  AC_MSG_NOTICE([Have found SPGPU])
	  SPGPULIBNAME="libpsbgpu.a";
	  SPGPU_DIR="$pac_cv_spgpudir";
	  SPGPU_DEFINES="-DHAVE_SPGPU";
	  SPGPU_INCDIR="$SPGPU_DIR/include";
	  SPGPU_INCLUDES="-I$SPGPU_INCDIR";
	  SPGPU_LIBS="-lspgpu -L$SPGPU_DIR/lib";
	  LGPU=-lpsb_gpu
	  CUDA_DIR="$pac_cv_cuda_dir";
	  CUDA_DEFINES="-DHAVE_CUDA";
	  CUDA_INCLUDES="-I$pac_cv_cuda_dir/include"
	  CUDA_LIBDIR="-L$pac_cv_cuda_dir/lib64 -L$pac_cv_cuda_dir/lib"
	  FDEFINES="$psblas_cv_define_prepend-DHAVE_GPU $psblas_cv_define_prepend-DHAVE_SPGPU $psblas_cv_define_prepend-DHAVE_CUDA $FDEFINES";
	  CDEFINES="-DHAVE_SPGPU -DHAVE_CUDA $CDEFINES" ;
	  fi
  fi
fi
LIBS="$SAVE_LIBS"
CPPFLAGS="$SAVE_CPPFLAGS"
])dnl 




dnl @synopsis PAC_CHECK_CUDA
dnl
dnl Will try to find the cuda library and headers.
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_CUDA,
[AC_ARG_WITH(cuda, AC_HELP_STRING([--with-cuda=DIR], [Specify the directory for CUDA library and includes.]),
        [pac_cv_cuda_dir=$withval],
        [pac_cv_cuda_dir=''])

AC_LANG([C])
SAVE_LIBS="$LIBS"
SAVE_CPPFLAGS="$CPPFLAGS"
if test "x$pac_cv_cuda_dir" != "x"; then 
   CUDA_DIR="$pac_cv_cuda_dir"
   LIBS="-L$pac_cv_cuda_dir/lib $LIBS"
   CUDA_INCLUDES="-I$pac_cv_cuda_dir/include"
   CUDA_DEFINES="-DHAVE_CUDA"
   CPPFLAGS="$CUDA_INCLUDES $CPPFLAGS"
   CUDA_LIBDIR="-L$pac_cv_cuda_dir/lib64 -L$pac_cv_cuda_dir/lib"
   if test -f "$pac_cv_cuda_dir/bin/nvcc"; then
     CUDA_NVCC="$pac_cv_cuda_dir/bin/nvcc"
   else
     CUDA_NVCC="nvcc"
   fi
fi
AC_MSG_CHECKING([cuda dir $pac_cv_cuda_dir])
AC_CHECK_HEADER([cuda_runtime.h],
 [pac_cuda_header_ok=yes],
 [pac_cuda_header_ok=no; CUDA_INCLUDES=""])

if test "x$pac_cuda_header_ok" == "xyes" ; then 
 CUDA_LIBS="-lcusparse -lcublas -lcudart $CUDA_LIBDIR"
 LIBS="$CUDA_LIBS -lm $LIBS";
 AC_MSG_CHECKING([for cudaMemcpy in $CUDA_LIBS])
 AC_TRY_LINK_FUNC(cudaMemcpy, 
		  [pac_cv_have_cuda=yes;pac_cuda_lib_ok=yes; ],
		  [pac_cv_have_cuda=no;pac_cuda_lib_ok=no; CUDA_LIBS=""])
 AC_MSG_RESULT($pac_cuda_lib_ok)

fi
LIBS="$SAVE_LIBS"
CPPFLAGS="$SAVE_CPPFLAGS"
])dnl 

dnl @synopsis PAC_ARG_WITH_CUDACC
dnl
dnl Test for --with-cudacc="set_of_cc".
dnl 
dnl Defines the CC to compile for
dnl
dnl
dnl Example use:
dnl
dnl PAC_ARG_WITH_CUDACC
dnl 
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN([PAC_ARG_WITH_CUDACC],
[
AC_ARG_WITH(cudacc,
AC_HELP_STRING([--with-cudacc], [A comma-separated list of CCs to compile to, for example,
 --with-cudacc=30,35,37,50,60]),
[pac_cv_cudacc=$withval],
[pac_cv_cudacc=''])
])

AC_DEFUN(PAC_ARG_WITH_LIBRSB,
	 [SAVE_LIBS="$LIBS"
	  SAVE_CPPFLAGS="$CPPFLAGS"

	  AC_ARG_WITH(librsb,
	  AC_HELP_STRING([--with-librsb], [The directory for LIBRSB, for example,
 	  --with-librsb=/opt/packages/librsb]),
	  [pac_cv_librsb_dir=$withval],
	  [pac_cv_librsb_dir=''])
	  
	  if test "x$pac_cv_librsb_dir" != "x"; then 
	  LIBS="-L$pac_cv_librsb_dir $LIBS"
	  RSB_INCLUDES="-I$pac_cv_librsb_dir"
	  # CPPFLAGS="$GPU_INCLUDES $CUDA_INCLUDES $CPPFLAGS"
	  RSB_LIBDIR="-L$pac_cv_librsb_dir"
	  fi
	  #AC_MSG_CHECKING([librsb dir $pac_cv_librsb_dir])
	  AC_CHECK_HEADER([$pac_cv_librsb_dir/rsb.h],
			  [pac_rsb_header_ok=yes],
			  [pac_rsb_header_ok=no; RSB_INCLUDES=""])
	  
	  if test "x$pac_rsb_header_ok" == "xyes" ; then 
	  RSB_LIBS="-lrsb $RSB_LIBDIR"
	  # LIBS="$GPU_LIBS $CUDA_LIBS -lm $LIBS";
	  # AC_MSG_CHECKING([for spgpuCreate in $GPU_LIBS])
	  # AC_TRY_LINK_FUNC(spgpuCreate, 
	  # 		   [pac_cv_have_spgpu=yes;pac_gpu_lib_ok=yes; ],
	  # 		   [pac_cv_have_spgpu=no;pac_gpu_lib_ok=no; GPU_LIBS=""])
	  # AC_MSG_RESULT($pac_gpu_lib_ok)
	  # if test "x$pac_cv_have_spgpu" == "xyes" ; then 
	  # AC_MSG_NOTICE([Have found SPGPU])
	  RSBLIBNAME="librsb.a";
	  LIBRSB_DIR="$pac_cv_librsb_dir";
	  # SPGPU_DEFINES="-DHAVE_SPGPU";
	  LIBRSB_INCDIR="$LIBRSB_DIR";
	  LIBRSB_INCLUDES="-I$LIBRSB_INCDIR";
	  LIBRSB_LIBS="-lrsb -L$LIBRSB_DIR";
	  # CUDA_DIR="$pac_cv_cuda_dir";
	  LIBRSB_DEFINES="-DHAVE_RSB";
	  LRSB=-lpsb_rsb
	  # CUDA_INCLUDES="-I$pac_cv_cuda_dir/include"
	  # CUDA_LIBDIR="-L$pac_cv_cuda_dir/lib64 -L$pac_cv_cuda_dir/lib"
	  FDEFINES="$LIBRSB_DEFINES $psblas_cv_define_prepend $FDEFINES";
	  CDEFINES="$LIBRSB_DEFINES $CDEFINES";#CDEFINES="-DHAVE_SPGPU -DHAVE_CUDA $CDEFINES";
	  fi
#  fi
LIBS="$SAVE_LIBS"
CPPFLAGS="$SAVE_CPPFLAGS"
])
dnl

dnl @synopsis PAC_CHECK_CUDA_VERSION
dnl
dnl Will try to find the cuda version
dnl
dnl Will use $CC
dnl
dnl If the test passes, will execute ACTION-IF-FOUND. Otherwise, ACTION-IF-NOT-FOUND.
dnl Note : This file will be likely to induce the compiler to create a module file
dnl (for a module called conftest).
dnl Depending on the compiler flags, this could cause a conftest.mod file to appear
dnl in the present directory, or in another, or with another name. So be warned!
dnl
dnl @author Salvatore Filippone <salvatore.filippone@uniroma2.it>
dnl
AC_DEFUN(PAC_CHECK_CUDA_VERSION,
[AC_LANG_PUSH([C])
SAVE_LIBS="$LIBS"
SAVE_CPPFLAGS="$CPPFLAGS"
if test "x$pac_cv_have_cuda" == "x"; then  
        PAC_CHECK_CUDA()
fi
if test "x$pac_cv_have_cuda" == "xyes"; then
   CUDA_DIR="$pac_cv_cuda_dir"
   LIBS="-L$pac_cv_cuda_dir/lib $LIBS"
   CUDA_INCLUDES="-I$pac_cv_cuda_dir/include"
   CUDA_DEFINES="-DHAVE_CUDA"
   CPPFLAGS="$CUDA_INCLUDES $CPPFLAGS"
   CUDA_LIBDIR="-L$pac_cv_cuda_dir/lib64 -L$pac_cv_cuda_dir/lib"
  CUDA_LIBS="-lcusparse -lcublas -lcudart $CUDA_LIBDIR"
  LIBS="$CUDA_LIBS -lm $LIBS";
  AC_MSG_CHECKING([for CUDA version])
  AC_LINK_IFELSE([AC_LANG_SOURCE([
#include <stdio.h>
#include <cuda.h>

int main(int argc, char *argv[])
{
  printf("%d",CUDA_VERSION);
  return(0);
} ])],
	[pac_cv_cuda_version=`./conftest${ac_exeext} | sed 's/^ *//'`;],
	[pac_cv_cuda_version="unknown";])
 
 AC_MSG_RESULT($pac_cv_cuda_version)
 fi
AC_LANG_POP([C]) 
LIBS="$SAVE_LIBS"
CPPFLAGS="$SAVE_CPPFLAGS"
])dnl 


