# AC_OPENACC
# ---------
# Check which options need to be passed to the C compiler to support Openacc.
# Set the OPENACC_CFLAGS / OPENACC_CXXFLAGS / OPENACC_FFLAGS variable to these
# options.
# The options are necessary at compile time (so the #pragmas are understood)
# and at link time (so the appropriate library is linked with).
# This macro takes care to not produce redundant options if $CC $CFLAGS already
# supports Openacc.
#
# For each candidate option, we do a compile test first, then a link test;
# if the compile test succeeds but the link test fails, that means we have
# found the correct option but it doesn't work because the libraries are
# broken.  (This can happen, for instance, with SunPRO C and a bad combination
# of operating system patches.)
#
# Several of the options in our candidate list can be misinterpreted by
# compilers that don't use them to activate Openacc support; for example,
# many compilers understand "-openacc" to mean "write output to a file
# named 'penmp'" rather than "enable Openacc".  We can't completely avoid
# the possibility of clobbering files named 'penmp' or 'mp' in configure's
# working directory; therefore, this macro will bomb out if any such file
# already exists when it's invoked.
AC_DEFUN([AX_C_OPENACC],
[AC_REQUIRE([_AX_OPENACC_SAFE_WD])]dnl
[AC_ARG_ENABLE([openacc],
   [AS_HELP_STRING([--disable-openacc], [do not use Openacc])])]dnl
[
  OPENACC_[]_AC_LANG_PREFIX[]FLAGS=
  if test "$enable_openacc" != no; then
      AC_LANG_PUSH([C])
      AC_CACHE_CHECK([for $[]_AC_CC[] option to support Openacc],
      [ax_cv_prog_[]_AC_LANG_ABBREV[]_openacc],
      [ax_cv_prog_[]_AC_LANG_ABBREV[]_openacc='not found'
      dnl Try these flags:
      dnl   (on by default)      ''
      dnl   GCC >= 4.2           -fopenacc
      dnl   SunPRO C             -xopenacc
      dnl   Intel C              -openacc
      dnl   SGI C, PGI C         -mp
      dnl   Tru64 Compaq C       -omp
      dnl   IBM XL C (AIX, Linux) -qsmp=omp
      dnl   Cray CCE             -homp
      dnl   NEC SX               -Popenacc
      dnl   Lahey Fortran (Linux)  --openacc
      for ac_option in '' -fopenacc -openacc -acc; do

        ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
        _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $ac_option"
        AC_COMPILE_IFELSE([
#ifndef _OPENACC
#error "OpenACC not supported"
#endif
#include <openacc.h>
  int main (void) { acc_init (0); return 0;}
],
          [AC_LINK_IFELSE([
#ifndef _OPENACC
#error "OpenACC not supported"
#endif
#include <openacc.h>
 int main (void) { acc_init (0); return 0;}
],
            [ax_cv_prog_[]_AC_LANG_ABBREV[]_openacc=$ac_option],
            [ax_cv_prog_[]_AC_LANG_ABBREV[]_openacc='unsupported'])])
        _AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS

        if test "$ax_cv_prog_[]_AC_LANG_ABBREV[]_openacc" != 'not found'; then
          break
        fi
      done
      if test "$ax_cv_prog_[]_AC_LANG_ABBREV[]_openacc" = 'not found'; then
        ac_cv_prog_[]_AC_LANG_ABBREV[]_openacc='unsupported'
      elif test "$ax_cv_prog_[]_AC_LANG_ABBREV[]_openacc" = ''; then
        ac_cv_prog_[]_AC_LANG_ABBREV[]_openacc='none needed'
      fi
      dnl _AX_OPENACC_SAFE_WD checked that these files did not exist before we
      dnl started probing for Openacc support, so if they exist now, they were
      dnl created by the probe loop and it's safe to delete them.
      rm -f penmp mp])
    if test "$ax_cv_prog_[]_AC_LANG_ABBREV[]_openacc" != 'unsupported' && \
       test "$ax_cv_prog_[]_AC_LANG_ABBREV[]_openacc" != 'none needed'; then
      OPENACC_[]_AC_LANG_PREFIX[]FLAGS="$ax_cv_prog_[]_AC_LANG_ABBREV[]_openacc"
    fi
   AC_LANG_POP([C])
  fi
])

# _AC_OPENACC_SAFE_WD
# ------------------
# AC_REQUIREd by AC_OPENACC.  Checks both at autoconf time and at
# configure time for files that AC_OPENACC clobbers.
AC_DEFUN([_AX_OPENACC_SAFE_WD],
[m4_syscmd([test ! -e penmp && test ! -e mp])]dnl
[m4_if(sysval, [0], [], [m4_fatal(m4_normalize(
  [AX_OPENACC clobbers files named 'mp' and 'penmp'.
   To use AX_OPENACC you must not have either of these files
   at the top level of your source tree.]))])]dnl
[if test -e penmp || test -e mp; then
  AC_MSG_ERROR(m4_normalize(
    [AX@&t@_OPENACC clobbers files named 'mp' and 'penmp'.
     Aborting configure because one of these files already exists.]))
fi])

