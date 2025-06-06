!> Test program for y = alpha * x + betha * y  psb_geaxbpy routine
!! Check the README.md to see all details about the tests.
!!
!! Authors: Luca Pepé Sciarria, Staccone Simone (Tor Vergata University)
!! 
!! Type: Synchronous.
!! 
!! ROUTINE PARAMETERS
!! alpha    Description: the scalar α.
!!          Scope: global
!!          Type: required
!!          Intent: in
!!          Specified as: a number of the data type indicated in Table 1.
!!
!! x        Description: the local portion of global dense matrix x.
!!          Scope: local
!!          Type: required
!!          Intent: in
!!          Specified as: a rank one or two array or an object of type psb T vect type
!!          containing numbers of type specified in Table 1. The rank of x must be
!!          the same of y.
!!
!! beta     Description: the scalar β.
!!          Scope: global
!!          Type: required
!!          Intent: in.
!!          Specified as: a number of the data type indicated in Table 1.
!! 
!! y        Description: the local portion of the global dense matrix y.
!!          Scope: local
!!          Type: required
!!          Intent: inout
!!          Specified as: a rank one or two array or an object of type psb T vect type
!!          containing numbers of the type indicated in Table 1. The rank of y must
!!          be the same of x.
!!
!! desc_a   Description: contains data structures for communications.
!!          Scope: local
!!          Type: required
!!          Intent: in
!!          Specified as: an object of type psb desc type.
!!
!! info     Description: Error code.
!!          Scope: local
!!          Type: required
!!          Intent: out.
!!          Specified as: An integer value; 0 means no error has been detected.
!!
module psb_geaxpby_test
    contains
    subroutine psb_spmm_kernel()
    end subroutine

end module psb_geaxpby_test