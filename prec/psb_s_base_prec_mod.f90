!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006, 2010, 2015, 2017
!        Salvatore Filippone    Cranfield University
!        Alfredo Buttari        CNRS-IRIT, Toulouse
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	Module to   define PREC_DATA,           !!
!!      structure for preconditioning.          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module psb_s_base_prec_mod

  ! Reduces size of .mod file.
  use psb_base_mod, only : psb_spk_, psb_ipk_, psb_long_int_k_,&
       & psb_desc_type, psb_sizeof, psb_free, psb_cdfree, psb_errpush, psb_act_abort_,&
       & psb_sizeof_int, psb_sizeof_long_int, psb_sizeof_sp, psb_sizeof_dp, &
       & psb_erractionsave, psb_erractionrestore, psb_error, psb_get_errstatus, psb_success_,&
       & psb_s_base_sparse_mat, psb_sspmat_type, psb_s_csr_sparse_mat,& 
       & psb_s_base_vect_type, psb_s_vect_type, psb_i_base_vect_type

  use psb_prec_const_mod

  type, abstract :: psb_s_base_prec_type
    integer(psb_ipk_) :: ictxt
  contains
    procedure, pass(prec) :: set_ctxt   => psb_s_base_set_ctxt
    procedure, pass(prec) :: get_ctxt   => psb_s_base_get_ctxt
    procedure, pass(prec) :: get_nzeros => psb_s_base_get_nzeros
    procedure, pass(prec) :: precseti   => psb_s_base_precseti
    procedure, pass(prec) :: precsetr   => psb_s_base_precsetr
    procedure, pass(prec) :: precsetc   => psb_s_base_precsetc
    generic, public       :: precset    => precseti, precsetr, precsetc
    procedure(psb_s_base_apply_vect), pass(prec), deferred :: s_apply_v  
    procedure(psb_s_base_apply), pass(prec), deferred :: s_apply    
    generic, public       :: apply     => s_apply, s_apply_v
    generic, public       :: build     => precbld
    generic, public       :: descr     => precdescr
    procedure(psb_s_base_precbld), pass(prec), deferred :: precbld    
    procedure(psb_s_base_sizeof), pass(prec), deferred :: sizeof     
    procedure(psb_s_base_precinit), pass(prec), deferred :: precinit   
    procedure(psb_s_base_precfree), pass(prec), deferred :: free   
    procedure(psb_s_base_precdescr), pass(prec), deferred :: precdescr  
    procedure(psb_s_base_precdump), pass(prec), deferred :: dump       
    procedure(psb_s_base_precclone), pass(prec), deferred :: clone       
  end type psb_s_base_prec_type

  private :: psb_s_base_set_ctxt, psb_s_base_get_ctxt, &
       & psb_s_base_get_nzeros

  abstract interface 
    subroutine psb_s_base_apply_vect(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import psb_ipk_, psb_spk_, psb_desc_type, psb_s_vect_type, &
           & psb_s_base_vect_type, psb_sspmat_type, psb_s_base_prec_type,&
           & psb_s_base_sparse_mat
      implicit none 
      type(psb_desc_type),intent(in)        :: desc_data
      class(psb_s_base_prec_type), intent(inout)  :: prec
      real(psb_spk_),intent(in)          :: alpha, beta
      type(psb_s_vect_type),intent(inout)   :: x
      type(psb_s_vect_type),intent(inout)   :: y
      integer(psb_ipk_), intent(out)                  :: info
      character(len=1), optional            :: trans
      real(psb_spk_),intent(inout), optional, target :: work(:)

    end subroutine psb_s_base_apply_vect
  end interface

  abstract interface 
    subroutine psb_s_base_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import psb_ipk_, psb_spk_, psb_desc_type, psb_s_vect_type, &
           & psb_s_base_vect_type, psb_sspmat_type, psb_s_base_prec_type,&
           & psb_s_base_sparse_mat
      implicit none 
      type(psb_desc_type),intent(in)       :: desc_data
      class(psb_s_base_prec_type), intent(inout)  :: prec
      real(psb_spk_),intent(in)         :: alpha, beta
      real(psb_spk_),intent(inout)      :: x(:)
      real(psb_spk_),intent(inout)      :: y(:)
      integer(psb_ipk_), intent(out)                 :: info
      character(len=1), optional           :: trans
      real(psb_spk_),intent(inout), optional, target :: work(:)

    end subroutine psb_s_base_apply
  end interface


  abstract interface 
    subroutine psb_s_base_precinit(prec,info)
      import psb_ipk_, psb_spk_, psb_desc_type, psb_s_vect_type, &
           & psb_s_base_vect_type, psb_sspmat_type, psb_s_base_prec_type,&
           & psb_s_base_sparse_mat
      Implicit None

      class(psb_s_base_prec_type),intent(inout) :: prec
      integer(psb_ipk_), intent(out)                     :: info

    end subroutine psb_s_base_precinit
  end interface


  abstract interface 
    subroutine psb_s_base_precbld(a,desc_a,prec,info,amold,vmold,imold)
      import psb_ipk_, psb_spk_, psb_desc_type, psb_s_vect_type, &
           & psb_s_base_vect_type, psb_sspmat_type, psb_s_base_prec_type,&
           & psb_s_base_sparse_mat, psb_i_base_vect_type
      Implicit None

      type(psb_sspmat_type), intent(in), target :: a
      type(psb_desc_type), intent(inout), target   :: desc_a
      class(psb_s_base_prec_type),intent(inout) :: prec
      integer(psb_ipk_), intent(out)                      :: info
      class(psb_s_base_sparse_mat), intent(in), optional :: amold
      class(psb_s_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine psb_s_base_precbld
  end interface


  abstract interface 
    subroutine psb_s_base_precfree(prec,info)
      import psb_ipk_, psb_spk_, psb_desc_type, psb_s_vect_type, &
           & psb_s_base_vect_type, psb_sspmat_type, psb_s_base_prec_type,&
           & psb_s_base_sparse_mat
      Implicit None

      class(psb_s_base_prec_type), intent(inout) :: prec
      integer(psb_ipk_), intent(out)                :: info

    end subroutine psb_s_base_precfree
  end interface


  abstract interface 
    subroutine psb_s_base_precdescr(prec,iout)
      import psb_ipk_, psb_spk_, psb_desc_type, psb_s_vect_type, &
           & psb_s_base_vect_type, psb_sspmat_type, psb_s_base_prec_type,&
           & psb_s_base_sparse_mat
      Implicit None

      class(psb_s_base_prec_type), intent(in) :: prec
      integer(psb_ipk_), intent(in), optional    :: iout

    end subroutine psb_s_base_precdescr
  end interface

  abstract interface   
    subroutine psb_s_base_precdump(prec,info,prefix,head)
      import psb_ipk_, psb_spk_, psb_desc_type, psb_s_vect_type, &
           & psb_s_base_vect_type, psb_sspmat_type, psb_s_base_prec_type,&
           & psb_s_base_sparse_mat
      implicit none 
      class(psb_s_base_prec_type), intent(in) :: prec
      integer(psb_ipk_), intent(out)             :: info
      character(len=*), intent(in), optional :: prefix,head

    end subroutine psb_s_base_precdump
  end interface
 
  abstract interface   
    subroutine psb_s_base_precclone(prec,precout,info)
      import psb_ipk_, psb_spk_, psb_desc_type, psb_s_vect_type, &
           & psb_s_base_vect_type, psb_sspmat_type, psb_s_base_prec_type,&
           & psb_s_base_sparse_mat
      implicit none 
      class(psb_s_base_prec_type), intent(inout)              :: prec
      class(psb_s_base_prec_type), allocatable, intent(inout) :: precout
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_s_base_precclone
  end interface

contains

  subroutine psb_s_base_precseti(prec,what,val,info)
    Implicit None

    class(psb_s_base_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(in)                      :: what 
    integer(psb_ipk_), intent(in)                      :: val 
    integer(psb_ipk_), intent(out)                     :: info
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='s_base_precseti'

    !
    ! Base version does nothing.

    info = psb_success_ 

    return

  end subroutine psb_s_base_precseti

  subroutine psb_s_base_precsetr(prec,what,val,info)
    Implicit None

    class(psb_s_base_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(in)                      :: what 
    real(psb_spk_), intent(in)               :: val 
    integer(psb_ipk_), intent(out)                     :: info
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='s_base_precsetr'


    !
    ! Base version does nothing.
    !

    info = psb_success_ 

    return

  end subroutine psb_s_base_precsetr

  subroutine psb_s_base_precsetc(prec,what,val,info)
    Implicit None

    class(psb_s_base_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(in)                      :: what 
    character(len=*), intent(in)             :: val
    integer(psb_ipk_), intent(out)                     :: info
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='s_base_precsetc'


    !
    ! Base version does nothing.
    !

    info = psb_success_ 

    return

  end subroutine psb_s_base_precsetc

  subroutine psb_s_base_set_ctxt(prec,ictxt)
    implicit none 
    class(psb_s_base_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in)  :: ictxt

    prec%ictxt = ictxt

  end subroutine psb_s_base_set_ctxt

  function psb_s_base_sizeof(prec) result(val)
    class(psb_s_base_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val

    val = 0
    return
  end function psb_s_base_sizeof

  function psb_s_base_get_ctxt(prec) result(val)
    class(psb_s_base_prec_type), intent(in) :: prec
    integer(psb_ipk_) :: val

    val = prec%ictxt
    return
  end function psb_s_base_get_ctxt

  function psb_s_base_get_nzeros(prec) result(res)
    class(psb_s_base_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: res

    res = 0

  end function psb_s_base_get_nzeros

end module psb_s_base_prec_mod
