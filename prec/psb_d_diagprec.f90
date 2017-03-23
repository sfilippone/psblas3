!   
!                Parallel Sparse BLAS  version 3.4
!      (C) Copyright 2006, 2010, 2015
!                         Salvatore Filippone    University of Rome Tor Vergata
!                         Alfredo Buttari        CNRS-IRIT, Toulouse
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
module psb_d_diagprec

  use psb_d_base_prec_mod
  
  type, extends(psb_d_base_prec_type) :: psb_d_diag_prec_type
    real(psb_dpk_), allocatable     :: d(:)
    type(psb_d_vect_type), allocatable :: dv
  contains
    procedure, pass(prec) :: d_apply_v => psb_d_diag_apply_vect
    procedure, pass(prec) :: d_apply   => psb_d_diag_apply
    procedure, pass(prec) :: precbld    => psb_d_diag_precbld
    procedure, pass(prec) :: precinit   => psb_d_diag_precinit  
    procedure, pass(prec) :: precdescr  => psb_d_diag_precdescr
    procedure, pass(prec) :: sizeof     => psb_d_diag_sizeof
    procedure, pass(prec) :: dump       => psb_d_diag_dump
    procedure, pass(prec) :: clone      => psb_d_diag_clone
    procedure, pass(prec) :: free       => psb_d_diag_precfree
    procedure, pass(prec) :: get_nzeros => psb_d_diag_get_nzeros
  end type psb_d_diag_prec_type

  private :: psb_d_diag_sizeof,&
       & psb_d_diag_precinit, psb_d_diag_precfree, psb_d_diag_precdescr,&
       & psb_d_diag_get_nzeros
  
  
  
  interface  
    subroutine psb_d_diag_apply_vect(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import :: psb_ipk_, psb_desc_type, psb_d_diag_prec_type, psb_d_vect_type, psb_dpk_
      type(psb_desc_type),intent(in)    :: desc_data
      class(psb_d_diag_prec_type), intent(inout)  :: prec
      type(psb_d_vect_type),intent(inout)   :: x
      real(psb_dpk_),intent(in)         :: alpha, beta
      type(psb_d_vect_type),intent(inout)   :: y
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), optional        :: trans
      real(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine psb_d_diag_apply_vect
  end interface
  
  interface  
    subroutine psb_d_diag_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import :: psb_ipk_, psb_desc_type, psb_d_diag_prec_type, psb_d_vect_type, psb_dpk_
      type(psb_desc_type),intent(in)    :: desc_data
      class(psb_d_diag_prec_type), intent(inout)  :: prec
      real(psb_dpk_),intent(inout)      :: x(:)
      real(psb_dpk_),intent(in)         :: alpha, beta
      real(psb_dpk_),intent(inout)      :: y(:)
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), optional        :: trans
      real(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine psb_d_diag_apply
  end interface
  
  interface
    subroutine psb_d_diag_precbld(a,desc_a,prec,info,upd,amold,afmt,vmold)
      import :: psb_ipk_, psb_desc_type, psb_d_diag_prec_type, psb_d_vect_type, psb_dpk_, &
           & psb_dspmat_type, psb_d_base_sparse_mat, psb_d_base_vect_type
      type(psb_dspmat_type), intent(in), target :: a
      type(psb_desc_type), intent(in), target   :: desc_a
      class(psb_d_diag_prec_type),intent(inout) :: prec
      integer(psb_ipk_), intent(out)                      :: info
      character, intent(in), optional           :: upd
      character(len=*), intent(in), optional    :: afmt
      class(psb_d_base_sparse_mat), intent(in), optional :: amold
      class(psb_d_base_vect_type), intent(in), optional  :: vmold
    end subroutine psb_d_diag_precbld
  end interface

  interface 
    subroutine psb_d_diag_dump(prec,info,prefix,head)
      import :: psb_ipk_, psb_desc_type, psb_d_diag_prec_type, psb_d_vect_type, psb_dpk_
      implicit none 
      class(psb_d_diag_prec_type), intent(in) :: prec
      integer(psb_ipk_), intent(out)                    :: info
      character(len=*), intent(in), optional  :: prefix,head
    end subroutine psb_d_diag_dump
  end interface
  

contains
  

  subroutine psb_d_diag_precinit(prec,info)
    Implicit None

    class(psb_d_diag_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(out)                     :: info
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='d_diag_precinit'

    call psb_erractionsave(err_act)

    info = psb_success_


    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  end subroutine psb_d_diag_precinit


  subroutine psb_d_diag_precfree(prec,info)

    Implicit None

    class(psb_d_diag_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)                :: info

    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='d_diag_precset'

    call psb_erractionsave(err_act)

    info = psb_success_

    if (allocated(prec%dv)) call prec%dv%free(info)

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_d_diag_precfree
  

  subroutine psb_d_diag_precdescr(prec,iout)
    Implicit None

    class(psb_d_diag_prec_type), intent(in) :: prec
    integer(psb_ipk_), intent(in), optional    :: iout

    integer(psb_ipk_) :: err_act, nrow, info
    character(len=20)  :: name='d_diag_precdescr'

    integer(psb_ipk_) :: iout_

    call psb_erractionsave(err_act)

    info = psb_success_

    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if

    write(iout_,*) 'Diagonal scaling'

    call psb_erractionsave(err_act)

    info = psb_success_

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_d_diag_precdescr

  function psb_d_diag_sizeof(prec) result(val)
    class(psb_d_diag_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    
    val = psb_sizeof_dp * prec%get_nzeros()
    return
  end function psb_d_diag_sizeof

  function psb_d_diag_get_nzeros(prec) result(val)
    class(psb_d_diag_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val

    val = 0
    if (allocated(prec%dv)) val = val + prec%dv%get_nrows()
    return
  end function psb_d_diag_get_nzeros


  subroutine psb_d_diag_clone(prec,precout,info)
    use psb_error_mod
    use psb_realloc_mod

    Implicit None

    class(psb_d_diag_prec_type), intent(inout) :: prec
    class(psb_d_base_prec_type), allocatable, intent(inout)  :: precout
    integer(psb_ipk_), intent(out)               :: info

    integer(psb_ipk_) :: err_act, i
    character(len=20)  :: name='d_diag_clone'

    call psb_erractionsave(err_act)

    info = psb_success_
    if (allocated(precout)) then
      call precout%free(info)
      if (info == psb_success_) deallocate(precout, stat=info)
    end if
    if (info == psb_success_) &
         & allocate(psb_d_diag_prec_type :: precout, stat=info)
    if (info /= 0) goto 9999
    select type(pout => precout)
    type is (psb_d_diag_prec_type) 
      call pout%set_ctxt(prec%get_ctxt())

      if (allocated(prec%dv)) then 
        allocate(pout%dv,stat=info)
        if (info == 0) call prec%dv%clone(pout%dv,info)
      end if
      if (info == 0) call psb_safe_ab_cpy(prec%d,pout%d,info)
      class default
      info = psb_err_internal_error_
    end select
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_d_diag_clone


end module psb_d_diagprec
