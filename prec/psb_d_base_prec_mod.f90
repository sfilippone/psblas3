!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
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
module psb_d_base_prec_mod

  ! Reduces size of .mod file.
  use psb_base_mod, only : psb_dpk_, psb_ipk_, psb_epk_, psb_ctxt_type, &
       & psb_desc_type, psb_sizeof, psb_free, psb_cdfree, psb_errpush, psb_act_abort_,&
       & psb_sizeof_ip, psb_sizeof_lp, psb_sizeof_sp, psb_sizeof_dp, &
       & psb_erractionsave, psb_erractionrestore, psb_error, &
       & psb_errstatus_fatal, psb_success_,&
       & psb_d_base_sparse_mat, psb_dspmat_type, psb_d_csr_sparse_mat,& 
       & psb_d_base_vect_type, psb_d_vect_type, psb_i_base_vect_type

  use psb_prec_const_mod

  type, abstract :: psb_d_base_prec_type
    type(psb_ctxt_type) :: ctxt
  contains
    procedure, pass(prec) :: set_ctxt   => psb_d_base_set_ctxt
    procedure, pass(prec) :: get_ctxt   => psb_d_base_get_ctxt
    procedure, pass(prec) :: get_nzeros => psb_d_base_get_nzeros
    procedure, pass(prec) :: precseti   => psb_d_base_precseti
    procedure, pass(prec) :: precsetr   => psb_d_base_precsetr
    procedure, pass(prec) :: precsetc   => psb_d_base_precsetc
    generic, public       :: precset    => precseti, precsetr, precsetc
    procedure(psb_d_base_apply_vect), pass(prec), deferred :: d_apply_v  
    procedure(psb_d_base_apply), pass(prec), deferred :: d_apply    
    generic, public       :: apply     => d_apply, d_apply_v
    generic, public       :: build     => precbld
    generic, public       :: descr     => precdescr
    procedure, pass(prec) :: desc_prefix => psb_d_base_desc_prefix
    procedure, pass(prec) :: allocate_wrk => psb_d_base_allocate_wrk
    procedure, pass(prec) :: free_wrk     => psb_d_base_free_wrk
    procedure, pass(prec) :: is_allocated_wrk => psb_d_base_is_allocated_wrk
    procedure(psb_d_base_precbld), pass(prec), deferred :: precbld    
    procedure(psb_d_base_sizeof), pass(prec), deferred :: sizeof     
    procedure(psb_d_base_precinit), pass(prec), deferred :: precinit   
    procedure(psb_d_base_precfree), pass(prec), deferred :: free   
    procedure(psb_d_base_precdescr), pass(prec), deferred :: precdescr  
    procedure(psb_d_base_precdump), pass(prec), deferred :: dump       
    procedure(psb_d_base_precclone), pass(prec), deferred :: clone       
  end type psb_d_base_prec_type

  private :: psb_d_base_set_ctxt, psb_d_base_get_ctxt, &
       & psb_d_base_get_nzeros

  abstract interface 
    subroutine psb_d_base_apply_vect(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import psb_ipk_, psb_dpk_, psb_desc_type, psb_d_vect_type, &
           & psb_d_base_vect_type, psb_dspmat_type, psb_d_base_prec_type,&
           & psb_d_base_sparse_mat
      implicit none 
      type(psb_desc_type),intent(in)        :: desc_data
      class(psb_d_base_prec_type), intent(inout)  :: prec
      real(psb_dpk_),intent(in)          :: alpha, beta
      type(psb_d_vect_type),intent(inout)   :: x
      type(psb_d_vect_type),intent(inout)   :: y
      integer(psb_ipk_), intent(out)                  :: info
      character(len=1), optional            :: trans
      real(psb_dpk_),intent(inout), optional, target :: work(:)

    end subroutine psb_d_base_apply_vect
  end interface

  abstract interface 
    subroutine psb_d_base_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import psb_ipk_, psb_dpk_, psb_desc_type, psb_d_vect_type, &
           & psb_d_base_vect_type, psb_dspmat_type, psb_d_base_prec_type,&
           & psb_d_base_sparse_mat
      implicit none 
      type(psb_desc_type),intent(in)       :: desc_data
      class(psb_d_base_prec_type), intent(inout)  :: prec
      real(psb_dpk_),intent(in)         :: alpha, beta
      real(psb_dpk_),intent(inout)      :: x(:)
      real(psb_dpk_),intent(inout)      :: y(:)
      integer(psb_ipk_), intent(out)                 :: info
      character(len=1), optional           :: trans
      real(psb_dpk_),intent(inout), optional, target :: work(:)

    end subroutine psb_d_base_apply
  end interface


  abstract interface 
    subroutine psb_d_base_precinit(prec,info)
      import psb_ipk_, psb_dpk_, psb_desc_type, psb_d_vect_type, &
           & psb_d_base_vect_type, psb_dspmat_type, psb_d_base_prec_type,&
           & psb_d_base_sparse_mat
      Implicit None

      class(psb_d_base_prec_type),intent(inout) :: prec
      integer(psb_ipk_), intent(out)                     :: info

    end subroutine psb_d_base_precinit
  end interface


  abstract interface 
    subroutine psb_d_base_precbld(a,desc_a,prec,info,amold,vmold,imold)
      import psb_ipk_, psb_dpk_, psb_desc_type, psb_d_vect_type, &
           & psb_d_base_vect_type, psb_dspmat_type, psb_d_base_prec_type,&
           & psb_d_base_sparse_mat, psb_i_base_vect_type
      Implicit None

      type(psb_dspmat_type), intent(in), target :: a
      type(psb_desc_type), intent(inout), target   :: desc_a
      class(psb_d_base_prec_type),intent(inout) :: prec
      integer(psb_ipk_), intent(out)                      :: info
      class(psb_d_base_sparse_mat), intent(in), optional :: amold
      class(psb_d_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine psb_d_base_precbld
  end interface


  abstract interface 
    subroutine psb_d_base_precfree(prec,info)
      import psb_ipk_, psb_dpk_, psb_desc_type, psb_d_vect_type, &
           & psb_d_base_vect_type, psb_dspmat_type, psb_d_base_prec_type,&
           & psb_d_base_sparse_mat
      Implicit None

      class(psb_d_base_prec_type), intent(inout) :: prec
      integer(psb_ipk_), intent(out)                :: info

    end subroutine psb_d_base_precfree
  end interface


  abstract interface 
    subroutine psb_d_base_precdescr(prec,iout,root, verbosity,prefix)
      import psb_ipk_, psb_dpk_, psb_desc_type, psb_d_vect_type, &
           & psb_d_base_vect_type, psb_dspmat_type, psb_d_base_prec_type,&
           & psb_d_base_sparse_mat
      Implicit None

      class(psb_d_base_prec_type), intent(in) :: prec
      integer(psb_ipk_), intent(in), optional   :: iout
      integer(psb_ipk_), intent(in), optional   :: root
      integer(psb_ipk_), intent(in), optional   :: verbosity
      character(len=*), intent(in), optional    :: prefix
    end subroutine psb_d_base_precdescr
  end interface

  abstract interface   
    subroutine psb_d_base_precdump(prec,info,prefix,head)
      import psb_ipk_, psb_dpk_, psb_desc_type, psb_d_vect_type, &
           & psb_d_base_vect_type, psb_dspmat_type, psb_d_base_prec_type,&
           & psb_d_base_sparse_mat
      implicit none 
      class(psb_d_base_prec_type), intent(in) :: prec
      integer(psb_ipk_), intent(out)             :: info
      character(len=*), intent(in), optional :: prefix,head

    end subroutine psb_d_base_precdump
  end interface
 
  abstract interface   
    subroutine psb_d_base_precclone(prec,precout,info)
      import psb_ipk_, psb_dpk_, psb_desc_type, psb_d_vect_type, &
           & psb_d_base_vect_type, psb_dspmat_type, psb_d_base_prec_type,&
           & psb_d_base_sparse_mat
      implicit none 
      class(psb_d_base_prec_type), intent(inout)              :: prec
      class(psb_d_base_prec_type), allocatable, intent(inout) :: precout
      integer(psb_ipk_), intent(out)               :: info
    end subroutine psb_d_base_precclone
  end interface

contains

  subroutine psb_d_base_precseti(prec,what,val,info)
    Implicit None

    class(psb_d_base_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(in)                      :: what 
    integer(psb_ipk_), intent(in)                      :: val 
    integer(psb_ipk_), intent(out)                     :: info
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='d_base_precseti'

    !
    ! Base version does nothing.

    info = psb_success_ 

    return

  end subroutine psb_d_base_precseti

  subroutine psb_d_base_precsetr(prec,what,val,info)
    Implicit None

    class(psb_d_base_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(in)                      :: what 
    real(psb_dpk_), intent(in)               :: val 
    integer(psb_ipk_), intent(out)                     :: info
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='d_base_precsetr'


    !
    ! Base version does nothing.
    !

    info = psb_success_ 

    return

  end subroutine psb_d_base_precsetr

  subroutine psb_d_base_precsetc(prec,what,val,info)
    Implicit None

    class(psb_d_base_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(in)                      :: what 
    character(len=*), intent(in)             :: val
    integer(psb_ipk_), intent(out)                     :: info
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='d_base_precsetc'


    !
    ! Base version does nothing.
    !

    info = psb_success_ 

    return

  end subroutine psb_d_base_precsetc

  subroutine psb_d_base_allocate_wrk(prec,info,vmold,desc)
    use psb_base_mod
    implicit none
    
    ! Arguments
    class(psb_d_base_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)        :: info
    class(psb_d_base_vect_type), intent(in), optional  :: vmold
    type(psb_desc_type), intent(in), optional :: desc

    ! Local variables
    integer(psb_ipk_) :: err_act
    character(len=20)   :: name
    
    info=psb_success_
    name = 'psb_d_allocate_wrk'
    call psb_erractionsave(err_act)
    
    if (psb_get_errstatus().ne.0) goto 9999

    !
    ! Base version does nothing.
    !

    info = psb_success_ 

    call psb_erractionrestore(err_act)
    return
    
9999 call psb_error_handler(err_act)
    return
    
  end subroutine psb_d_base_allocate_wrk

  subroutine psb_d_base_free_wrk(prec,info)
    use psb_base_mod
    implicit none
    
    ! Arguments
    class(psb_d_base_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)        :: info

    ! Local variables
    integer(psb_ipk_) :: err_act
    character(len=20)   :: name
    
    info=psb_success_
    name = 'psb_d_allocate_wrk'
    call psb_erractionsave(err_act)
    
    if (psb_get_errstatus().ne.0) goto 9999

    !
    ! Base version does nothing.
    !

    info = psb_success_ 

    call psb_erractionrestore(err_act)
    return
    
9999 call psb_error_handler(err_act)
    return
    
  end subroutine psb_d_base_free_wrk

  function psb_d_base_is_allocated_wrk(prec) result(res)
    use psb_base_mod
    implicit none
    
    ! Arguments
    class(psb_d_base_prec_type), intent(in) :: prec
    logical :: res

    ! In the base version we can say yes, because 
    ! there is nothing to allocate

    res = .true.
    
  end function psb_d_base_is_allocated_wrk
  
  subroutine psb_d_base_set_ctxt(prec,ctxt)
    implicit none 
    class(psb_d_base_prec_type), intent(inout) :: prec
    type(psb_ctxt_type) :: ctxt

    prec%ctxt = ctxt

  end subroutine psb_d_base_set_ctxt

  function psb_d_base_sizeof(prec) result(val)
    class(psb_d_base_prec_type), intent(in) :: prec
    integer(psb_epk_) :: val

    val = 0
    return
  end function psb_d_base_sizeof

  function psb_d_base_get_ctxt(prec) result(val)
    class(psb_d_base_prec_type), intent(in) :: prec
    type(psb_ctxt_type) :: val

    val = prec%ctxt
    return
  end function psb_d_base_get_ctxt

  function psb_d_base_get_nzeros(prec) result(res)
    class(psb_d_base_prec_type), intent(in) :: prec
    integer(psb_epk_) :: res

    res = 0

  end function psb_d_base_get_nzeros

  function psb_d_base_desc_prefix(prec) result(res)
    use psb_base_mod, only : psb_info, psb_root_
    implicit none 
    class(psb_d_base_prec_type), intent(in) :: prec
    character(len=32) :: res 
    !
    character(len=32) :: frmtv
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: ni, iam, np

    ctxt = prec%ctxt
    call psb_info(ctxt,iam,np)
    
    res = ''
    if (iam /= psb_root_) then
      ni  = floor(log10(1.0*np)) + 1
      write(frmtv,'(a,i8.8,a)') '(a,i',ni,',a)'
      write(res,frmtv) 'Process ',iam,': Preconditioner: '
    else
      write(res,'(a)') 'Preconditioner: '
    end if

  end function psb_d_base_desc_prefix

end module psb_d_base_prec_mod
