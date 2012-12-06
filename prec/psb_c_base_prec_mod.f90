!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	Module to   define PREC_DATA,           !!
!!      structure for preconditioning.          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module psb_c_base_prec_mod

  ! Reduces size of .mod file.
  use psb_base_mod, only : psb_spk_, psb_ipk_, psb_long_int_k_,&
       & psb_desc_type, psb_sizeof, psb_free, psb_cdfree, psb_errpush, psb_act_abort_,&
       & psb_sizeof_int, psb_sizeof_long_int, psb_sizeof_sp, psb_sizeof_dp, &
       & psb_erractionsave, psb_erractionrestore, psb_error, psb_get_errstatus, psb_success_,&
       & psb_c_base_sparse_mat, psb_cspmat_type, psb_c_csr_sparse_mat,& 
       & psb_c_base_vect_type, psb_c_vect_type

  use psb_prec_const_mod

  type psb_c_base_prec_type
    integer(psb_ipk_) :: ictxt
  contains
    procedure, pass(prec) :: set_ctxt  => psb_c_base_set_ctxt
    procedure, pass(prec) :: get_ctxt  => psb_c_base_get_ctxt
    procedure, pass(prec) :: c_apply_v => psb_c_base_apply_vect
    procedure, pass(prec) :: c_apply   => psb_c_base_apply
    generic, public       :: apply     => c_apply, c_apply_v
    procedure, pass(prec) :: precbld   => psb_c_base_precbld
    procedure, pass(prec) :: precseti  => psb_c_base_precseti
    procedure, pass(prec) :: precsetr  => psb_c_base_precsetr
    procedure, pass(prec) :: precsetc  => psb_c_base_precsetc
    procedure, pass(prec) :: sizeof    => psb_c_base_sizeof
    generic, public       :: precset   => precseti, precsetr, precsetc
    procedure, pass(prec) :: precinit  => psb_c_base_precinit
    procedure, pass(prec) :: precfree  => psb_c_base_precfree
    procedure, pass(prec) :: precdescr => psb_c_base_precdescr
    procedure, pass(prec) :: dump      => psb_c_base_precdump
    procedure, pass(prec) :: get_nzeros => psb_c_base_get_nzeros
  end type psb_c_base_prec_type
  
  private :: psb_c_base_apply, psb_c_base_precbld, psb_c_base_precseti,&
       & psb_c_base_precsetr, psb_c_base_precsetc, psb_c_base_sizeof,&
       & psb_c_base_precinit, psb_c_base_precfree, psb_c_base_precdescr,&
       & psb_c_base_precdump, psb_c_base_set_ctxt, psb_c_base_get_ctxt, &
       & psb_c_base_apply_vect, psb_c_base_get_nzeros
  
contains

  subroutine psb_c_base_apply_vect(alpha,prec,x,beta,y,desc_data,info,trans,work)
    implicit none 
    type(psb_desc_type),intent(in)        :: desc_data
    class(psb_c_base_prec_type), intent(inout)  :: prec
    complex(psb_spk_),intent(in)          :: alpha, beta
    type(psb_c_vect_type),intent(inout)   :: x
    type(psb_c_vect_type),intent(inout)   :: y
    integer(psb_ipk_), intent(out)                  :: info
    character(len=1), optional            :: trans
    complex(psb_spk_),intent(inout), optional, target :: work(:)
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='c_base_prec_apply'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine psb_c_base_apply_vect

  subroutine psb_c_base_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
    implicit none 
    type(psb_desc_type),intent(in)       :: desc_data
    class(psb_c_base_prec_type), intent(in)  :: prec
    complex(psb_spk_),intent(in)         :: alpha, beta
    complex(psb_spk_),intent(inout)      :: x(:)
    complex(psb_spk_),intent(inout)      :: y(:)
    integer(psb_ipk_), intent(out)                 :: info
    character(len=1), optional           :: trans
    complex(psb_spk_),intent(inout), optional, target :: work(:)
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='c_base_prec_apply'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine psb_c_base_apply

  subroutine psb_c_base_precinit(prec,info)
    Implicit None
   
    class(psb_c_base_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(out)                     :: info
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='c_base_precinit'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_c_base_precinit

  subroutine psb_c_base_precbld(a,desc_a,prec,info,upd,amold,afmt,vmold)
    Implicit None
    
    type(psb_cspmat_type), intent(in), target :: a
    type(psb_desc_type), intent(in), target   :: desc_a
    class(psb_c_base_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(out)                      :: info
    character, intent(in), optional           :: upd
    character(len=*), intent(in), optional    :: afmt
    class(psb_c_base_sparse_mat), intent(in), optional :: amold
    class(psb_c_base_vect_type), intent(in), optional  :: vmold
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='c_base_precbld'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_c_base_precbld

  subroutine psb_c_base_precseti(prec,what,val,info)
    Implicit None
    
    class(psb_c_base_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(in)                      :: what 
    integer(psb_ipk_), intent(in)                      :: val 
    integer(psb_ipk_), intent(out)                     :: info
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='c_base_precseti'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_c_base_precseti

  subroutine psb_c_base_precsetr(prec,what,val,info)
    Implicit None
    
    class(psb_c_base_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(in)                      :: what 
    real(psb_spk_), intent(in)               :: val 
    integer(psb_ipk_), intent(out)                     :: info
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='c_base_precsetr'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_c_base_precsetr

  subroutine psb_c_base_precsetc(prec,what,val,info)
    Implicit None
    
    class(psb_c_base_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(in)                      :: what 
    character(len=*), intent(in)             :: val
    integer(psb_ipk_), intent(out)                     :: info
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='c_base_precsetc'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_c_base_precsetc

  subroutine psb_c_base_precfree(prec,info)
    Implicit None

    class(psb_c_base_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)                :: info
    
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='c_base_precfree'
    
    call psb_erractionsave(err_act)
    
    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine psb_c_base_precfree
  

  subroutine psb_c_base_precdescr(prec,iout)
    Implicit None

    class(psb_c_base_prec_type), intent(in) :: prec
    integer(psb_ipk_), intent(in), optional    :: iout

    integer(psb_ipk_) :: err_act, nrow, info
    character(len=20)  :: name='c_base_precdescr'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine psb_c_base_precdescr
  
  subroutine psb_c_base_precdump(prec,info,prefix,head)
    implicit none 
    class(psb_c_base_prec_type), intent(in) :: prec
    integer(psb_ipk_), intent(out)             :: info
    character(len=*), intent(in), optional :: prefix,head
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='c_base_precdump'

    call psb_erractionsave(err_act)

    !
    ! This is the base version and we should throw an error. 
    ! Or should it be the NULL preonditioner???
    !
    info = 700
    call psb_errpush(info,name)
    goto 9999 
    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine psb_c_base_precdump

  subroutine psb_c_base_set_ctxt(prec,ictxt)
    implicit none 
    class(psb_c_base_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(in)  :: ictxt

    prec%ictxt = ictxt

  end subroutine psb_c_base_set_ctxt

  function psb_c_base_sizeof(prec) result(val)
    class(psb_c_base_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    
    val = 0
    return
  end function psb_c_base_sizeof

  function psb_c_base_get_ctxt(prec) result(val)
    class(psb_c_base_prec_type), intent(in) :: prec
    integer(psb_ipk_) :: val
    
    val = prec%ictxt
    return
  end function psb_c_base_get_ctxt

  function psb_c_base_get_nzeros(prec) result(res)
    class(psb_c_base_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: res
    
    res = 0

  end function psb_c_base_get_nzeros

end module psb_c_base_prec_mod
