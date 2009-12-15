!!$ 
!!$              Parallel Sparse BLAS  version 2.2
!!$    (C) Copyright 2006/2007/2008
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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
module psb_s_prec_type

  ! Reduces size of .mod file.
  use psb_sparse_mod, only : psb_dpk_, psb_spk_, psb_long_int_k_,&
       & psb_desc_type, psb_sizeof, psb_free, psb_cdfree,&
       & psb_erractionsave, psb_erractionrestore, psb_error, psb_get_errstatus,&
       & psb_s_sparse_mat
  
  use psb_prec_const_mod


  type psb_s_base_prec_type
  contains
    procedure, pass(prec) :: apply     => s_base_apply
    procedure, pass(prec) :: precbld   => s_base_precbld
    procedure, pass(prec) :: precseti  => s_base_precseti
    procedure, pass(prec) :: precsetr  => s_base_precsetr
    procedure, pass(prec) :: precsetc  => s_base_precsetc
    procedure, pass(prec) :: sizeof    => s_base_sizeof
    generic, public       :: precset   => precseti, precsetr, precsetc
    procedure, pass(prec) :: precinit  => s_base_precinit
    procedure, pass(prec) :: precfree  => s_base_precfree
    procedure, pass(prec) :: precdescr => s_base_precdescr
  end type psb_s_base_prec_type
  
  type psb_sprec_type
    class(psb_s_base_prec_type), allocatable :: prec
  contains
    procedure, pass(prec)               :: s_apply2v
    procedure, pass(prec)               :: s_apply1v
    generic, public                     :: apply => s_apply2v, s_apply1v
  end type psb_sprec_type


  interface psb_precfree
    module procedure psb_s_precfree
  end interface

  interface psb_nullify_prec
    module procedure psb_nullify_sprec
  end interface

  interface psb_precdescr
    module procedure  psb_sfile_prec_descr
  end interface

  interface psb_sizeof
    module procedure psb_sprec_sizeof
  end interface



  interface psb_precaply
    subroutine psb_sprc_aply(prec,x,y,desc_data,info,trans,work)
      use psb_sparse_mod, only  : psb_desc_type, psb_spk_
      import psb_sprec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_sprec_type), intent(in)  :: prec
      real(psb_spk_),intent(in)         :: x(:)
      real(psb_spk_),intent(inout)      :: y(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
      real(psb_spk_),intent(inout), optional, target :: work(:)
    end subroutine psb_sprc_aply
    subroutine psb_sprc_aply1(prec,x,desc_data,info,trans)
      use psb_sparse_mod, only  : psb_desc_type, psb_spk_
      import psb_sprec_type
      type(psb_desc_type),intent(in)    :: desc_data
      type(psb_sprec_type), intent(in)  :: prec
      real(psb_spk_),intent(inout)      :: x(:)
      integer, intent(out)              :: info
      character(len=1), optional        :: trans
    end subroutine psb_sprc_aply1
  end interface


contains

  
  subroutine psb_sfile_prec_descr(p,iout)
    use psb_sparse_mod
    type(psb_sprec_type), intent(in) :: p
    integer, intent(in), optional    :: iout
    integer :: iout_,info
    character(len=20) :: name='prec_descr' 
    
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if

    if (.not.allocated(p%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
    end if
    call p%prec%precdescr(iout)
    
  end subroutine psb_sfile_prec_descr


  subroutine psb_s_precfree(p,info)
    use psb_sparse_mod
    type(psb_sprec_type), intent(inout) :: p 
    integer, intent(out) ::  info
    integer :: me, err_act,i 
    character(len=20) :: name
    if(psb_get_errstatus() /= 0) return
    info=0
    name = 'psb_precfree'
    call psb_erractionsave(err_act)

    me=-1

    if (allocated(p%prec)) then 
      call p%prec%precfree(info)
      if (info /= 0) goto 9999
      deallocate(p%prec,stat=info)
      if (info /= 0) goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine psb_s_precfree

  subroutine psb_nullify_sprec(p)
    type(psb_sprec_type), intent(inout) :: p

  end subroutine psb_nullify_sprec

  function psb_sprec_sizeof(prec) result(val)
    use psb_sparse_mod
    type(psb_sprec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    integer             :: i
    
    val = 0
    
    if (allocated(prec%prec)) then 
      val = val + prec%prec%sizeof()
    end if
    
  end function psb_sprec_sizeof
 
  subroutine s_apply2v(prec,x,y,desc_data,info,trans,work)
    use psb_sparse_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_sprec_type), intent(in)  :: prec
    real(psb_spk_),intent(in)       :: x(:)
    real(psb_spk_),intent(inout)    :: y(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    real(psb_spk_),intent(inout), optional, target :: work(:)
    
    character     :: trans_ 
    real(psb_spk_), pointer :: work_(:)
    integer :: ictxt,np,me,err_act
    character(len=20)   :: name
    
    name='s_apply2v'
    info = 0
    call psb_erractionsave(err_act)
    
    ictxt = psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)
    
    if (present(trans)) then 
      trans_=trans
    else
      trans_='N'
    end if
    
    if (present(work)) then 
      work_ => work
    else
      allocate(work_(4*psb_cd_get_local_cols(desc_data)),stat=info)
      if (info /= 0) then 
        info = 4010
        call psb_errpush(info,name,a_err='Allocate')
        goto 9999      
      end if
      
    end if
    
    if (.not.allocated(prec%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if
    call prec%prec%apply(sone,x,szero,y,desc_data,info,trans_,work=work_)
    if (present(work)) then 
    else
      deallocate(work_,stat=info)
      if (info /= 0) then 
        info = 4010
        call psb_errpush(info,name,a_err='DeAllocate')
        goto 9999      
      end if
    end if
    
    call psb_erractionrestore(err_act)
    return
    
9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine s_apply2v

  subroutine s_apply1v(prec,x,desc_data,info,trans)
    use psb_sparse_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_sprec_type), intent(in)  :: prec
    real(psb_spk_),intent(inout)    :: x(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans

    character     :: trans_
    integer :: ictxt,np,me, err_act
    real(psb_spk_), pointer :: WW(:), w1(:)
    character(len=20)   :: name
    name='s_apply1v'
    info = 0
    call psb_erractionsave(err_act)
    
    
    ictxt=psb_cd_get_context(desc_data)
    call psb_info(ictxt, me, np)
    if (present(trans)) then 
      trans_=psb_toupper(trans)
    else
      trans_='N'
    end if
    
    if (.not.allocated(prec%prec)) then 
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if
    allocate(ww(size(x)),w1(size(x)),stat=info)
    if (info /= 0) then 
      info = 4010
      call psb_errpush(info,name,a_err='Allocate')
      goto 9999      
    end if
    call prec%prec%apply(sone,x,szero,ww,desc_data,info,trans_,work=w1)
    if(info /=0) goto 9999
    x(:) = ww(:)
    deallocate(ww,W1,stat=info)
    if (info /= 0) then 
      info = 4010
      call psb_errpush(info,name,a_err='DeAllocate')
      goto 9999      
    end if

    
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_errpush(info,name)
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return
    
  end subroutine s_apply1v



  subroutine s_base_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
    use psb_sparse_mod
    type(psb_desc_type),intent(in)    :: desc_data
    class(psb_s_base_prec_type), intent(in)  :: prec
    real(psb_spk_),intent(in)         :: alpha, beta
    real(psb_spk_),intent(in)         :: x(:)
    real(psb_spk_),intent(inout)      :: y(:)
    integer, intent(out)              :: info
    character(len=1), optional        :: trans
    real(psb_spk_),intent(inout), optional, target :: work(:)
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_prec_apply'

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

  end subroutine s_base_apply

  subroutine s_base_precinit(prec,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_s_base_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_precinit'

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
  end subroutine s_base_precinit

  subroutine s_base_precbld(a,desc_a,prec,info,upd)
    
    use psb_sparse_mod
    Implicit None
    
    type(psb_s_sparse_mat), intent(in), target :: a
    type(psb_desc_type), intent(in), target  :: desc_a
    class(psb_s_base_prec_type),intent(inout) :: prec
    integer, intent(out)                     :: info
    character, intent(in), optional          :: upd
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_precbld'

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
  end subroutine s_base_precbld

  subroutine s_base_precseti(prec,what,val,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_s_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    integer, intent(in)                      :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_precseti'

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
  end subroutine s_base_precseti

  subroutine s_base_precsetr(prec,what,val,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_s_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    real(psb_spk_), intent(in)               :: val 
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_precsetr'

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
  end subroutine s_base_precsetr

  subroutine s_base_precsetc(prec,what,val,info)
    
    use psb_sparse_mod
    Implicit None
    
    class(psb_s_base_prec_type),intent(inout) :: prec
    integer, intent(in)                      :: what 
    character(len=*), intent(in)             :: val
    integer, intent(out)                     :: info
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_precsetc'

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
  end subroutine s_base_precsetc

  subroutine s_base_precfree(prec,info)
    
    use psb_sparse_mod
    Implicit None

    class(psb_s_base_prec_type), intent(inout) :: prec
    integer, intent(out)                :: info
    
    Integer :: err_act, nrow
    character(len=20)  :: name='s_base_precfree'
    
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
    
  end subroutine s_base_precfree
  

  subroutine s_base_precdescr(prec,iout)
    
    use psb_sparse_mod
    Implicit None

    class(psb_s_base_prec_type), intent(in) :: prec
    integer, intent(in), optional    :: iout

    Integer :: err_act, nrow, info
    character(len=20)  :: name='s_base_precdescr'

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
    
  end subroutine s_base_precdescr
  

  function s_base_sizeof(prec) result(val)
    use psb_sparse_mod
    class(psb_s_base_prec_type), intent(in) :: prec
    integer(psb_long_int_k_) :: val
    
    val = 0
    return
  end function s_base_sizeof


end module psb_s_prec_type
