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
module psb_z_nullprec

  use psb_z_base_prec_mod
  
  type, extends(psb_z_base_prec_type) :: psb_z_null_prec_type
  contains
    procedure, pass(prec) :: z_apply_v => psb_z_null_apply_vect
    procedure, pass(prec) :: z_apply   => psb_z_null_apply
    procedure, pass(prec) :: precbld   => psb_z_null_precbld
    procedure, pass(prec) :: precinit  => psb_z_null_precinit
    procedure, pass(prec) :: precdescr => psb_z_null_precdescr
    procedure, pass(prec) :: sizeof    => psb_z_null_sizeof
    procedure, pass(prec) :: dump      => psb_z_null_dump
    procedure, pass(prec) :: clone     => psb_z_null_clone
    procedure, pass(prec) :: free      => psb_z_null_precfree
  end type psb_z_null_prec_type

  private :: psb_z_null_precbld, psb_z_null_sizeof,&
       & psb_z_null_precinit, psb_z_null_precfree, psb_z_null_precdescr
  

  interface
    subroutine psb_z_null_apply_vect(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import :: psb_ipk_, psb_desc_type, psb_z_null_prec_type, psb_z_vect_type, psb_dpk_
      type(psb_desc_type),intent(in)       :: desc_data
      class(psb_z_null_prec_type), intent(inout)  :: prec
      type(psb_z_vect_type),intent(inout)  :: x
      complex(psb_dpk_),intent(in)         :: alpha, beta
      type(psb_z_vect_type),intent(inout)  :: y
      integer(psb_ipk_), intent(out)                 :: info
      character(len=1), optional           :: trans
      complex(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine psb_z_null_apply_vect
  end interface
  
  interface
    subroutine psb_z_null_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import :: psb_ipk_, psb_desc_type, psb_z_null_prec_type, psb_dpk_
      type(psb_desc_type),intent(in)       :: desc_data
      class(psb_z_null_prec_type), intent(inout)  :: prec
      complex(psb_dpk_),intent(inout)      :: x(:)
      complex(psb_dpk_),intent(in)         :: alpha, beta
      complex(psb_dpk_),intent(inout)      :: y(:)
      integer(psb_ipk_), intent(out)                 :: info
      character(len=1), optional           :: trans
      complex(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine psb_z_null_apply
  end interface
  
  
contains
  

  subroutine psb_z_null_precinit(prec,info)
    
    Implicit None
    
    class(psb_z_null_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(out)                     :: info
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='z_null_precinit'

    call psb_erractionsave(err_act)

    info = psb_success_

    
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  end subroutine psb_z_null_precinit

  subroutine psb_z_null_precbld(a,desc_a,prec,info,amold,vmold,imold)
    
    Implicit None
    
    type(psb_zspmat_type), intent(in), target :: a
    type(psb_desc_type), intent(inout), target   :: desc_a
    class(psb_z_null_prec_type),intent(inout) :: prec
    integer(psb_ipk_), intent(out)                      :: info
    class(psb_z_base_sparse_mat), intent(in), optional :: amold
    class(psb_z_base_vect_type), intent(in), optional  :: vmold
    class(psb_i_base_vect_type), intent(in), optional  :: imold
    
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='z_null_precbld'

    call psb_erractionsave(err_act)

    info = psb_success_
    call prec%set_ctxt(desc_a%get_ctxt())
    
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    
    return
  end subroutine psb_z_null_precbld

  subroutine psb_z_null_precfree(prec,info)
    
    Implicit None

    class(psb_z_null_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)                :: info
    
    integer(psb_ipk_) :: err_act, nrow
    character(len=20)  :: name='z_null_precset'
    
    call psb_erractionsave(err_act)
    
    info = psb_success_
    
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
    
  end subroutine psb_z_null_precfree
  

  !
  !
  ! verbosity:
  !        <0: suppress all messages
  !         0: normal
  !        >1: increased details 
  !
  subroutine psb_z_null_precdescr(prec,iout,root, verbosity,prefix)
    use psb_penv_mod
    use psb_error_mod
    
    Implicit None

    class(psb_z_null_prec_type), intent(in)   :: prec
    integer(psb_ipk_), intent(in), optional   :: iout
    integer(psb_ipk_), intent(in), optional   :: root
    integer(psb_ipk_), intent(in), optional   :: verbosity
    character(len=*), intent(in), optional  :: prefix
      
    integer(psb_ipk_) :: err_act, nrow, info
    character(len=20)  :: name='z_null_precset'
    character(len=32) :: dprefix, frmtv
    integer(psb_ipk_) :: ni
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iout_, iam, np, root_, verbosity_
    character(1024)    :: prefix_

    call psb_erractionsave(err_act)

    info = psb_success_
   
    if (present(iout)) then 
      iout_ = iout
    else
      iout_ = 6 
    end if
    if (present(root)) then 
      root_ = root
    else
      root_ = psb_root_
    end if
    if (present(prefix)) then
      prefix_ = prefix
    else
      prefix_ = ""
    end if
    if (present(verbosity)) then
      verbosity_ = verbosity
    else
      verbosity_ = 0
    end if
    if (verbosity_ < 0) goto 9998

    ctxt = prec%ctxt
    call psb_info(ctxt,iam,np)
    if (root_ == -1) root_ = iam

    if (iam == root_) &
         &  write(iout_,*) trim(prefix_),' ', trim(prec%desc_prefix()),' ',&
         & 'No preconditioning'

9998 continue
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
    
  end subroutine psb_z_null_precdescr


  subroutine psb_z_null_dump(prec,info,prefix,head)
    use psb_base_mod, only : psb_info
    implicit none 
    class(psb_z_null_prec_type), intent(in) :: prec
    integer(psb_ipk_), intent(out)          :: info
    character(len=*), intent(in), optional  :: prefix,head
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: iout, iam, np, lname
    logical :: isopen
    character(len=80)  :: prefix_
    character(len=120) :: fname ! len should be at least 20 more than
    
    !  len of prefix_ 
    
    info = 0
    ctxt = prec%get_ctxt()
    call psb_info(ctxt,iam,np)
    
    if (present(prefix)) then 
      prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
    else
      prefix_ = "dump_null_z"
    end if
    
    lname = len_trim(prefix_)
    fname = trim(prefix_)
    write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
    lname = lname + 5
    ! Search for an unused unit to write
    iout = 7
    do 
      inquire(unit=iout, opened=isopen)
      if (.not.isopen) exit
      iout = iout + 1
      if (iout > 99) exit
    end do
    if (iout > 99) then 
      write(psb_err_unit,*) 'Error: could not find a free unit for I/O'
      return
    end if
    open(iout,file=fname,iostat=info)
    write(iout,*) 'Null (Identity) Preconditioner. Nothing to be printed, really!'

  end subroutine psb_z_null_dump

  function psb_z_null_sizeof(prec) result(val)

    class(psb_z_null_prec_type), intent(in) :: prec
    integer(psb_epk_) :: val
    
    val = 0

    return
  end function psb_z_null_sizeof


  subroutine psb_z_null_clone(prec,precout,info)
    use psb_const_mod
    use psb_error_mod

    Implicit None

    class(psb_z_null_prec_type), intent(inout) :: prec
    class(psb_z_base_prec_type), allocatable, intent(inout)  :: precout
    integer(psb_ipk_), intent(out)               :: info

    integer(psb_ipk_) :: err_act, i
    character(len=20)  :: name='z_null_clone'

    call psb_erractionsave(err_act)

    info = psb_success_
    if (allocated(precout)) then
      call precout%free(info)
      if (info == psb_success_) deallocate(precout, stat=info)
    end if
    if (info == psb_success_) &
       & allocate(psb_z_null_prec_type :: precout, stat=info)
    if (info /= 0) goto 9999
    select type(pout => precout)
    type is (psb_z_null_prec_type) 
      call pout%set_ctxt(prec%get_ctxt())

    class default
      info = psb_err_internal_error_
    end select
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_z_null_clone


end module psb_z_nullprec
