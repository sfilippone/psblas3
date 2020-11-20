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
module psb_z_bjacprec

  use psb_z_base_prec_mod
  use psb_z_ilu_fact_mod
  use psb_z_ainv_fact_mod
  use psb_z_invk_fact_mod
  use psb_z_invt_fact_mod

  type, extends(psb_z_base_prec_type)   :: psb_z_bjac_prec_type
    integer(psb_ipk_), allocatable      :: iprcparm(:)
    real(psb_dpk_), allocatable       :: rprcparm(:)
    type(psb_zspmat_type), allocatable  :: av(:)
    type(psb_z_vect_type), allocatable  :: dv, wrk(:)
  contains
    procedure, pass(prec) :: z_apply_v => psb_z_bjac_apply_vect
    procedure, pass(prec) :: z_apply   => psb_z_bjac_apply
    procedure, pass(prec) :: precbld   => psb_z_bjac_precbld
    procedure, pass(prec) :: precinit  => psb_z_bjac_precinit
    procedure, pass(prec) :: precseti  => psb_z_bjac_precseti
    procedure, pass(prec) :: precsetr  => psb_z_bjac_precsetr
    procedure, pass(prec) :: precdescr => psb_z_bjac_precdescr
    procedure, pass(prec) :: dump      => psb_z_bjac_dump
    procedure, pass(prec) :: clone     => psb_z_bjac_clone
    procedure, pass(prec) :: free      => psb_z_bjac_precfree
    procedure, pass(prec) :: sizeof    => psb_z_bjac_sizeof
    procedure, pass(prec) :: get_nzeros => psb_z_bjac_get_nzeros
    procedure, pass(prec) :: allocate_wrk => psb_z_bjac_allocate_wrk
    procedure, pass(prec) :: free_wrk     => psb_z_bjac_free_wrk
    procedure, pass(prec) :: is_allocated_wrk => psb_z_bjac_is_allocated_wrk
  end type psb_z_bjac_prec_type

  character(len=15), parameter, private :: &
       &  fact_names(0:6)=(/'None          ','ILU(0)        ',&
       &  'ILU(n)        ','ILU(eps)      ','AINV(eps)     ',&
       &  'INVT          ','INVK          '/)

  private :: psb_z_bjac_sizeof, psb_z_bjac_precdescr, psb_z_bjac_get_nzeros

  interface
    subroutine psb_z_bjac_dump(prec,info,prefix,head)
      import :: psb_ipk_, psb_desc_type, psb_z_bjac_prec_type, psb_z_vect_type, psb_dpk_
      class(psb_z_bjac_prec_type), intent(in) :: prec
      integer(psb_ipk_), intent(out)                    :: info
      character(len=*), intent(in), optional  :: prefix,head
    end subroutine psb_z_bjac_dump
  end interface

  interface
    subroutine psb_z_bjac_apply_vect(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import :: psb_ipk_, psb_desc_type, psb_z_bjac_prec_type, psb_z_vect_type, psb_dpk_
      type(psb_desc_type),intent(in)    :: desc_data
      class(psb_z_bjac_prec_type), intent(inout)  :: prec
      complex(psb_dpk_),intent(in)         :: alpha,beta
      type(psb_z_vect_type),intent(inout)   :: x
      type(psb_z_vect_type),intent(inout)   :: y
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), optional        :: trans
      complex(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine psb_z_bjac_apply_vect
  end interface

  interface
    subroutine psb_z_bjac_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
      import :: psb_ipk_, psb_desc_type, psb_z_bjac_prec_type, psb_z_vect_type, psb_dpk_

      type(psb_desc_type),intent(in)    :: desc_data
      class(psb_z_bjac_prec_type), intent(inout)  :: prec
      complex(psb_dpk_),intent(in)         :: alpha,beta
      complex(psb_dpk_),intent(inout)      :: x(:)
      complex(psb_dpk_),intent(inout)      :: y(:)
      integer(psb_ipk_), intent(out)              :: info
      character(len=1), optional        :: trans
      complex(psb_dpk_),intent(inout), optional, target :: work(:)
    end subroutine psb_z_bjac_apply
  end interface

  interface
    subroutine psb_z_bjac_precinit(prec,info)
      import :: psb_ipk_, psb_desc_type, psb_z_bjac_prec_type, psb_z_vect_type, psb_dpk_
      class(psb_z_bjac_prec_type),intent(inout) :: prec
      integer(psb_ipk_), intent(out)                     :: info
    end subroutine psb_z_bjac_precinit
  end interface

  interface
    subroutine psb_z_bjac_precbld(a,desc_a,prec,info,amold,vmold,imold)
      import :: psb_ipk_, psb_desc_type, psb_z_bjac_prec_type, psb_z_vect_type, psb_dpk_, &
           & psb_zspmat_type, psb_z_base_sparse_mat, psb_z_base_vect_type, &
           & psb_i_base_vect_type
      type(psb_zspmat_type), intent(in), target :: a
      type(psb_desc_type), intent(inout), target   :: desc_a
      class(psb_z_bjac_prec_type),intent(inout) :: prec
      integer(psb_ipk_), intent(out)                      :: info
      class(psb_z_base_sparse_mat), intent(in), optional :: amold
      class(psb_z_base_vect_type), intent(in), optional  :: vmold
      class(psb_i_base_vect_type), intent(in), optional  :: imold
    end subroutine psb_z_bjac_precbld
  end interface

  interface
    subroutine psb_z_bjac_precseti(prec,what,val,info)
      import :: psb_ipk_, psb_desc_type, psb_z_bjac_prec_type, psb_z_vect_type, psb_dpk_
      class(psb_z_bjac_prec_type),intent(inout) :: prec
      integer(psb_ipk_), intent(in)                      :: what
      integer(psb_ipk_), intent(in)                      :: val
      integer(psb_ipk_), intent(out)                     :: info
    end subroutine psb_z_bjac_precseti
  end interface

  interface
    subroutine psb_z_bjac_precsetr(prec,what,val,info)
      import :: psb_ipk_, psb_desc_type, psb_z_bjac_prec_type, psb_z_vect_type, psb_dpk_
      class(psb_z_bjac_prec_type),intent(inout) :: prec
      integer(psb_ipk_), intent(in)                      :: what
      real(psb_dpk_), intent(in)                      :: val
      integer(psb_ipk_), intent(out)                     :: info
    end subroutine psb_z_bjac_precsetr
  end interface


contains

  subroutine psb_z_bjac_precdescr(prec,iout,root)
    use psb_penv_mod
    use psb_error_mod
    implicit none

    class(psb_z_bjac_prec_type), intent(in) :: prec
    integer(psb_ipk_), intent(in), optional    :: iout
    integer(psb_ipk_), intent(in), optional    :: root

    integer(psb_ipk_) :: err_act, nrow, info
    character(len=20) :: name='z_bjac_precdescr'
    integer(psb_ipk_) :: iout_, ictxt, iam, np, root_

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

    if (.not.allocated(prec%iprcparm)) then
      info = 1124
      call psb_errpush(info,name,a_err="preconditioner")
      goto 9999
    end if

    ictxt = prec%ictxt
    call psb_info(ictxt,iam,np)
    if (root_ == -1) root_ = iam

    if (iam == root_) &
         &  write(iout_,*) trim(prec%desc_prefix()),' ',&
         & 'Block Jacobi with: ',&
         &  fact_names(prec%iprcparm(psb_f_type_))

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_z_bjac_precdescr


  function psb_z_bjac_sizeof(prec) result(val)
    class(psb_z_bjac_prec_type), intent(in) :: prec
    integer(psb_epk_) :: val

    val = 0
    if (allocated(prec%dv)) then
      val = val + (2*psb_sizeof_dp) * prec%dv%get_nrows()
    endif
    if (allocated(prec%av)) then
      val = val + prec%av(psb_l_pr_)%sizeof()
      val = val + prec%av(psb_u_pr_)%sizeof()
    endif
    return
  end function psb_z_bjac_sizeof

  function psb_z_bjac_get_nzeros(prec) result(val)

    class(psb_z_bjac_prec_type), intent(in) :: prec
    integer(psb_epk_) :: val

    val = 0
    if (allocated(prec%dv)) then
      val = val + prec%dv%get_nrows()
    endif
    if (allocated(prec%av)) then
      val = val + prec%av(psb_l_pr_)%get_nzeros()
      val = val + prec%av(psb_u_pr_)%get_nzeros()
    endif
    return
  end function psb_z_bjac_get_nzeros


  subroutine psb_z_bjac_precfree(prec,info)

    Implicit None

    class(psb_z_bjac_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)                :: info

    integer(psb_ipk_) :: err_act, i
    character(len=20)  :: name='z_bjac_precfree'

    call psb_erractionsave(err_act)

    info = psb_success_
    if (allocated(prec%av)) then
      do i=1,size(prec%av)
        call prec%av(i)%free()
      enddo
      deallocate(prec%av,stat=info)
    end if

    if (allocated(prec%dv)) then
      call prec%dv%free(info)
      if (info == 0) deallocate(prec%dv,stat=info)
    end if
    if (allocated(prec%iprcparm)) then
      deallocate(prec%iprcparm,stat=info)
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_z_bjac_precfree


  subroutine psb_z_bjac_clone(prec,precout,info)
    use psb_error_mod
    use psb_realloc_mod
    Implicit None

    class(psb_z_bjac_prec_type), intent(inout)              :: prec
    class(psb_z_base_prec_type), allocatable, intent(inout) :: precout
    integer(psb_ipk_), intent(out)               :: info

    integer(psb_ipk_) :: err_act, i
    character(len=20)  :: name='z_bjac_clone'

    call psb_erractionsave(err_act)

    info = psb_success_
    if (allocated(precout)) then
      call precout%free(info)
      if (info == psb_success_) deallocate(precout, stat=info)
    end if
    if (info == psb_success_) &
         & allocate(psb_z_bjac_prec_type :: precout, stat=info)
    if (info /= 0) goto 9999
    select type(pout => precout)
    type is (psb_z_bjac_prec_type)
      call pout%set_ctxt(prec%get_ctxt())

      if (allocated(prec%av)) then
        allocate(pout%av(size(prec%av)),stat=info)
        do i=1,size(prec%av)
          if (info /= psb_success_) exit
          call prec%av(i)%clone(pout%av(i),info)
        enddo
        if (info /= psb_success_) goto 9999
      end if

      if (allocated(prec%dv)) then
        allocate(pout%dv,stat=info)
        if (info == 0) call prec%dv%clone(pout%dv,info)
      end if
      class default
      info = psb_err_internal_error_
    end select
    if (info /= 0) goto 9999

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_z_bjac_clone

  subroutine psb_z_bjac_allocate_wrk(prec,info,vmold,desc)
    use psb_base_mod
    implicit none

    ! Arguments
    class(psb_z_bjac_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)        :: info
    class(psb_z_base_vect_type), intent(in), optional  :: vmold
    type(psb_desc_type), intent(in), optional :: desc

    ! Local variables
    integer(psb_ipk_) :: err_act, i
    character(len=20)   :: name

    info=psb_success_
    name = 'psb_z_allocate_wrk'
    call psb_erractionsave(err_act)

    if (psb_get_errstatus().ne.0) goto 9999
    if (allocated(prec%wrk)) then
      if (size(prec%wrk)<2) then
        do i=1,size(prec%wrk)
          if (info == 0) call prec%wrk(i)%free(info)
        end do
        if (info == 0) deallocate(prec%wrk,stat=info)
      end if
    end if
    if (info /= 0) then
      info = psb_err_internal_error_; call psb_errpush(info,name,a_err="deallocate");   goto 9999
    end if

    if (.not.allocated(prec%wrk)) then
      if (.not.present(desc)) then
        info = psb_err_internal_error_; call psb_errpush(info,name,a_err="no desc?");   goto 9999
      end if
      allocate(prec%wrk(2),stat=info)
      do i=1, 2
        if (info == 0)   call psb_geall(prec%wrk(i),desc,info)
        if (info == 0)   call psb_geasb(prec%wrk(i),desc,info,mold=vmold,scratch=.true.)
      end do
    end if
    if (info /= 0) then
      info = psb_err_internal_error_; call psb_errpush(info,name,a_err="allocate");   goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_z_bjac_allocate_wrk

  subroutine psb_z_bjac_free_wrk(prec,info)
    use psb_base_mod
    implicit none

    ! Arguments
    class(psb_z_bjac_prec_type), intent(inout) :: prec
    integer(psb_ipk_), intent(out)        :: info

    ! Local variables
    integer(psb_ipk_) :: err_act
    integer(psb_ipk_) :: i
    character(len=20) :: name

    info=psb_success_
    name = 'psb_z_allocate_wrk'
    call psb_erractionsave(err_act)

    if (psb_get_errstatus().ne.0) goto 9999

    info = psb_success_
    if (allocated(prec%wrk)) then
      do i=1,size(prec%wrk)
        if (info == 0) call prec%wrk(i)%free(info)
      end do
      if (info == 0) deallocate(prec%wrk,stat=info)
    end if
    if (info /= 0) then
      info = psb_err_internal_error_; call psb_errpush(info,name,a_err="deallocate");   goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_z_bjac_free_wrk

  function psb_z_bjac_is_allocated_wrk(prec) result(res)
    use psb_base_mod
    implicit none

    ! Arguments
    class(psb_z_bjac_prec_type), intent(in) :: prec
    logical :: res

    ! In the base version we can say yes, because
    ! there is nothing to allocate

    res = allocated(prec%wrk)

  end function psb_z_bjac_is_allocated_wrk


end module psb_z_bjacprec
