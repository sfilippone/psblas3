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
!
! File: psb_zkrylovsubspace_mod.F90
!  Module contains Krylov subspace constructors
!
! type psb_zpolkrylov : Polynomial Krylov Space
! v(:)    psb_z_vect_type orthogonal basis of the Krylov space
! h(:,:)  complex(psb_dpk_) projected matrix
! maxsize integer maximum size
! k       current vector
module psb_zkrylovsubspace_mod

  use psb_base_mod
  use psb_util_mod
  use psb_krylov_mod
  implicit none

  type :: psb_zpolkrylov
    type(psb_z_vect_type), allocatable, dimension(:) :: v
    complex(psb_dpk_), allocatable, dimension(:,:) :: h
    integer(psb_ipk_) :: maxsize
    integer(psb_ipk_) :: k
  contains
    procedure, pass(kryl) :: allocate => psb_zallocate_kryl
    procedure, pass(kryl) :: free => psb_zfree_kryl
    procedure, pass(kryl) :: arnoldi => psb_zarnoldi_kryl
    procedure, pass(kryl) :: write => psb_zwrite_kryl
  end type

contains

  ! Allocate the storage for the Krylov space basis and projected matrix
  subroutine psb_zallocate_kryl(kryl,desc_a,maxsize,info)
    use psb_base_mod
    implicit none
    class(psb_zpolkrylov), intent(inout) :: kryl
    type(psb_desc_type), intent(in)    :: desc_a
    integer(psb_ipk_), intent(in)         :: maxsize
    integer(psb_ipk_), intent(out)        :: info

    ! Local
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: err_act
    character(len=20) :: name, ch_err
    integer(psb_lpk_) :: m


    name = 'polkrylov_allocate'
    if (psb_errstatus_fatal()) return
    info = psb_success_
    call psb_erractionsave(err_act)
    ctxt = desc_a%get_context()

    if (.not.allocated(kryl%h)) then
      allocate(kryl%h(maxsize+1,maxsize), stat=info)
    else
      info = -1
    end if
    if ( info /= 0) then
      info=psb_err_from_subroutine_
      ch_err='allocate h'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if (info == psb_success_) call psb_geall(kryl%v,desc_a,info,n=maxsize)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      ch_err='allocate v'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    kryl%k = 0_psb_ipk_
    m = desc_a%get_global_rows()
    if( maxsize > m ) then
      write(psb_out_unit,'("Warning: maxsize > matrix size: defaulting")')
      kryl%maxsize = m-1
    else
      kryl%maxsize = maxsize
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return

  end subroutine psb_zallocate_kryl

  ! Frees the Krylov space basis and projected matrix
  subroutine psb_zfree_kryl(kryl,desc_a,info)
    use psb_base_mod
    implicit none
    class(psb_zpolkrylov), intent(inout) :: kryl
    type(psb_desc_type), intent(in)    :: desc_a
    integer(psb_ipk_), intent(out)        :: info

    ! Local
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: err_act
    character(len=20) :: name, ch_err


    name = 'polkrylov_free'
    if (psb_errstatus_fatal()) return
    info = psb_success_
    call psb_erractionsave(err_act)
    ctxt = desc_a%get_context()

    if (allocated(kryl%h)) deallocate(kryl%h, stat=info)
    if ( info /= 0) then
      info=psb_err_from_subroutine_
      ch_err='free h'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if (info == psb_success_) call psb_gefree(kryl%v,desc_a,info)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      ch_err='free v'
      call psb_errpush(info,name)
      goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)

    return

  end subroutine psb_zfree_kryl

  ! Expand the Krylov space by doing one sweep of Arnoldi
  subroutine psb_zarnoldi_kryl(kryl,a,desc_a,x,l)
    use psb_base_mod
    implicit none
    class(psb_zpolkrylov), intent(inout)  :: kryl
    type(psb_zspmat_type), intent(in)     :: a
    type(psb_z_vect_type), Intent(inout)  :: x
    type(psb_desc_type), intent(in)         :: desc_a
    integer(psb_ipk_), optional, intent(in) :: l

    ! Local
    integer(psb_ipk_) :: l_,i,i1,j,repeat
    integer(psb_lpk_) :: mglob
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: err_act, info
    character(len=20) :: name, ch_err
    real(psb_dpk_) :: nrm0,scal

    info = psb_success_
    name = 'psb_zarnoldi_kryl'
    call psb_erractionsave(err_act)
    ctxt = desc_a%get_context()

    ! If no specificn number of vector is required we just perform a single
    ! Arnoldi step
    if(present(l)) then
      l_ = l
    else
      l_ = 1
    endif

    ! Check if the vector makes sense with the descriptor
    if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      call psb_errpush(info,name)
      goto 9999
    endif
    call psb_chkvect(mglob,lone,x%get_nrows(),lone,lone,desc_a,info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      call psb_errpush(info,name,a_err='psb_chkvect on X')
      goto 9999
    end if

    if (kryl%k == 0) then
      ! We have never done the first step
      nrm0 = psb_genrm2(x,desc_a,info)
      call psb_geaxpby(zone/nrm0,x,zzero,kryl%v(1),desc_a,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_non_
        call psb_errpush(info,name)
        goto 9999
      end if
      kryl%k = 1
    end if

    i = kryl%k
    if ((i+l_) > kryl%maxsize) then
      write(psb_out_unit,'("Warning: Request exceed allocated size, reducing")')
      l_ = kryl%maxsize - i
    end if
    do repeat=i,i+l_
      i1 = i + 1
      ! new candidate vector
      call psb_spmm(zone,a,kryl%v(i),zzero,kryl%v(i1),desc_a,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_non_
        call psb_errpush(info,name)
        goto 9999
      end if
      ! Appy MGS procedure
      do j=1,i
        kryl%h(j,i) = psb_gedot(kryl%v(j),kryl%v(i1),desc_a,info)
        call psb_geaxpby(-kryl%h(j,i),kryl%v(j),zone,kryl%v(i1),desc_a,info)
        if (info /= psb_success_) then
          info=psb_err_from_subroutine_non_
          call psb_errpush(info,name)
          goto 9999
        end if
      end do
      kryl%h(i1,i) = psb_genrm2(kryl%v(i1),desc_a,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_non_
        call psb_errpush(info,name)
        goto 9999
      end if
      scal = zone/kryl%h(i1,i)
      call psb_gescal(kryl%v(i1),scal,kryl%v(i1),desc_a,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_non_
        call psb_errpush(info,name)
        goto 9999
      end if
      i = i1
    end do
    kryl%k = i + 1

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine psb_zarnoldi_kryl

  ! Write Krylov subspace to file
  subroutine psb_zwrite_kryl(kryl,desc_a,filename)
    use psb_base_mod
    use psb_util_mod
    implicit none
    class(psb_zpolkrylov), intent(inout)  :: kryl
    type(psb_desc_type), intent(inout)      :: desc_a
    character(len=*), intent(in)            :: filename

    ! local
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: err_act, info, iam, np, i
    integer(psb_lpk_) :: mglob
    character(len=20) :: name, ch_err
    complex(psb_dpk_), allocatable :: vglobal(:,:)


    info = psb_success_
    name = 'psb_zwrite_kryl'
    call psb_erractionsave(err_act)
    ctxt = desc_a%get_context()

    call psb_info(ctxt,iam,np)

    if (iam == psb_root_) then
      call mm_array_write(kryl%h,"Projected Matrix",info,filename=filename//"_h")
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_non_
        call psb_errpush(info,name)
        goto 9999
      end if
    end if

    ! mglob = desc_a%get_global_rows()
    ! allocate(vglobal(mglob,kryl%k), stat=info)
    ! if (info /= psb_success_) then
    !   info=psb_err_from_subroutine_non_
    !   call psb_errpush(info,name)
    !   goto 9999
    ! end if
    ! do i=1,kryl%k
    !   call psb_gather(vglobal(:,i),kryl%v(i),desc_a,info)
    ! end do
    ! if (iam == psb_root_) then
    !   call mm_array_write(vglobal,"Krylov basis",info,filename=filename//"_v")
    !   if (info /= psb_success_) then
    !     info=psb_err_from_subroutine_non_
    !     call psb_errpush(info,name)
    !     goto 9999
    !   end if
    ! end if
    ! if (allocated(vglobal)) deallocate(vglobal, stat=info)
    ! if (info /= psb_success_) then
    !   info=psb_err_from_subroutine_non_
    !   call psb_errpush(info,name)
    !   goto 9999
    ! end if


    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)
    return

  end subroutine


end module psb_zkrylovsubspace_mod
