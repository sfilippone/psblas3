!
!                Parallel Sparse BLAS  version 3.5.1
!      (C) Copyright 2015
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
! File: vecoperation.f90
!
program vecoperation
  use psb_base_mod
  use psb_util_mod
  implicit none

  ! input parameters
  integer(psb_ipk_) :: idim = 100

  ! miscellaneous
  real(psb_dpk_), parameter :: one = 1.d0
  real(psb_dpk_) :: t1, t2

  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! vector
  type(psb_d_vect_type)  :: x,y,z
  ! blacs parameters
  integer(psb_ipk_) :: ictxt, iam, np
  ! other variables
  integer(psb_ipk_) :: nr, nlr, info, i, ii, ib=1
  integer(psb_lpk_) :: nt
  integer(psb_lpk_), allocatable     :: myidx(:)
  real(psb_dpk_)    :: zt(1), dotresult, norm2, norm1, norminf
  character(len=20) :: name,ch_err,readinput
  real(psb_dpk_), allocatable :: vx(:), vy(:), vz(:)
  real(psb_dpk_) :: c

  info=psb_success_


  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then
    call psb_exit(ictxt) ! This should not happen, but just in case
    stop
  endif
  if(psb_get_errstatus() /= 0) goto 9999
  name='vecoperation'
  call psb_set_errverbosity(itwo)
  !
  ! Hello world
  !
  if (iam == psb_root_) then
    write(*,*) 'Welcome to PSBLAS version: ',psb_version_string_
    write(*,*) 'This is the ',trim(name),' sample program'
  end if

  call get_command_argument(1,readinput)
  if (len_trim(readinput) /= 0) read(readinput,*)idim

  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  if (iam == psb_root_) write(psb_out_unit,'("Vector size",I10)')idim

  !
  ! Data distribution
  ! We assume a simple contiguous block distribution
  !
  nt = (idim+np-1)/np
  nr = max(0,min(nt,idim-(iam*nt)))
  nt = nr

  call psb_sum(ictxt,nt)
  if (nt /= idim) then
    write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,idim
    info = -1
    call psb_barrier(ictxt)
    call psb_abort(ictxt)
    return
  end if

  call psb_cdall(ictxt,desc_a,info,nl=nr)
  myidx = desc_a%get_global_indices()
  nlr = size(myidx)

  !
  !  allocate and fill in the vectors
  !
  call psb_barrier(ictxt)
  t1 = psb_wtime()
  ! Allocate memory
  call psb_geall(x,desc_a,info)
  call psb_geall(y,desc_a,info)
  call psb_geall(z,desc_a,info)
  ! Put entries into the vectors
  do ii=1,nlr
    zt(1) = 1.0
    call psb_geins(ib,myidx(ii:),zt(1:),x,desc_a,info)
    zt(1) = 2.0
    call psb_geins(ib,myidx(ii:),zt(1:),y,desc_a,info)
  end do
  ! Assemble
  call psb_cdasb(desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='comm asb rout.'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  if (info == psb_success_) call psb_geasb(x,desc_a,info)
  if (info == psb_success_) call psb_geasb(y,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='vec asb rout.'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_barrier(ictxt)
  t2 = psb_wtime() - t1

  if (iam == psb_root_) write(psb_out_unit,'("Overall vector creation time : ",es12.5)')t2
  if (iam == psb_root_) then
    vx = x%get_vect()
    write(psb_out_unit,'("x = ",es12.1)')vx(:)
    vy = y%get_vect()
    write(psb_out_unit,'("y = ",es12.1)')vy(:)
  end if

  !
  ! Vector operations
  !
  dotresult = psb_gedot(x,y,desc_a,info) ! Dot-product
  if (iam == psb_root_) write(psb_out_unit,'("Dot product result : ",es12.5)')dotresult
  norm1 = psb_norm1(x,desc_a,info)
  norm2 = psb_norm2(x,desc_a,info)
  norminf = psb_normi(x,desc_a,info)
  if (iam == psb_root_) write(psb_out_unit,'("\|x\|_inf : ",es12.5," \|x\|_1 :",es12.5," \|x\|_2",es12.5)')norminf,norm1,norm2
  call psb_geaxpby(1.0_psb_dpk_, x, 1.0_psb_dpk_, y, desc_a, info)  ! \alpha x + \beta y

  if (iam == psb_root_) then
    write(psb_out_unit,'("axpby : x + y")')
    vx = x%get_vect()
    write(psb_out_unit,'("x = ",es12.1)')vx(:)
    vy = y%get_vect()
    write(psb_out_unit,'("y = ",es12.1)')vy(:)
  end if

  call psb_gemlt(x,y,desc_a,info)

  if (iam == psb_root_) then
    write(psb_out_unit,'("mlt : y = x*y ")')
    vx = x%get_vect()
    write(psb_out_unit,'("x = ",es12.1)')vx(:)
    vy = y%get_vect()
    write(psb_out_unit,'("y = ",es12.1)')vy(:)
  end if

  call psb_gediv(x,y,desc_a,info)

  if (iam == psb_root_) then
    write(psb_out_unit,'("div : x = x/y")')
    vx = x%get_vect()
    write(psb_out_unit,'("x = ",es12.1)')vx(:)
    vy = y%get_vect()
    write(psb_out_unit,'("y = ",es12.1)')vy(:)
  end if

  call psb_geinv(x,z,desc_a,info)

  if (iam == psb_root_) then
    write(psb_out_unit,'("inv : z = 1/x")')
    vx = x%get_vect()
    write(psb_out_unit,'("x = ",es12.1)')vx(:)
    vz = z%get_vect()
    write(psb_out_unit,'("z = ",es12.1)')vz(:)
  end if

  c = 1.0/2.0;
  call psb_gecmp(x,c,z,desc_a,info);

  if (iam == psb_root_) then
    write(psb_out_unit,'("|z(i)| >=",es12.1)')c
    vx = x%get_vect()
    write(psb_out_unit,'("x = ",es12.1)')vx(:)
    vz = z%get_vect()
    write(psb_out_unit,'("z = ",es12.1)')vz(:)
  end if

  !
  !  cleanup storage and exit
  !
  call psb_gefree(x,desc_a,info)
  call psb_gefree(y,desc_a,info)
  call psb_cdfree(desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='free routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_exit(ictxt)
  stop

9999 call psb_error(ictxt)

  stop
end program vecoperation
