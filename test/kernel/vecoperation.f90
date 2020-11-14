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
module unittestvector_mod

  use psb_base_mod, only : psb_dpk_, psb_ipk_, psb_desc_type,&
       &  psb_dspmat_type, psb_d_vect_type, dzero, psb_ctxt_type,&
       &  psb_d_base_sparse_mat, psb_d_base_vect_type, psb_i_base_vect_type

  interface psb_gen_const
   module procedure  psb_d_gen_const
  end interface psb_gen_const

contains

  function psb_check_ans(v,val,ictxt) result(ans)
    use psb_base_mod

    implicit none

    type(psb_d_vect_type) :: v
    real(psb_dpk_)        :: val
    type(psb_ctxt_type) :: ictxt
    logical               :: ans

    ! Local variables
    integer(psb_ipk_) :: np, iam, info
    real(psb_dpk_)    :: check
    real(psb_dpk_), allocatable    :: va(:)

    call psb_info(ictxt,iam,np)

    va = v%get_vect()
    va = va - val;

    check = maxval(va);

    call psb_sum(ictxt,check)

    if(check == 0.d0) then
      ans = .true.
    else
      ans = .false.
    end if

  end function psb_check_ans
  !
  !  subroutine to fill a vector with constant entries
  !
  subroutine psb_d_gen_const(v,val,idim,ictxt,desc_a,info)
    use psb_base_mod
    implicit none

    type(psb_d_vect_type) :: v
    type(psb_desc_type)   :: desc_a
    integer(psb_lpk_)     :: idim
    type(psb_ctxt_type) :: ictxt
    integer(psb_ipk_)     :: info
    real(psb_dpk_)        :: val

    ! Local variables
    integer(psb_ipk_), parameter    :: nb=20
    real(psb_dpk_)                  :: zt(nb)
    character(len=20)               :: name, ch_err
    integer(psb_ipk_)               :: np, iam, nr, nt
    integer(psb_ipk_)               :: n,nlr,ib,ii
    integer(psb_ipk_)               :: err_act
    integer(psb_lpk_), allocatable  :: myidx(:)


    info = psb_success_
    name = 'create_constant_vector'
    call psb_erractionsave(err_act)

    call psb_info(ictxt, iam, np)

    n = idim*np         ! The global dimension is the number of process times
                        ! the input size

    ! We use a simple minded block distribution
    nt = (n+np-1)/np
    nr = max(0,min(nt,n-(iam*nt)))
    nt = nr

    call psb_sum(ictxt,nt)
    if (nt /= n) then
      write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,n
      info = -1
      call psb_barrier(ictxt)
      call psb_abort(ictxt)
      return
    end if
    ! Allocate the descriptor with simple minded data distribution
    call psb_cdall(ictxt,desc_a,info,nl=nr)
    ! Allocate the vector on the recently build descriptor
    if (info == psb_success_) call psb_geall(v,desc_a,info)
    ! Check that allocation has gone good
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='allocation rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    myidx = desc_a%get_global_indices()
    nlr = size(myidx)

    do ii=1,nlr,nb
      ib = min(nb,nlr-ii+1)
      zt(:) = val
      call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib),v,desc_a,info)
      if(info /= psb_success_) exit
    end do

    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='insert rout.'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    ! Assembly of communicator and vector
    call psb_cdasb(desc_a,info)
    if (info == psb_success_) call psb_geasb(v,desc_a,info)

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ictxt,err_act)

    return
  end subroutine psb_d_gen_const

end module unittestvector_mod


program vecoperation
  use psb_base_mod
  use psb_util_mod
  use unittestvector_mod
  implicit none

  ! input parameters
  integer(psb_lpk_) :: idim = 100

  ! miscellaneous
  real(psb_dpk_), parameter :: one = 1.d0
  real(psb_dpk_), parameter :: two = 2.d0
  real(psb_dpk_), parameter :: onehalf = 0.5_psb_dpk_
  real(psb_dpk_), parameter :: negativeone = -1.d0
  real(psb_dpk_), parameter :: negativetwo = -2.d0
  real(psb_dpk_), parameter :: negativeonehalf = -0.5_psb_dpk_
  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! vector
  type(psb_d_vect_type)  :: x,y,z
  ! blacs parameters
  type(psb_ctxt_type) :: ictxt
  integer(psb_ipk_) :: iam, np
  ! auxiliary parameters
  integer(psb_ipk_) :: info
  character(len=20) :: name,ch_err,readinput
  real(psb_dpk_)    :: ans
  logical           :: hasitnotfailed
  integer(psb_lpk_), allocatable  :: myidx(:)
  integer(psb_ipk_) :: ib = 1
  real(psb_dpk_)                  :: zt(1)

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
  if (iam == psb_root_) write(psb_out_unit,'("Local vector size",I10)')idim
  if (iam == psb_root_) write(psb_out_unit,'("Global vector size",I10)')np*idim

  !
  ! Test of standard vector operation
  !
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  if (iam == psb_root_) write(psb_out_unit,'("Standard Vector Operations")')
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  ! X = 1
  call psb_d_gen_const(x,one,idim,ictxt,desc_a,info)
  hasitnotfailed = psb_check_ans(x,one,ictxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> Constant vector ")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- Constant vector ")')
  end if
  ! X = 1 , Y = -2, Y = X + Y = 1 -2 = -1
  call psb_d_gen_const(x,one,idim,ictxt,desc_a,info)
  call psb_d_gen_const(y,negativetwo,idim,ictxt,desc_a,info)
  call psb_geaxpby(one,x,one,y,desc_a,info)
  hasitnotfailed = psb_check_ans(y,negativeone,ictxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = X + Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = X + Y ")')
  end if
  ! X = 1 , Y =  2, Y = -X + Y = -1 +2 = 1
  call psb_d_gen_const(x,one,idim,ictxt,desc_a,info)
  call psb_d_gen_const(y,two,idim,ictxt,desc_a,info)
  call psb_geaxpby(negativeone,x,one,y,desc_a,info)
  hasitnotfailed = psb_check_ans(y,one,ictxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = -X + Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = -X + Y ")')
  end if
  ! X = 2 , Y =  -2, Y = 0.5*X + Y = 1 - 2 = -1
  call psb_d_gen_const(x,two,idim,ictxt,desc_a,info)
  call psb_d_gen_const(y,negativetwo,idim,ictxt,desc_a,info)
  call psb_geaxpby(onehalf,x,one,y,desc_a,info)
  hasitnotfailed = psb_check_ans(y,negativeone,ictxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = 0.5 X + Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = 0.5 X + Y ")')
  end if
  ! X = -2 , Y =  1, Z = 0, Z = X + Y = -2 + 1 = -1
  call psb_d_gen_const(x,negativetwo,idim,ictxt,desc_a,info)
  call psb_d_gen_const(y,one,idim,ictxt,desc_a,info)
  call psb_d_gen_const(z,dzero,idim,ictxt,desc_a,info)
  call psb_geaxpby(one,x,one,y,z,desc_a,info)
  hasitnotfailed = psb_check_ans(z,negativeone,ictxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Z = X + Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Z = X + Y ")')
  end if
  ! X = 2 , Y =  1, Z = 0, Z = X - Y = 2 - 1 = 1
  call psb_d_gen_const(x,two,idim,ictxt,desc_a,info)
  call psb_d_gen_const(y,one,idim,ictxt,desc_a,info)
  call psb_d_gen_const(z,dzero,idim,ictxt,desc_a,info)
  call psb_geaxpby(one,x,negativeone,y,z,desc_a,info)
  hasitnotfailed = psb_check_ans(z,one,ictxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Z = X - Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Z = X - Y ")')
  end if
  ! X = 2 , Y =  1, Z = 0, Z = -X + Y = -2 + 1 = -1
  call psb_d_gen_const(x,two,idim,ictxt,desc_a,info)
  call psb_d_gen_const(y,one,idim,ictxt,desc_a,info)
  call psb_d_gen_const(z,dzero,idim,ictxt,desc_a,info)
  call psb_geaxpby(negativeone,x,one,y,z,desc_a,info)
  hasitnotfailed = psb_check_ans(z,negativeone,ictxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Z = -X + Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Z = -X + Y ")')
  end if
  ! X = 2 , Y =  -0.5, Z = 0, Z = X*Y = 2*(-0.5) = -1
  call psb_d_gen_const(x,two,idim,ictxt,desc_a,info)
  call psb_d_gen_const(y,negativeonehalf,idim,ictxt,desc_a,info)
  call psb_d_gen_const(z,dzero,idim,ictxt,desc_a,info)
  call psb_gemlt(one,x,y,dzero,z,desc_a,info)
  hasitnotfailed = psb_check_ans(z,negativeone,ictxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> mlt Z = X*Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- mlt Z = X*Y ")')
  end if
  ! X = 1 , Y =  2, Z = 0, Z = X/Y = 1/2 = 0.5
  call psb_d_gen_const(x,one,idim,ictxt,desc_a,info)
  call psb_d_gen_const(y,two,idim,ictxt,desc_a,info)
  call psb_d_gen_const(z,dzero,idim,ictxt,desc_a,info)
  call psb_gediv(x,y,z,desc_a,info)
  hasitnotfailed = psb_check_ans(z,onehalf,ictxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> div Z = X/Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- div Z = X/Y ")')
  end if
  ! X = -1 , Z =  0, Z = |X| = |-1| = 1
  call psb_d_gen_const(x,negativeone,idim,ictxt,desc_a,info)
  call psb_d_gen_const(z,dzero,idim,ictxt,desc_a,info)
  call psb_geabs(x,z,desc_a,info)
  hasitnotfailed = psb_check_ans(z,one,ictxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> abs Z = |X|")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- abs Z = |X| ")')
  end if
  ! X = 2 , Z =  0, Z = 1/X = 1/2 = 0.5
  call psb_d_gen_const(x,two,idim,ictxt,desc_a,info)
  call psb_d_gen_const(z,dzero,idim,ictxt,desc_a,info)
  call psb_geinv(x,z,desc_a,info)
  hasitnotfailed = psb_check_ans(z,onehalf,ictxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> inv Z = 1/X")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- inv Z = 1/X ")')
  end if
  ! X = 1, Z = 0, c = -2, Z = X + c = -1
  call psb_d_gen_const(x,one,idim,ictxt,desc_a,info)
  call psb_d_gen_const(z,dzero,idim,ictxt,desc_a,info)
  call psb_geaddconst(x,negativetwo,z,desc_a,info)
  hasitnotfailed = psb_check_ans(z,negativeone,ictxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> Add constant Z = X + c")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- Add constant Z = X + c")')
  end if

  !
  ! Vector to field operation
  !
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  if (iam == psb_root_) write(psb_out_unit,'("Vector to Field Operations")')
  if (iam == psb_root_) write(psb_out_unit,'(" ")')

  ! Dot product
  call psb_d_gen_const(x,two,idim,ictxt,desc_a,info)
  call psb_d_gen_const(y,onehalf,idim,ictxt,desc_a,info)
  ans = psb_gedot(x,y,desc_a,info)
  if (iam == psb_root_) then
    if(ans == np*idim) write(psb_out_unit,'("TEST PASSED >>> Dot product")')
    if(ans /= np*idim) write(psb_out_unit,'("TEST FAILED --- Dot product")')
  end if
  ! MaxNorm
  call psb_d_gen_const(x,negativeonehalf,idim,ictxt,desc_a,info)
  ans = psb_geamax(x,desc_a,info)
  if (iam == psb_root_) then
    if(ans == onehalf) write(psb_out_unit,'("TEST PASSED >>> MaxNorm")')
    if(ans /= onehalf) write(psb_out_unit,'("TEST FAILED --- MaxNorm")')
  end if

  call psb_gefree(x,desc_a,info)
  call psb_gefree(y,desc_a,info)
  call psb_gefree(z,desc_a,info)
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
