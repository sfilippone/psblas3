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

  use psb_base_mod

  interface psb_gen_const
   module procedure  psb_d_gen_const, psb_s_gen_const, &
       &  psb_c_gen_const, psb_z_gen_const, &
       &  psb_d_gen_const_multi, psb_s_gen_const_multi
  end interface psb_gen_const

  interface psb_check_ans
   module procedure  psb_d_check_ans_v, psb_c_check_ans_v, &
       &  psb_z_check_ans_v, psb_s_check_ans_v, &
       &  psb_d_check_ans_mv, psb_s_check_ans_mv, &
       &  psb_c_check_ans_mv, psb_z_check_ans_mv
  end interface psb_check_ans

contains

  function psb_d_check_ans_v(v,val,ctxt) result(ans)
    use psb_base_mod

    implicit none

    type(psb_d_vect_type) :: v
    real(psb_dpk_)        :: val
    type(psb_ctxt_type) :: ctxt
    logical               :: ans

    ! Local variables
    integer(psb_ipk_) :: np, iam, info
    real(psb_dpk_)    :: check
    real(psb_dpk_), allocatable    :: va(:)

    call psb_info(ctxt,iam,np)

    va = v%get_vect()
    va = va - val;

    check = maxval(va);

    call psb_sum(ctxt,check)

    if(check == 0.d0) then
      ans = .true.
    else
      ans = .false.
    end if

  end function psb_d_check_ans_v

  function psb_s_check_ans_v(v,val,ctxt) result(ans)
    use psb_base_mod

    implicit none

    type(psb_s_vect_type) :: v
    real(psb_spk_)        :: val
    type(psb_ctxt_type) :: ctxt
    logical               :: ans

    ! Local variables
    integer(psb_ipk_) :: np, iam, info
    real(psb_spk_)    :: check
    real(psb_spk_), allocatable    :: va(:)

    call psb_info(ctxt,iam,np)

    va = v%get_vect()
    va = va - val;

    check = maxval(va);

    call psb_sum(ctxt,check)

    if(check == 0.e0) then
      ans = .true.
    else
      ans = .false.
    end if

  end function psb_s_check_ans_v

  function psb_c_check_ans_v(v,val,ctxt) result(ans)
    use psb_base_mod

    implicit none

    type(psb_c_vect_type) :: v
    complex(psb_spk_)     :: val
    type(psb_ctxt_type) :: ctxt
    logical               :: ans

    ! Local variables
    integer(psb_ipk_) :: np, iam, info
    real(psb_spk_)    :: check
    complex(psb_spk_), allocatable    :: va(:)

    call psb_info(ctxt,iam,np)

    va = v%get_vect()
    va = va - val;

    check = maxval(abs(va));

    call psb_sum(ctxt,check)

    if(check == 0.e0) then
      ans = .true.
    else
      ans = .false.
    end if

  end function psb_c_check_ans_v

  function psb_z_check_ans_v(v,val,ctxt) result(ans)
    use psb_base_mod

    implicit none

    type(psb_z_vect_type) :: v
    complex(psb_dpk_)     :: val
    type(psb_ctxt_type) :: ctxt
    logical               :: ans

    ! Local variables
    integer(psb_ipk_) :: np, iam, info
    real(psb_dpk_)    :: check
    complex(psb_dpk_), allocatable    :: va(:)

    call psb_info(ctxt,iam,np)

    va = v%get_vect()
    va = va - val;

    check = maxval(abs(va));

    call psb_sum(ctxt,check)

    if(check == 0.d0) then
      ans = .true.
    else
      ans = .false.
    end if

  end function psb_z_check_ans_v

  function psb_d_check_ans_mv(v,val,ctxt) result(ans)
    use psb_base_mod

    implicit none

    type(psb_d_multivect_type) :: v
    real(psb_dpk_)        :: val
    type(psb_ctxt_type) :: ctxt
    logical               :: ans

    ! Local variables
    integer(psb_ipk_) :: np, iam, info
    real(psb_dpk_)    :: check
    real(psb_dpk_), allocatable    :: va(:,:)

    call psb_info(ctxt,iam,np)

    va = v%get_vect()
    va = va - val;

    check = maxval(va);

    call psb_sum(ctxt,check)

    if(check == 0.d0) then
      ans = .true.
    else
      ans = .false.
    end if

  end function psb_d_check_ans_mv

  function psb_s_check_ans_mv(v,val,ctxt) result(ans)
    use psb_base_mod

    implicit none

    type(psb_s_multivect_type) :: v
    real(psb_spk_)        :: val
    type(psb_ctxt_type) :: ctxt
    logical               :: ans

    ! Local variables
    integer(psb_ipk_) :: np, iam, info
    real(psb_spk_)    :: check
    real(psb_spk_), allocatable    :: va(:,:)

    call psb_info(ctxt,iam,np)

    va = v%get_vect()
    va = va - val;

    check = maxval(va);

    call psb_sum(ctxt,check)

    if(check == 0.e0) then
      ans = .true.
    else
      ans = .false.
    end if

  end function psb_s_check_ans_mv

  function psb_c_check_ans_mv(v,val,ctxt) result(ans)
    use psb_base_mod

    implicit none

    type(psb_c_multivect_type) :: v
    complex(psb_spk_)     :: val
    type(psb_ctxt_type) :: ctxt
    logical               :: ans

    ! Local variables
    integer(psb_ipk_) :: np, iam, info
    real(psb_spk_)    :: check
    complex(psb_spk_), allocatable    :: va(:,:)

    call psb_info(ctxt,iam,np)

    va = v%get_vect()
    va = va - val;

    check = maxval(abs(va));

    call psb_sum(ctxt,check)

    if(check == 0.e0) then
      ans = .true.
    else
      ans = .false.
    end if

  end function psb_c_check_ans_mv

  function psb_z_check_ans_mv(v,val,ctxt) result(ans)
    use psb_base_mod

    implicit none

    type(psb_z_multivect_type) :: v
    complex(psb_dpk_)     :: val
    type(psb_ctxt_type) :: ctxt
    logical               :: ans

    ! Local variables
    integer(psb_ipk_) :: np, iam, info
    real(psb_dpk_)    :: check
    complex(psb_dpk_), allocatable    :: va(:,:)

    call psb_info(ctxt,iam,np)

    va = v%get_vect()
    va = va - val;

    check = maxval(abs(va));

    call psb_sum(ctxt,check)

    if(check == 0.d0) then
      ans = .true.
    else
      ans = .false.
    end if

  end function psb_z_check_ans_mv
  !
  !  subroutine to fill a vector with constant entries
  !
  subroutine psb_z_gen_const(v,val,idim,ctxt,desc_a,info)
    use psb_base_mod
    implicit none

    type(psb_z_vect_type) :: v
    type(psb_desc_type)   :: desc_a
    integer(psb_lpk_)     :: idim
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)     :: info
    complex(psb_dpk_)        :: val

    ! Local variables
    integer(psb_ipk_), parameter    :: nb=20
    complex(psb_dpk_)                  :: zt(nb)
    character(len=20)               :: name, ch_err
    integer(psb_ipk_)               :: np, iam, nr, nt
    integer(psb_ipk_)               :: n,nlr,ib,ii
    integer(psb_ipk_)               :: err_act
    integer(psb_lpk_), allocatable  :: myidx(:)


    info = psb_success_
    name = 'create_constant_vector'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)

    n = idim*np         ! The global dimension is the number of process times
                        ! the input size

    ! We use a simple minded block distribution
    nt = (n+np-1)/np
    nr = max(0,min(nt,n-(iam*nt)))
    nt = nr

    call psb_sum(ctxt,nt)
    if (nt /= n) then
      write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,n
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end if
    ! Allocate the descriptor with simple minded data distribution
    call psb_cdall(ctxt,desc_a,info,nl=nr)
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

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psb_z_gen_const
  !
  !  subroutine to fill a vector with constant entries
  !
  subroutine psb_c_gen_const(v,val,idim,ctxt,desc_a,info)
    use psb_base_mod
    implicit none

    type(psb_c_vect_type) :: v
    type(psb_desc_type)   :: desc_a
    integer(psb_lpk_)     :: idim
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)     :: info
    complex(psb_spk_)        :: val

    ! Local variables
    integer(psb_ipk_), parameter    :: nb=20
    complex(psb_spk_)                  :: zt(nb)
    character(len=20)               :: name, ch_err
    integer(psb_ipk_)               :: np, iam, nr, nt
    integer(psb_ipk_)               :: n,nlr,ib,ii
    integer(psb_ipk_)               :: err_act
    integer(psb_lpk_), allocatable  :: myidx(:)


    info = psb_success_
    name = 'create_constant_vector'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)

    n = idim*np         ! The global dimension is the number of process times
                        ! the input size

    ! We use a simple minded block distribution
    nt = (n+np-1)/np
    nr = max(0,min(nt,n-(iam*nt)))
    nt = nr

    call psb_sum(ctxt,nt)
    if (nt /= n) then
      write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,n
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end if
    ! Allocate the descriptor with simple minded data distribution
    call psb_cdall(ctxt,desc_a,info,nl=nr)
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

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psb_c_gen_const
  !
  !  subroutine to fill a vector with constant entries
  !
  subroutine psb_s_gen_const(v,val,idim,ctxt,desc_a,info)
    use psb_base_mod
    implicit none

    type(psb_s_vect_type) :: v
    type(psb_desc_type)   :: desc_a
    integer(psb_lpk_)     :: idim
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)     :: info
    real(psb_spk_)        :: val

    ! Local variables
    integer(psb_ipk_), parameter    :: nb=20
    real(psb_spk_)                  :: zt(nb)
    character(len=20)               :: name, ch_err
    integer(psb_ipk_)               :: np, iam, nr, nt
    integer(psb_ipk_)               :: n,nlr,ib,ii
    integer(psb_ipk_)               :: err_act
    integer(psb_lpk_), allocatable  :: myidx(:)


    info = psb_success_
    name = 'create_constant_vector'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)

    n = idim*np         ! The global dimension is the number of process times
                        ! the input size

    ! We use a simple minded block distribution
    nt = (n+np-1)/np
    nr = max(0,min(nt,n-(iam*nt)))
    nt = nr

    call psb_sum(ctxt,nt)
    if (nt /= n) then
      write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,n
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end if
    ! Allocate the descriptor with simple minded data distribution
    call psb_cdall(ctxt,desc_a,info,nl=nr)
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

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psb_s_gen_const
  !
  !  subroutine to fill a vector with constant entries
  !
  subroutine psb_d_gen_const(v,val,idim,ctxt,desc_a,info)
    use psb_base_mod
    implicit none

    type(psb_d_vect_type) :: v
    type(psb_desc_type)   :: desc_a
    integer(psb_lpk_)     :: idim
    type(psb_ctxt_type) :: ctxt
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

    call psb_info(ctxt, iam, np)

    n = idim*np         ! The global dimension is the number of process times
                        ! the input size

    ! We use a simple minded block distribution
    nt = (n+np-1)/np
    nr = max(0,min(nt,n-(iam*nt)))
    nt = nr

    call psb_sum(ctxt,nt)
    if (nt /= n) then
      write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,n
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end if
    ! Allocate the descriptor with simple minded data distribution
    call psb_cdall(ctxt,desc_a,info,nl=nr)
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

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psb_d_gen_const
  !
  !  subroutine to fill a multivectorvector with constant entries
  !
  subroutine psb_d_gen_const_multi(v,val,idim,jdim,ctxt,desc_a,info)
    use psb_base_mod
    implicit none

    type(psb_d_multivect_type) :: v
    type(psb_desc_type)   :: desc_a
    integer(psb_lpk_)     :: idim
    integer(psb_ipk_)     :: jdim !! Number of columns of the multivector
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)     :: info
    real(psb_dpk_)        :: val

    ! Local variables
    integer(psb_ipk_), parameter    :: nb=20
    real(psb_dpk_)                  :: zt(nb,jdim) ! Temporary array to fill the vector
    character(len=40)               :: name, ch_err
    integer(psb_ipk_)               :: np, iam, nr, nt
    integer(psb_ipk_)               :: n,nlr,ib,ii
    integer(psb_ipk_)               :: err_act
    integer(psb_lpk_), allocatable  :: myidx(:)


    info = psb_success_
    name = 'create_constant_multivector'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)

    n = idim*np         ! The global dimension is the number of process times
                        ! the input size

    ! We use a simple minded block distribution
    nt = (n+np-1)/np
    nr = max(0,min(nt,n-(iam*nt)))
    nt = nr

    call psb_sum(ctxt,nt)
    if (nt /= n) then
      write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,n
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end if
    ! Allocate the descriptor with simple minded data distribution
    call psb_cdall(ctxt,desc_a,info,nl=nr)
    ! Allocate the vector on the recently build descriptor
    if (info == psb_success_) call psb_geall(v,desc_a,info,n=jdim)
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
      zt(:,:) = val
      call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib,1:jdim),v,desc_a,info)
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

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psb_d_gen_const_multi
  !
  !  subroutine to fill a multivectorvector with constant entries
  !
  subroutine psb_s_gen_const_multi(v,val,idim,jdim,ctxt,desc_a,info)
    use psb_base_mod
    implicit none

    type(psb_s_multivect_type) :: v
    type(psb_desc_type)   :: desc_a
    integer(psb_lpk_)     :: idim
    integer(psb_ipk_)     :: jdim !! Number of columns of the multivector
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)     :: info
    real(psb_spk_)        :: val

    ! Local variables
    integer(psb_ipk_), parameter    :: nb=20
    real(psb_spk_)                  :: zt(nb,jdim) ! Temporary array to fill the vector
    character(len=40)               :: name, ch_err
    integer(psb_ipk_)               :: np, iam, nr, nt
    integer(psb_ipk_)               :: n,nlr,ib,ii
    integer(psb_ipk_)               :: err_act
    integer(psb_lpk_), allocatable  :: myidx(:)


    info = psb_success_
    name = 'create_constant_multivector'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)

    n = idim*np         ! The global dimension is the number of process times
                        ! the input size

    ! We use a simple minded block distribution
    nt = (n+np-1)/np
    nr = max(0,min(nt,n-(iam*nt)))
    nt = nr

    call psb_sum(ctxt,nt)
    if (nt /= n) then
      write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,n
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end if
    ! Allocate the descriptor with simple minded data distribution
    call psb_cdall(ctxt,desc_a,info,nl=nr)
    ! Allocate the vector on the recently build descriptor
    if (info == psb_success_) call psb_geall(v,desc_a,info,n=jdim)
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
      zt(:,:) = val
      call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib,1:jdim),v,desc_a,info)
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

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psb_s_gen_const_multi
    !
  !  subroutine to fill a multivectorvector with constant entries
  !
  subroutine psb_c_gen_const_multi(v,val,idim,jdim,ctxt,desc_a,info)
    use psb_base_mod
    implicit none

    type(psb_c_multivect_type) :: v
    type(psb_desc_type)   :: desc_a
    integer(psb_lpk_)     :: idim
    integer(psb_ipk_)     :: jdim !! Number of columns of the multivector
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)     :: info
    complex(psb_spk_)        :: val

    ! Local variables
    integer(psb_ipk_), parameter    :: nb=20
    complex(psb_spk_)                  :: zt(nb,jdim) ! Temporary array to fill the vector
    character(len=40)               :: name, ch_err
    integer(psb_ipk_)               :: np, iam, nr, nt
    integer(psb_ipk_)               :: n,nlr,ib,ii
    integer(psb_ipk_)               :: err_act
    integer(psb_lpk_), allocatable  :: myidx(:)


    info = psb_success_
    name = 'create_constant_multivector'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)

    n = idim*np         ! The global dimension is the number of process times
                        ! the input size

    ! We use a simple minded block distribution
    nt = (n+np-1)/np
    nr = max(0,min(nt,n-(iam*nt)))
    nt = nr

    call psb_sum(ctxt,nt)
    if (nt /= n) then
      write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,n
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end if
    ! Allocate the descriptor with simple minded data distribution
    call psb_cdall(ctxt,desc_a,info,nl=nr)
    ! Allocate the vector on the recently build descriptor
    if (info == psb_success_) call psb_geall(v,desc_a,info,n=jdim)
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
      zt(:,:) = val
      call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib,1:jdim),v,desc_a,info)
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

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psb_c_gen_const_multi
    !
  !  subroutine to fill a multivectorvector with constant entries
  !
  subroutine psb_z_gen_const_multi(v,val,idim,jdim,ctxt,desc_a,info)
    use psb_base_mod
    implicit none

    type(psb_z_multivect_type) :: v
    type(psb_desc_type)   :: desc_a
    integer(psb_lpk_)     :: idim
    integer(psb_ipk_)     :: jdim !! Number of columns of the multivector
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_)     :: info
    complex(psb_dpk_)        :: val

    ! Local variables
    integer(psb_ipk_), parameter    :: nb=20
    complex(psb_dpk_)               :: zt(nb,jdim) ! Temporary array to fill the vector
    character(len=40)               :: name, ch_err
    integer(psb_ipk_)               :: np, iam, nr, nt
    integer(psb_ipk_)               :: n,nlr,ib,ii
    integer(psb_ipk_)               :: err_act
    integer(psb_lpk_), allocatable  :: myidx(:)


    info = psb_success_
    name = 'create_constant_multivector'
    call psb_erractionsave(err_act)

    call psb_info(ctxt, iam, np)

    n = idim*np         ! The global dimension is the number of process times
                        ! the input size

    ! We use a simple minded block distribution
    nt = (n+np-1)/np
    nr = max(0,min(nt,n-(iam*nt)))
    nt = nr

    call psb_sum(ctxt,nt)
    if (nt /= n) then
      write(psb_err_unit,*) iam, 'Initialization error ',nr,nt,n
      info = -1
      call psb_barrier(ctxt)
      call psb_abort(ctxt)
      return
    end if
    ! Allocate the descriptor with simple minded data distribution
    call psb_cdall(ctxt,desc_a,info,nl=nr)
    ! Allocate the vector on the recently build descriptor
    if (info == psb_success_) call psb_geall(v,desc_a,info,n=jdim)
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
      zt(:,:) = val
      call psb_geins(ib,myidx(ii:ii+ib-1),zt(1:ib,1:jdim),v,desc_a,info)
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

9999 call psb_error_handler(ctxt,err_act)

    return
  end subroutine psb_z_gen_const_multi

end module unittestvector_mod


program vecoperation
  use psb_base_mod
  use psb_util_mod
  use unittestvector_mod
  implicit none

  ! input parameters
  integer(psb_lpk_) :: idim = 100 ! Local vector size
  integer(psb_ipk_) :: nmv = 10 ! Number of columns of the multivector

  ! miscellaneous
  real(psb_dpk_), parameter :: one = 1.d0
  real(psb_dpk_), parameter :: two = 2.d0
  real(psb_dpk_), parameter :: onehalf = 0.5_psb_dpk_
  real(psb_dpk_), parameter :: negativeone = -1.d0
  real(psb_dpk_), parameter :: negativetwo = -2.d0
  real(psb_dpk_), parameter :: negativeonehalf = -0.5_psb_dpk_

  real(psb_spk_), parameter :: stwo = 2.e0
  real(psb_spk_), parameter :: sonehalf = 0.5_psb_spk_
  real(psb_spk_), parameter :: snegativeone = -1.e0
  real(psb_spk_), parameter :: snegativetwo = -2.e0
  real(psb_spk_), parameter :: snegativeonehalf = -0.5_psb_spk_

  complex(psb_spk_), parameter :: ctwo = (2.e0, 0.e0)
  complex(psb_spk_), parameter :: cnegativeone = (-1.e0, 0.e0)
  complex(psb_spk_), parameter :: cnegativetwo = (-2.e0, 0.e0)
  complex(psb_spk_), parameter :: conehalf = (0.5e0, 0.e0)

  complex(psb_dpk_), parameter :: ztwo = (2.d0, 0.d0)
  complex(psb_dpk_), parameter :: znegativeone = (-1.d0, 0.d0)
  complex(psb_dpk_), parameter :: znegativetwo = (-2.d0, 0.d0)
  complex(psb_dpk_), parameter :: zonehalf = (0.5d0, 0.d0)


  ! descriptor
  type(psb_desc_type)   :: desc_a
  ! vector
  type(psb_d_vect_type)  :: x, y, z
  type(psb_s_vect_type)  :: sx, sy, sz
  type(psb_c_vect_type)  :: cx, cy, cz
  type(psb_z_vect_type)  :: zx, zy, zz
  ! multivector
  type(psb_d_multivect_type) :: mv1, mv2
  type(psb_s_multivect_type) :: smv1, smv2
  type(psb_c_multivect_type) :: cmv1, cmv2
  type(psb_z_multivect_type) :: zmv1, zmv2
  ! blacs parameters
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: iam, np
  ! auxiliary parameters
  integer(psb_ipk_) :: ii
  integer(psb_ipk_) :: info
  character(len=20) :: name,ch_err,readinput
  real(psb_dpk_)    :: ans
  real(psb_dpk_), allocatable :: ansmv(:)
  logical           :: hasitnotfailed
  integer(psb_lpk_), allocatable  :: myidx(:)
  integer(psb_ipk_) :: ib = 1
  real(psb_dpk_)                  :: zt(1)

  info=psb_success_
  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)

  if (iam < 0) then
    call psb_exit(ctxt) ! This should not happen, but just in case
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
  if (iam == psb_root_) write(psb_out_unit,'("---- Standard Vector Operations ----")')
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  if (iam == psb_root_) write(psb_out_unit,'(" Real Double Precision")')
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  ! X = 1
  call psb_d_gen_const(x,one,idim,ctxt,desc_a,info)
  hasitnotfailed = psb_check_ans(x,one,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> Constant vector ")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- Constant vector ")')
  end if
  ! X = 1 , Y = -2, Y = X + Y = 1 -2 = -1
  call psb_d_gen_const(x,one,idim,ctxt,desc_a,info)
  call psb_d_gen_const(y,negativetwo,idim,ctxt,desc_a,info)
  call psb_geaxpby(one,x,one,y,desc_a,info)
  hasitnotfailed = psb_check_ans(y,negativeone,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = X + Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = X + Y ")')
  end if
  ! X = 1 , Y =  2, Y = -X + Y = -1 +2 = 1
  call psb_d_gen_const(x,one,idim,ctxt,desc_a,info)
  call psb_d_gen_const(y,two,idim,ctxt,desc_a,info)
  call psb_geaxpby(negativeone,x,one,y,desc_a,info)
  hasitnotfailed = psb_check_ans(y,one,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = -X + Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = -X + Y ")')
  end if
  ! X = 2 , Y =  -2, Y = 0.5*X + Y = 1 - 2 = -1
  call psb_d_gen_const(x,two,idim,ctxt,desc_a,info)
  call psb_d_gen_const(y,negativetwo,idim,ctxt,desc_a,info)
  call psb_geaxpby(onehalf,x,one,y,desc_a,info)
  hasitnotfailed = psb_check_ans(y,negativeone,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = 0.5 X + Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = 0.5 X + Y ")')
  end if
  ! X = -2 , Y =  1, Z = 0, Z = X + Y = -2 + 1 = -1
  call psb_d_gen_const(x,negativetwo,idim,ctxt,desc_a,info)
  call psb_d_gen_const(y,one,idim,ctxt,desc_a,info)
  call psb_d_gen_const(z,dzero,idim,ctxt,desc_a,info)
  call psb_geaxpby(one,x,one,y,z,desc_a,info)
  hasitnotfailed = psb_check_ans(z,negativeone,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Z = X + Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Z = X + Y ")')
  end if
  ! X = 2 , Y =  1, Z = 0, Z = X - Y = 2 - 1 = 1
  call psb_d_gen_const(x,two,idim,ctxt,desc_a,info)
  call psb_d_gen_const(y,one,idim,ctxt,desc_a,info)
  call psb_d_gen_const(z,dzero,idim,ctxt,desc_a,info)
  call psb_geaxpby(one,x,negativeone,y,z,desc_a,info)
  hasitnotfailed = psb_check_ans(z,one,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Z = X - Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Z = X - Y ")')
  end if
  ! X = 2 , Y =  -0.5, Z = 0, Z = X*Y = 2*(-0.5) = -1
  call psb_d_gen_const(x,two,idim,ctxt,desc_a,info)
  call psb_d_gen_const(y,negativeonehalf,idim,ctxt,desc_a,info)
  call psb_d_gen_const(z,dzero,idim,ctxt,desc_a,info)
  call psb_gemlt(one,x,y,dzero,z,desc_a,info)
  hasitnotfailed = psb_check_ans(z,negativeone,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> mlt Z = X*Y")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- mlt Z = X*Y ")')
  end if

  ! Single Precision Real
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  if (iam == psb_root_) write(psb_out_unit,'(" Real Single Precision")')
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  ! X = 1
  call psb_s_gen_const(sx,sone,idim,ctxt,desc_a,info)
  hasitnotfailed = psb_check_ans(sx,sone,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> Constant vector (Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- Constant vector (Single Precision)")')
  end if
  ! X = 2 , Y =  -2, Y = 0.5*X + Y = 1 - 2 = -1
  call psb_s_gen_const(sx,stwo,idim,ctxt,desc_a,info)
  call psb_s_gen_const(sy,snegativetwo,idim,ctxt,desc_a,info)
  call psb_geaxpby(sonehalf,sx,sone,sy,desc_a,info)
  hasitnotfailed = psb_check_ans(sy,snegativeone,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = 0.5 X + Y (Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = 0.5 X + Y (Single Precision)")')
  end if
  ! X = 2 , Y =  -2, Y = X + Y = 0
  call psb_s_gen_const(sx,stwo,idim,ctxt,desc_a,info)
  call psb_s_gen_const(sy,snegativetwo,idim,ctxt,desc_a,info)
  call psb_geaxpby(sone,sx,sone,sy,desc_a,info)
  hasitnotfailed = psb_check_ans(sy,(0.0_psb_spk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = X + Y (Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = X + Y (Single Precision)")')
  end if
  ! X = 2 , Y =  -2, Y = 0.5*X + Y = 1 - 2 = -1
  call psb_s_gen_const(sx,stwo,idim,ctxt,desc_a,info)
  call psb_s_gen_const(sy,snegativetwo,idim,ctxt,desc_a,info)
  call psb_geaxpby(sonehalf,sx,sone,sy,desc_a,info)
  hasitnotfailed = psb_check_ans(sy,snegativeone,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = 0.5 X + Y (Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = 0.5 X + Y (Single Precision)")')
  end if
  ! X = 2 , Y =  -2, Y = 0.5*X + Y = 1 - 2 = -1
  call psb_s_gen_const(sx,stwo,idim,ctxt,desc_a,info)
  call psb_s_gen_const(sy,snegativetwo,idim,ctxt,desc_a,info)
  call psb_geaxpby(sonehalf,sx,sone,sy,desc_a,info)
  hasitnotfailed = psb_check_ans(sy,snegativeone,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = 0.5 X + Y (Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = 0.5 X + Y (Single Precision)")')
  end if
  ! X = 2 , Y =  1, Z = 0, Z = X - Y = 2 - 1 = 1
  call psb_s_gen_const(sx,stwo,idim,ctxt,desc_a,info)
  call psb_s_gen_const(sy,sone,idim,ctxt,desc_a,info)
  call psb_s_gen_const(sz,szero,idim,ctxt,desc_a,info)
  call psb_geaxpby(sone,sx,snegativeone,sy,sz,desc_a,info)
  hasitnotfailed = psb_check_ans(sz,sone,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Z = X - Y (Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Z = X - Y (Single Precision)")')
  end if
  ! X = 2 , Y =  -0.5, Z = 0, Z = X*Y = -1
  call psb_s_gen_const(sx,stwo,idim,ctxt,desc_a,info)
  call psb_s_gen_const(sy,snegativeonehalf,idim,ctxt,desc_a,info)
  call psb_s_gen_const(sz,szero,idim,ctxt,desc_a,info)
  call psb_gemlt(sone,sx,sy,szero,sz,desc_a,info)
  hasitnotfailed = psb_check_ans(sz,snegativeone,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> mlt Z = X*Y (Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- mlt Z = X*Y (Single Precision)")')
  end if

  ! Single Precision Complex
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  if (iam == psb_root_) write(psb_out_unit,'(" Complex Single Precision")')
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  ! X = 1 + 0i
  call psb_c_gen_const(cx,(1.0_psb_spk_, 0.0_psb_spk_),idim,ctxt,desc_a,info)
  hasitnotfailed = psb_check_ans(cx,(1.0_psb_spk_, 0.0_psb_spk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> Constant vector (Complex Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- Constant vector (Complex Single Precision)")')
  end if
  ! X = 2+0i , Y = -2+0i, Y = X + Y = 0+0i
  call psb_c_gen_const(cx,ctwo,idim,ctxt,desc_a,info)
  call psb_c_gen_const(cy,cnegativetwo,idim,ctxt,desc_a,info)
  call psb_geaxpby(cone,cx,cone,cy,desc_a,info)
  hasitnotfailed = psb_check_ans(cy,(0.0_psb_spk_, 0.0_psb_spk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = X + Y (Complex Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = X + Y (Complex Single Precision)")')
  end if
  ! X = 1 , Y = -2, Y = X + Y = 1 -2 = -1
  call psb_c_gen_const(cx,(1.0_psb_spk_, 0.0_psb_spk_),idim,ctxt,desc_a,info)
  call psb_c_gen_const(cy,cnegativetwo,idim,ctxt,desc_a,info)
  call psb_geaxpby(cone,cx,cone,cy,desc_a,info)
  hasitnotfailed = psb_check_ans(cy,(-1.0_psb_spk_, 0.0_psb_spk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = X + Y (Complex Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = X + Y (Complex Single Precision)")')
  end if
  ! X = 1 , Y =  2, Y = -X + Y = -1 +2 = 1
  call psb_c_gen_const(cx,(1.0_psb_spk_, 0.0_psb_spk_),idim,ctxt,desc_a,info)
  call psb_c_gen_const(cy,(2.0_psb_spk_, 0.0_psb_spk_),idim,ctxt,desc_a,info)
  call psb_geaxpby(cnegativeone,cx,cone,cy,desc_a,info)
  hasitnotfailed = psb_check_ans(cy,(1.0_psb_spk_, 0.0_psb_spk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = -X + Y (Complex Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = -X + Y (Complex Single Precision)")')
  end if
  ! X = 2 , Y =  -2, Y = 0.5*X + Y = 1 - 2 = -1
  call psb_c_gen_const(cx,ctwo,idim,ctxt,desc_a,info)
  call psb_c_gen_const(cy,cnegativetwo,idim,ctxt,desc_a,info)
  call psb_geaxpby(conehalf,cx,cone,cy,desc_a,info)
  hasitnotfailed = psb_check_ans(cy,(-1.0_psb_spk_, 0.0_psb_spk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = 0.5 X + Y (Complex Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = 0.5 X + Y (Complex Single Precision)")')
  end if
  ! X = -2 , Y =  1, Z = 0, Z = X + Y = -2 + 1 = -1
  call psb_c_gen_const(cx,cnegativetwo,idim,ctxt,desc_a,info)
  call psb_c_gen_const(cy,(1.0_psb_spk_, 0.0_psb_spk_),idim,ctxt,desc_a,info)
  call psb_c_gen_const(cz,(0.0_psb_spk_, 0.0_psb_spk_),idim,ctxt,desc_a,info)
  call psb_geaxpby(cone,cx,cone,cy,cz,desc_a,info)
  hasitnotfailed = psb_check_ans(cz,(-1.0_psb_spk_, 0.0_psb_spk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Z = X + Y (Complex Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Z = X + Y (Complex Single Precision)")')
  end if
  ! X = 2 , Y =  1, Z = 0, Z = X - Y = 2 - 1 = 1
  call psb_c_gen_const(cx,ctwo,idim,ctxt,desc_a,info)
  call psb_c_gen_const(cy,(1.0_psb_spk_, 0.0_psb_spk_),idim,ctxt,desc_a,info)
  call psb_c_gen_const(cz,(0.0_psb_spk_, 0.0_psb_spk_),idim,ctxt,desc_a,info)
  call psb_geaxpby(cone,cx,cnegativeone,cy,cz,desc_a,info)
  hasitnotfailed = psb_check_ans(cz,(1.0_psb_spk_, 0.0_psb_spk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Z = X - Y (Complex Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Z = X - Y (Complex Single Precision)")')
  end if
   ! X = 2 , Y =  -0.5, Z = 0, Z = X*Y = 2*(-0.5) = -1
  call psb_c_gen_const(cx,ctwo,idim,ctxt,desc_a,info)
  call psb_c_gen_const(cy,(snegativeonehalf,0.0_psb_spk_),idim,ctxt,desc_a,info)
  call psb_c_gen_const(cz,(0.0_psb_spk_, 0.0_psb_spk_),idim,ctxt,desc_a,info)
  call psb_gemlt(cone,cx,cy,(0.0_psb_spk_, 0.0_psb_spk_),cz,desc_a,info)
  hasitnotfailed = psb_check_ans(cz,(-1.0_psb_spk_, 0.0_psb_spk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> mlt Z = X*Y (Complex Single Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- mlt Z = X*Y (Complex Single Precision)")')
  end if

  ! Double Precision Complex
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  if (iam == psb_root_) write(psb_out_unit,'(" Complex Double Precision")')
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  ! X = 1 + 0i
  call psb_z_gen_const(zx,(1.0_psb_dpk_, 0.0_psb_dpk_),idim,ctxt,desc_a,info)
  hasitnotfailed = psb_check_ans(zx,(1.0_psb_dpk_, 0.0_psb_dpk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> Constant vector (Complex Double Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- Constant vector (Complex Double Precision)")')
  end if
  ! X = 2+0i , Y = -2+0i, Y = X + Y = 0+0i
  call psb_z_gen_const(zx,ztwo,idim,ctxt,desc_a,info)
  call psb_z_gen_const(zy,znegativetwo,idim,ctxt,desc_a,info)
  call psb_geaxpby(zone,zx,zone,zy,desc_a,info)
  hasitnotfailed = psb_check_ans(zy,(0.0_psb_dpk_, 0.0_psb_dpk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = X + Y (Complex Double Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = X + Y (Complex Double Precision)")')
  end if
  ! X = 1 , Y = -2, Y = X + Y = 1 -2 = -1
  call psb_z_gen_const(zx,(1.0_psb_dpk_, 0.0_psb_dpk_),idim,ctxt,desc_a,info)
  call psb_z_gen_const(zy,znegativetwo,idim,ctxt,desc_a,info)
  call psb_geaxpby(zone,zx,zone,zy,desc_a,info)
  hasitnotfailed = psb_check_ans(zy,(-1.0_psb_dpk_, 0.0_psb_dpk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = X + Y (Complex Double Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = X + Y (Complex Double Precision)")')
  end if
  ! X = 1 , Y =  2, Y = -X + Y = -1 +2 = 1
  call psb_z_gen_const(zx,(1.0_psb_dpk_, 0.0_psb_dpk_),idim,ctxt,desc_a,info)
  call psb_z_gen_const(zy,(2.0_psb_dpk_, 0.0_psb_dpk_),idim,ctxt,desc_a,info)
  call psb_geaxpby(znegativeone,zx,zone,zy,desc_a,info)
  hasitnotfailed = psb_check_ans(zy,(1.0_psb_dpk_, 0.0_psb_dpk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = -X + Y (Complex Double Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = -X + Y (Complex Double Precision)")')
  end if
  ! X = 2 , Y =  -2, Y = 0.5*X + Y = 1 - 2 = -1
  call psb_z_gen_const(zx,ztwo,idim,ctxt,desc_a,info)
  call psb_z_gen_const(zy,znegativetwo,idim,ctxt,desc_a,info)
  call psb_geaxpby(zonehalf,zx,zone,zy,desc_a,info)
  hasitnotfailed = psb_check_ans(zy,(-1.0_psb_dpk_, 0.0_psb_dpk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Y = 0.5 X + Y (Complex Double Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Y = 0.5 X + Y (Complex Double Precision)")')
  end if
  ! X = -2 , Y =  1, Z = 0, Z = X + Y = -2 + 1 = -1
  call psb_z_gen_const(zx,znegativetwo,idim,ctxt,desc_a,info)
  call psb_z_gen_const(zy,(1.0_psb_dpk_, 0.0_psb_dpk_),idim,ctxt,desc_a,info)
  call psb_z_gen_const(zz,(0.0_psb_dpk_, 0.0_psb_dpk_),idim,ctxt,desc_a,info)
  call psb_geaxpby(zone,zx,zone,zy,zz,desc_a,info)
  hasitnotfailed = psb_check_ans(zz,(-1.0_psb_dpk_, 0.0_psb_dpk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Z = X + Y (Complex Double Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Z = X + Y (Complex Double Precision)")')
  end if
  ! X = 2 , Y =  1, Z = 0, Z = X - Y = 2 - 1 = 1
  call psb_z_gen_const(zx,ztwo,idim,ctxt,desc_a,info)
  call psb_z_gen_const(zy,(1.0_psb_dpk_, 0.0_psb_dpk_),idim,ctxt,desc_a,info)
  call psb_z_gen_const(zz,(0.0_psb_dpk_, 0.0_psb_dpk_),idim,ctxt,desc_a,info)
  call psb_geaxpby(zone,zx,znegativeone,zy,zz,desc_a,info)
  hasitnotfailed = psb_check_ans(zz,(1.0_psb_dpk_, 0.0_psb_dpk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> axpby Z = X - Y (Complex Double Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- axpby Z = X - Y (Complex Double Precision)")')
  end if
   ! X = 2 , Y =  -0.5, Z = 0, Z = X*Y = 2*(-0.5) = -1
  call psb_z_gen_const(zx,ztwo,idim,ctxt,desc_a,info)
  call psb_z_gen_const(zy,(snegativeonehalf,0.0_psb_dpk_),idim,ctxt,desc_a,info)
  call psb_z_gen_const(zz,(0.0_psb_dpk_, 0.0_psb_dpk_),idim,ctxt,desc_a,info)
  call psb_gemlt(zone,zx,zy,(0.0_psb_dpk_, 0.0_psb_dpk_),zz,desc_a,info)
  hasitnotfailed = psb_check_ans(zz,(-1.0_psb_dpk_, 0.0_psb_dpk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> mlt Z = X*Y (Complex Double Precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- mlt Z = X*Y (Complex Double Precision)")')
  end if

  !
  ! Vector to field operation
  !
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  if (iam == psb_root_) write(psb_out_unit,'("Vector to Field Operations")')
  if (iam == psb_root_) write(psb_out_unit,'(" ")')

  ! Dot product (double real)
  call psb_d_gen_const(x,two,idim,ctxt,desc_a,info)
  call psb_d_gen_const(y,onehalf,idim,ctxt,desc_a,info)
  ans = psb_gedot(x,y,desc_a,info)
  if (iam == psb_root_) then
    if(ans == np*idim) write(psb_out_unit,'("TEST PASSED >>> Dot product (double real)")')
    if(ans /= np*idim) write(psb_out_unit,'("TEST FAILED --- Dot product (double real)")')
  end if
  ! Dot product (single precision real)
  call psb_s_gen_const(sx,stwo,idim,ctxt,desc_a,info)
  call psb_s_gen_const(sy,sonehalf,idim,ctxt,desc_a,info)
  ans = psb_gedot(sx,sy,desc_a,info)
  if (iam == psb_root_) then
    if(ans == np*idim) write(psb_out_unit,'("TEST PASSED >>> Dot product (single precision real)")')
    if(ans /= np*idim) write(psb_out_unit,'("TEST FAILED --- Dot product (single precision real)")')
  end if
  ! Dot product (single precision complex)
  call psb_c_gen_const(cx,(2.0_psb_spk_, 0.0_psb_spk_),idim,ctxt,desc_a,info)
  call psb_c_gen_const(cy,(0.5_psb_spk_, 0.0_psb_spk_),idim,ctxt,desc_a,info)
  ans = psb_gedot(cx,cy,desc_a,info)
  if (iam == psb_root_) then
    if(ans == np*idim) write(psb_out_unit,'("TEST PASSED >>> Dot product (single precision complex)")')
    if(ans /= np*idim) write(psb_out_unit,'("TEST FAILED --- Dot product (single precision complex)")')
  end if
  ! Dot product (double precision complex)
  call psb_z_gen_const(zx,(2.0_psb_dpk_, 0.0_psb_dpk_),idim,ctxt,desc_a,info)
  call psb_z_gen_const(zy,(0.5_psb_dpk_, 0.0_psb_dpk_),idim,ctxt,desc_a,info)
  ans = psb_gedot(zx,zy,desc_a,info)
  if (iam == psb_root_) then
    if(ans == np*idim) write(psb_out_unit,'("TEST PASSED >>> Dot product (double precision complex)")')
    if(ans /= np*idim) write(psb_out_unit,'("TEST FAILED --- Dot product (double precision complex)")')
  end if
  ! MaxNorm (double real)
  call psb_d_gen_const(x,negativeonehalf,idim,ctxt,desc_a,info)
  ans = psb_geamax(x,desc_a,info)
  if (iam == psb_root_) then
    if(ans == onehalf) write(psb_out_unit,'("TEST PASSED >>> MaxNorm (double real)")')
    if(ans /= onehalf) write(psb_out_unit,'("TEST FAILED --- MaxNorm (double real)")')
  end if
  ! MaxNorm (single real)
  call psb_s_gen_const(sx,snegativeonehalf,idim,ctxt,desc_a,info)
  ans = psb_geamax(sx,desc_a,info)
  if (iam == psb_root_) then
    if(ans == onehalf) write(psb_out_unit,'("TEST PASSED >>> MaxNorm (single precision real)")')
    if(ans /= onehalf) write(psb_out_unit,'("TEST FAILED --- MaxNorm (single precision real)")')
  end if
  ! MaxNorm (single complex)
  call psb_c_gen_const(cx,(snegativeonehalf,0.0_psb_spk_),idim,ctxt,desc_a,info)
  ans = psb_geamax(cx,desc_a,info)
  if (iam == psb_root_) then
    if(ans == onehalf) write(psb_out_unit,'("TEST PASSED >>> MaxNorm (single precision complex)")')
    if(ans /= onehalf) write(psb_out_unit,'("TEST FAILED --- MaxNorm (single precision complex)")')
  end if
  ! MaxNorm (double complex)
  call psb_z_gen_const(zx,(snegativeonehalf,0.0_psb_dpk_),idim,ctxt,desc_a,info)
  ans = psb_geamax(zx,desc_a,info)
  if (iam == psb_root_) then
    if(ans == onehalf) write(psb_out_unit,'("TEST PASSED >>> MaxNorm (double precision complex)")')
    if(ans /= onehalf) write(psb_out_unit,'("TEST FAILED --- MaxNorm (double precision complex)")')
  end if



  !
  ! Test of multivector operation
  !
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  if (iam == psb_root_) write(psb_out_unit,'("---- Multivector Operations ----")')
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  ! X = 1
  call psb_d_gen_const_multi(mv1,one,idim,nmv,ctxt,desc_a,info)
  hasitnotfailed = psb_check_ans(mv1,one,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> Constant multivector ")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- Constant multivector ")')
  end if
  ! X = 1 (single precision)
  call psb_s_gen_const_multi(smv1,sone,idim,nmv,ctxt,desc_a,info)
  hasitnotfailed = psb_check_ans(smv1,sone,ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> Constant multivector (single precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- Constant multivector (single precision)")')
  end if
  ! X = 1 (complex single precision)
  call psb_c_gen_const_multi(cmv1,(1.0_psb_spk_, 0.0_psb_spk_),idim,nmv,ctxt,desc_a,info)
  hasitnotfailed = psb_check_ans(cmv1,(1.0_psb_spk_, 0.0_psb_spk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> Constant multivector (complex single precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- Constant multivector (complex single precision)")')
  end if
  ! X = 1 (complex double precision)
  call psb_z_gen_const_multi(zmv1,(1.0_psb_dpk_, 0.0_psb_dpk_),idim,nmv,ctxt,desc_a,info)
  hasitnotfailed = psb_check_ans(zmv1,(1.0_psb_dpk_, 0.0_psb_dpk_),ctxt)
  if (iam == psb_root_) then
    if(hasitnotfailed) write(psb_out_unit,'("TEST PASSED >>> Constant multivector (complex double precision)")')
    if(.not.hasitnotfailed) write(psb_out_unit,'("TEST FAILED --- Constant multivector (complex double precision)")')
  end if

  ! 
  ! Multivector to field operation
  !
  if (iam == psb_root_) write(psb_out_unit,'(" ")')
  if (iam == psb_root_) write(psb_out_unit,'("Multivector to Field Operations")')
  if (iam == psb_root_) write(psb_out_unit,'(" ")')

  ! Dot product: multivector vs multivector (double real)
  call psb_d_gen_const_multi(mv1,two,idim,nmv,ctxt,desc_a,info)
  call psb_d_gen_const_multi(mv2,onehalf,idim,nmv,ctxt,desc_a,info)
  ansmv = psb_gedot(mv1,mv2,desc_a,info)
  if (iam == psb_root_) then
    ! write ansmv to check
    if(all(ansmv(:) == np*idim)) write(psb_out_unit,'("TEST PASSED >>> Dot product (mv vs mv) (double real)")')
    if(any(ansmv(:) /= np*idim)) write(psb_out_unit,'("TEST FAILED --- Dot product (mv vs mv) (double real)")')
  end if
  ! Dot product: multivector vs multivector (single real)
  call psb_s_gen_const_multi(smv1,stwo,idim,nmv,ctxt,desc_a,info)
  call psb_s_gen_const_multi(smv2,sonehalf,idim,nmv,ctxt,desc_a,info)
  ansmv = psb_gedot(smv1,smv2,desc_a,info)
  if (iam == psb_root_) then
    ! write ansmv to check
    if(all(ansmv(:) == np*idim)) write(psb_out_unit,'("TEST PASSED >>> Dot product (mv vs mv) (single real)")')
    if(any(ansmv(:) /= np*idim)) write(psb_out_unit,'("TEST FAILED --- Dot product (mv vs mv) (single real)")')
  end if
  ! Dot product: multivector vs multivector (single complex)
  call psb_c_gen_const_multi(cmv1,(2.0_psb_spk_, 0.0_psb_spk_),idim,nmv,ctxt,desc_a,info)
  call psb_c_gen_const_multi(cmv2,(0.5_psb_spk_, 0.0_psb_spk_),idim,nmv,ctxt,desc_a,info)
  ansmv = psb_gedot(cmv1,cmv2,desc_a,info)
  if (iam == psb_root_) then
    ! write ansmv to check
    if(all(ansmv(:) == np*idim)) write(psb_out_unit,'("TEST PASSED >>> Dot product (mv vs mv) (single complex)")')
    if(any(ansmv(:) /= np*idim)) write(psb_out_unit,'("TEST FAILED --- Dot product (mv vs mv) (single complex)")')
  end if
  ! Dot product: multivector vs multivector (double complex)
  call psb_z_gen_const_multi(zmv1,(2.0_psb_dpk_, 0.0_psb_dpk_),idim,nmv,ctxt,desc_a,info)
  call psb_z_gen_const_multi(zmv2,(0.5_psb_dpk_, 0.0_psb_dpk_),idim,nmv,ctxt,desc_a,info)
  ansmv = psb_gedot(zmv1,zmv2,desc_a,info)
  if (iam == psb_root_) then
    ! write ansmv to check
    if(all(ansmv(:) == np*idim)) write(psb_out_unit,'("TEST PASSED >>> Dot product (mv vs mv) (double complex)")')
    if(any(ansmv(:) /= np*idim)) write(psb_out_unit,'("TEST FAILED --- Dot product (mv vs mv) (double complex)")')
  end if
  ! Dot product: multivector vs vector (double real)
  call psb_d_gen_const_multi(mv1,two,idim,nmv,ctxt,desc_a,info)
  call psb_d_gen_const(x,onehalf,idim,ctxt,desc_a,info)
  ansmv = psb_gedot(mv1,x,desc_a,info)
  if (iam == psb_root_) then
    ! write ansmv to check
    if(all(ansmv(:) == np*idim)) write(psb_out_unit,'("TEST PASSED >>> Dot product (mv vs vector) (double real)")')
    if(any(ansmv(:) /= np*idim)) write(psb_out_unit,'("TEST FAILED --- Dot product (mv vs vector) (double real)")')
  end if
  ! Dot product: multivector vs vector (single real)
  call psb_s_gen_const_multi(smv1,stwo,idim,nmv,ctxt,desc_a,info)
  call psb_s_gen_const(sx,sonehalf,idim,ctxt,desc_a,info)
  ansmv = psb_gedot(smv1,sx,desc_a,info)
  if (iam == psb_root_) then
    ! write ansmv to check
    if(all(ansmv(:) == np*idim)) write(psb_out_unit,'("TEST PASSED >>> Dot product (mv vs vector) (single real)")')
    if(any(ansmv(:) /= np*idim)) write(psb_out_unit,'("TEST FAILED --- Dot product (mv vs vector) (single real)")')
  end if
  ! Dot product: multivector vs vector (single complex)
  call psb_c_gen_const_multi(cmv1,(2.0_psb_spk_, 0.0_psb_spk_),idim,nmv,ctxt,desc_a,info)
  call psb_c_gen_const(cx,(0.5_psb_spk_, 0.0_psb_spk_),idim,ctxt,desc_a,info)
  ansmv = psb_gedot(cmv1,cx,desc_a,info)
  if (iam == psb_root_) then
    ! write ansmv to check
    if(all(ansmv(:) == np*idim)) write(psb_out_unit,'("TEST PASSED >>> Dot product (mv vs vector) (single complex)")')
    if(any(ansmv(:) /= np*idim)) write(psb_out_unit,'("TEST FAILED --- Dot product (mv vs vector) (single complex)")')
  end if
  ! Dot product: multivector vs vector (double complex)
  call psb_z_gen_const_multi(zmv1,(2.0_psb_dpk_, 0.0_psb_dpk_),idim,nmv,ctxt,desc_a,info)
  call psb_z_gen_const(zx,(0.5_psb_dpk_, 0.0_psb_dpk_),idim,ctxt,desc_a,info)
  ansmv = psb_gedot(zmv1,zx,desc_a,info)
  if (iam == psb_root_) then
    ! write ansmv to check
    if(all(ansmv(:) == np*idim)) write(psb_out_unit,'("TEST PASSED >>> Dot product (mv vs vector) (double complex)")')
    if(any(ansmv(:) /= np*idim)) write(psb_out_unit,'("TEST FAILED --- Dot product (mv vs vector) (double complex)")')
  end if



  call psb_gefree(x,desc_a,info)
  call psb_gefree(y,desc_a,info)
  call psb_gefree(z,desc_a,info)
  call psb_gefree(sx,desc_a,info)
  call psb_gefree(sy,desc_a,info)
  call psb_gefree(sz,desc_a,info)
  call psb_gefree(cx,desc_a,info)
  call psb_gefree(cy,desc_a,info)
  call psb_gefree(cz,desc_a,info)
  call psb_gefree(zx,desc_a,info)
  call psb_gefree(zy,desc_a,info)
  call psb_gefree(zz,desc_a,info)
  call psb_gefree(mv1,desc_a,info)
  call psb_gefree(mv2,desc_a,info)
  call psb_gefree(smv1,desc_a,info)
  call psb_gefree(smv2,desc_a,info)
  call psb_gefree(cmv1,desc_a,info)
  call psb_gefree(cmv2,desc_a,info)
  call psb_gefree(zmv1,desc_a,info)
  call psb_gefree(zmv2,desc_a,info)
  call psb_cdfree(desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='free routine'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if



  call psb_exit(ctxt)
  stop

9999 call psb_error(ctxt)

  stop
end program vecoperation
