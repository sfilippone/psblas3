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
module psb_realloc_mod
  use psb_const_mod
  implicit none

  !
  ! psb_realloc will reallocate the input array to have exactly 
  ! the size specified, possibly shortening it. 
  !
  Interface psb_realloc
    module procedure psb_reallocate1i
    module procedure psb_reallocate2i
    module procedure psb_reallocate2i1d
    module procedure psb_reallocate2i1s
    module procedure psb_reallocate1d
    module procedure psb_reallocate1s
    module procedure psb_reallocated2
    module procedure psb_reallocates2
    module procedure psb_reallocatei2
#if ! defined(LONG_INTEGERS)
    module procedure psb_reallocate1i8
    module procedure psb_reallocatei8_2
#endif
    module procedure psb_reallocate2i1z
    module procedure psb_reallocate2i1c
    module procedure psb_reallocate1z
    module procedure psb_reallocate1c
    module procedure psb_reallocatez2
    module procedure psb_reallocatec2
#if defined(LONG_INTEGERS)
    module procedure psb_reallocate1i4
    module procedure psb_reallocate1i4_i8
    module procedure psb_reallocate2i4
    module procedure psb_reallocate2i4_i8
    module procedure psb_rp1i1
    module procedure psb_rp1i2i2
    module procedure psb_ri1p2i2
    module procedure psb_rp1p2i2

    module procedure psb_rp1s1
    module procedure psb_rp1i2s2
    module procedure psb_ri1p2s2
    module procedure psb_rp1p2s2

    module procedure psb_rp1d1
    module procedure psb_rp1i2d2
    module procedure psb_ri1p2d2
    module procedure psb_rp1p2d2

    module procedure psb_rp1c1
    module procedure psb_rp1i2c2
    module procedure psb_ri1p2c2
    module procedure psb_rp1p2c2

    module procedure psb_rp1z1
    module procedure psb_rp1i2z2
    module procedure psb_ri1p2z2
    module procedure psb_rp1p2z2

#endif
  end Interface psb_realloc

  interface psb_move_alloc
    module procedure psb_smove_alloc1d
    module procedure psb_smove_alloc2d
    module procedure psb_dmove_alloc1d
    module procedure psb_dmove_alloc2d
    module procedure psb_imove_alloc1d
    module procedure psb_imove_alloc2d
#if !defined(LONG_INTEGERS)
    module procedure psb_i8move_alloc1d
    module procedure psb_i8move_alloc2d
#else
    module procedure psb_i4move_alloc1d
    module procedure psb_i4move_alloc2d
    module procedure psb_i4move_alloc1d_i8
    module procedure psb_i4move_alloc2d_i8
#endif
    module procedure psb_cmove_alloc1d
    module procedure psb_cmove_alloc2d
    module procedure psb_zmove_alloc1d
    module procedure psb_zmove_alloc2d
  end interface psb_move_alloc

  Interface psb_safe_ab_cpy
    module procedure psb_i_ab_cpy1d,psb_i_ab_cpy2d, &
         & psb_s_ab_cpy1d, psb_s_ab_cpy2d,&
         & psb_c_ab_cpy1d, psb_c_ab_cpy2d,&
         & psb_d_ab_cpy1d, psb_d_ab_cpy2d,&
         & psb_z_ab_cpy1d, psb_z_ab_cpy2d
  end Interface psb_safe_ab_cpy

  Interface psb_safe_cpy
    module procedure psb_i_cpy1d,psb_i_cpy2d, &
         & psb_s_cpy1d, psb_s_cpy2d,&
         & psb_c_cpy1d, psb_c_cpy2d,&
         & psb_d_cpy1d, psb_d_cpy2d,&
         & psb_z_cpy1d, psb_z_cpy2d
  end Interface psb_safe_cpy

  !
  ! psb_ensure_size will reallocate the input array if necessary
  ! to guarantee that its size is at least as large as the 
  ! value required, usually with some room to spare.
  !
  interface psb_ensure_size
    module procedure psb_icksz1d,&
#if !defined(LONG_INTEGERS)
         & psb_i8cksz1d, &
#endif
         & psb_scksz1d, psb_ccksz1d, &
         & psb_dcksz1d, psb_zcksz1d
  end Interface psb_ensure_size

  interface psb_size
    module procedure psb_isize1d, psb_isize2d,&
#if !defined(LONG_INTEGERS)
         & psb_i8size1d, psb_i8size2d,&
#endif
         & psb_ssize1d, psb_ssize2d,&
         & psb_csize1d, psb_csize2d,&
         & psb_dsize1d, psb_dsize2d,&
         & psb_zsize1d, psb_zsize2d
  end interface psb_size

  logical, private :: do_maybe_free_buffer = .true.

Contains
  
  function psb_get_maybe_free_buffer() result(res)
    logical :: res
    res = do_maybe_free_buffer
  end function psb_get_maybe_free_buffer

  subroutine psb_set_maybe_free_buffer(val)
    logical, intent(in) :: val
    do_maybe_free_buffer = val
  end subroutine psb_set_maybe_free_buffer
  
  subroutine psb_i_ab_cpy1d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),allocatable, intent(in)  :: vin(:)
    integer(psb_ipk_),allocatable, intent(out) :: vout(:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz,err_act,lb
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_ab_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_

    if (psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (allocated(vin)) then 
      isz = size(vin)
      lb  = lbound(vin,1)
      call psb_realloc(isz,vout,info,lb=lb)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:) = vin(:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_i_ab_cpy1d

  subroutine psb_i_ab_cpy2d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_), allocatable, intent(in)  :: vin(:,:)
    integer(psb_ipk_), allocatable, intent(out) :: vout(:,:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz1, isz2,err_act, lb1, lb2
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_ab_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_

    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if
    if (allocated(vin)) then 
      isz1 = size(vin,1)
      isz2 = size(vin,2)
      lb1  = lbound(vin,1)
      lb2  = lbound(vin,2)
      call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:,:) = vin(:,:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_i_ab_cpy2d

  subroutine psb_s_ab_cpy1d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    real(psb_spk_), allocatable, intent(in)  :: vin(:)
    real(psb_spk_), allocatable, intent(out) :: vout(:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz,err_act,lb
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_ab_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (allocated(vin)) then 
      isz = size(vin)
      lb  = lbound(vin,1)
      call psb_realloc(isz,vout,info,lb=lb)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:) = vin(:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_s_ab_cpy1d

  subroutine psb_s_ab_cpy2d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    real(psb_spk_), allocatable, intent(in)  :: vin(:,:)
    real(psb_spk_), allocatable, intent(out) :: vout(:,:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz1, isz2,err_act, lb1, lb2 
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_ab_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (allocated(vin)) then 
      isz1 = size(vin,1)
      isz2 = size(vin,2)
      lb1  = lbound(vin,1)
      lb2  = lbound(vin,2)
      call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:,:) = vin(:,:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_s_ab_cpy2d

  subroutine psb_d_ab_cpy1d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    real(psb_dpk_), allocatable, intent(in)  :: vin(:)
    real(psb_dpk_), allocatable, intent(out) :: vout(:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz,err_act,lb
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_ab_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (allocated(vin)) then 
      isz = size(vin)
      lb  = lbound(vin,1)
      call psb_realloc(isz,vout,info,lb=lb)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:) = vin(:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_d_ab_cpy1d

  subroutine psb_d_ab_cpy2d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    real(psb_dpk_), allocatable, intent(in)  :: vin(:,:)
    real(psb_dpk_), allocatable, intent(out) :: vout(:,:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz1, isz2,err_act, lb1, lb2 
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_ab_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (allocated(vin)) then 
      isz1 = size(vin,1)
      isz2 = size(vin,2)
      lb1  = lbound(vin,1)
      lb2  = lbound(vin,2)
      call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:,:) = vin(:,:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_d_ab_cpy2d

  subroutine psb_c_ab_cpy1d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    complex(psb_spk_), allocatable, intent(in)  :: vin(:)
    complex(psb_spk_), allocatable, intent(out) :: vout(:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz,err_act,lb
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_ab_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (allocated(vin)) then 
      isz = size(vin)
      lb  = lbound(vin,1)
      call psb_realloc(isz,vout,info,lb=lb)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:) = vin(:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_c_ab_cpy1d

  subroutine psb_c_ab_cpy2d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    complex(psb_spk_), allocatable, intent(in)  :: vin(:,:)
    complex(psb_spk_), allocatable, intent(out) :: vout(:,:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz1, isz2,err_act, lb1, lb2 
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_ab_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (allocated(vin)) then 
      isz1 = size(vin,1)
      isz2 = size(vin,2)
      lb1  = lbound(vin,1)
      lb2  = lbound(vin,2)
      call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:,:) = vin(:,:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_c_ab_cpy2d

  subroutine psb_z_ab_cpy1d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    complex(psb_dpk_), allocatable, intent(in)  :: vin(:)
    complex(psb_dpk_), allocatable, intent(out) :: vout(:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz,err_act,lb
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_ab_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if
    if (allocated(vin)) then 
      isz = size(vin)
      lb  = lbound(vin,1)
      call psb_realloc(isz,vout,info,lb=lb)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:) = vin(:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_z_ab_cpy1d

  subroutine psb_z_ab_cpy2d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    complex(psb_dpk_), allocatable, intent(in)  :: vin(:,:)
    complex(psb_dpk_), allocatable, intent(out) :: vout(:,:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz1, isz2,err_act, lb1, lb2 
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_ab_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if
    if (allocated(vin)) then 
      isz1 = size(vin,1)
      isz2 = size(vin,2)
      lb1  = lbound(vin,1)
      lb2  = lbound(vin,2)
      call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
      if (info /= psb_success_) then     
        info=psb_err_from_subroutine_
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
      else
        vout(:,:) = vin(:,:)
      endif
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_z_ab_cpy2d


  subroutine psb_i_cpy1d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_), intent(in)               :: vin(:)
    integer(psb_ipk_), allocatable, intent(out) :: vout(:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz,err_act,lb
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if
    isz = size(vin)
    lb  = lbound(vin,1)
    call psb_realloc(isz,vout,info,lb=lb)
    if (info /= psb_success_) then     
      info=psb_err_from_subroutine_
      char_err='psb_realloc'
      call psb_errpush(info,name,a_err=char_err)
      goto 9999
    else
      vout(:) = vin(:)
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_i_cpy1d

  subroutine psb_i_cpy2d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_), intent(in)               :: vin(:,:)
    integer(psb_ipk_), allocatable, intent(out) :: vout(:,:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz1, isz2,err_act, lb1, lb2
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    isz1 = size(vin,1)
    isz2 = size(vin,2)
    lb1  = lbound(vin,1)
    lb2  = lbound(vin,2)
    call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
    if (info /= psb_success_) then     
      info=psb_err_from_subroutine_
      char_err='psb_realloc'
      call psb_errpush(info,name,a_err=char_err)
      goto 9999
    else
      vout(:,:) = vin(:,:)
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_i_cpy2d

  subroutine psb_s_cpy1d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    real(psb_spk_), intent(in)               :: vin(:)
    real(psb_spk_), allocatable, intent(out) :: vout(:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz,err_act,lb
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    isz = size(vin)
    lb  = lbound(vin,1)
    call psb_realloc(isz,vout,info,lb=lb)
    if (info /= psb_success_) then     
      info=psb_err_from_subroutine_
      char_err='psb_realloc'
      call psb_errpush(info,name,a_err=char_err)
      goto 9999
    else
      vout(:) = vin(:)
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_s_cpy1d

  subroutine psb_s_cpy2d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    real(psb_spk_), intent(in)               :: vin(:,:)
    real(psb_spk_), allocatable, intent(out) :: vout(:,:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz1, isz2,err_act, lb1, lb2 
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    isz1 = size(vin,1)
    isz2 = size(vin,2)
    lb1  = lbound(vin,1)
    lb2  = lbound(vin,2)
    call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
    if (info /= psb_success_) then     
      info=psb_err_from_subroutine_
      char_err='psb_realloc'
      call psb_errpush(info,name,a_err=char_err)
      goto 9999
    else
      vout(:,:) = vin(:,:)
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_s_cpy2d

  subroutine psb_d_cpy1d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    real(psb_dpk_), intent(in)               :: vin(:)
    real(psb_dpk_), allocatable, intent(out) :: vout(:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz,err_act,lb
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if
    isz = size(vin)
    lb  = lbound(vin,1)
    call psb_realloc(isz,vout,info,lb=lb)
    if (info /= psb_success_) then     
      info=psb_err_from_subroutine_
      char_err='psb_realloc'
      call psb_errpush(info,name,a_err=char_err)
      goto 9999
    else
      vout(:) = vin(:)
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_d_cpy1d

  subroutine psb_d_cpy2d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    real(psb_dpk_), intent(in)               :: vin(:,:)
    real(psb_dpk_), allocatable, intent(out) :: vout(:,:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz1, isz2,err_act, lb1, lb2 
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    isz1 = size(vin,1)
    isz2 = size(vin,2)
    lb1  = lbound(vin,1)
    lb2  = lbound(vin,2)
    call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
    if (info /= psb_success_) then     
      info=psb_err_from_subroutine_
      char_err='psb_realloc'
      call psb_errpush(info,name,a_err=char_err)
      goto 9999
    else
      vout(:,:) = vin(:,:)
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_d_cpy2d

  subroutine psb_c_cpy1d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    complex(psb_spk_), intent(in)               :: vin(:)
    complex(psb_spk_), allocatable, intent(out) :: vout(:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz,err_act,lb
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    isz = size(vin)
    lb  = lbound(vin,1)
    call psb_realloc(isz,vout,info,lb=lb)
    if (info /= psb_success_) then     
      info=psb_err_from_subroutine_
      char_err='psb_realloc'
      call psb_errpush(info,name,a_err=char_err)
      goto 9999
    else
      vout(:) = vin(:)
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_c_cpy1d

  subroutine psb_c_cpy2d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    complex(psb_spk_), intent(in)               :: vin(:,:)
    complex(psb_spk_), allocatable, intent(out) :: vout(:,:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz1, isz2,err_act, lb1, lb2 
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    isz1 = size(vin,1)
    isz2 = size(vin,2)
    lb1  = lbound(vin,1)
    lb2  = lbound(vin,2)
    call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
    if (info /= psb_success_) then     
      info=psb_err_from_subroutine_
      char_err='psb_realloc'
      call psb_errpush(info,name,a_err=char_err)
      goto 9999
    else
      vout(:,:) = vin(:,:)
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_c_cpy2d

  subroutine psb_z_cpy1d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    complex(psb_dpk_), intent(in)               :: vin(:)
    complex(psb_dpk_), allocatable, intent(out) :: vout(:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz,err_act,lb
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_cpy'
    call psb_erractionsave(err_act)
    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if
    isz = size(vin)
    lb  = lbound(vin,1)
    call psb_realloc(isz,vout,info,lb=lb)
    if (info /= psb_success_) then     
      info=psb_err_from_subroutine_
      char_err='psb_realloc'
      call psb_errpush(info,name,a_err=char_err)
      goto 9999
    else
      vout(:) = vin(:)
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_z_cpy1d

  subroutine psb_z_cpy2d(vin,vout,info) 
    use psb_error_mod

    ! ...Subroutine Arguments  
    complex(psb_dpk_), intent(in)               :: vin(:,:)
    complex(psb_dpk_), allocatable, intent(out) :: vout(:,:)
    integer(psb_ipk_) :: info
    ! ...Local Variables

    integer(psb_ipk_) :: isz1, isz2,err_act, lb1, lb2 
    character(len=20)  :: name, char_err
    logical, parameter :: debug=.false.

    name='psb_safe_cpy'
    call psb_erractionsave(err_act)

    info=psb_success_
    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if
    isz1 = size(vin,1)
    isz2 = size(vin,2)
    lb1  = lbound(vin,1)
    lb2  = lbound(vin,2)
    call psb_realloc(isz1,isz2,vout,info,lb1=lb1,lb2=lb2)
    if (info /= psb_success_) then     
      info=psb_err_from_subroutine_
      char_err='psb_realloc'
      call psb_errpush(info,name,a_err=char_err)
      goto 9999
    else
      vout(:,:) = vin(:,:)
    endif

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  end subroutine psb_z_cpy2d


  function psb_isize1d(vin)
    integer(psb_ipk_) :: psb_isize1d
    integer(psb_ipk_), allocatable, intent(in) :: vin(:)

    if (.not.allocated(vin)) then 
      psb_isize1d = 0
    else
      psb_isize1d = size(vin)
    end if
  end function psb_isize1d

  function psb_isize2d(vin,dim)
    integer(psb_ipk_) :: psb_isize2d
    integer(psb_ipk_), allocatable, intent(in) :: vin(:,:)
    integer(psb_ipk_), optional :: dim
    integer(psb_ipk_) :: dim_

    if (.not.allocated(vin)) then 
      psb_isize2d = 0
    else
      if (present(dim)) then 
        dim_= dim
        psb_isize2d = size(vin,dim=dim_)
      else
        psb_isize2d = size(vin)
      end if
    end if
  end function psb_isize2d

#if !defined(LONG_INTEGERS)  
  function psb_i8size1d(vin)
    integer(psb_ipk_) :: psb_i8size1d
    integer(psb_long_int_k_), allocatable, intent(in) :: vin(:)

    if (.not.allocated(vin)) then 
      psb_i8size1d = 0
    else
      psb_i8size1d = size(vin)
    end if
  end function psb_i8size1d

  function psb_i8size2d(vin,dim)
    integer(psb_ipk_) :: psb_i8size2d
    integer(psb_long_int_k_), allocatable, intent(in) :: vin(:,:)
    integer(psb_ipk_), optional :: dim
    integer(psb_ipk_) :: dim_

    if (.not.allocated(vin)) then 
      psb_i8size2d = 0
    else
      if (present(dim)) then 
        dim_= dim
        psb_i8size2d = size(vin,dim=dim_)
      else
        psb_i8size2d = size(vin)
      end if
    end if
  end function psb_i8size2d
#endif

  function psb_ssize1d(vin)
    integer(psb_ipk_) :: psb_ssize1d
    real(psb_spk_), allocatable, intent(in) :: vin(:)

    if (.not.allocated(vin)) then 
      psb_ssize1d = 0
    else
      psb_ssize1d = size(vin)
    end if
  end function psb_ssize1d

  function psb_ssize2d(vin,dim)
    integer(psb_ipk_) :: psb_ssize2d
    real(psb_spk_), allocatable, intent(in) :: vin(:,:)
    integer(psb_ipk_), optional :: dim
    integer(psb_ipk_) :: dim_


    if (.not.allocated(vin)) then 
      psb_ssize2d = 0
    else
      if (present(dim)) then 
        dim_= dim
        psb_ssize2d = size(vin,dim=dim_)
      else
        psb_ssize2d = size(vin)
      end if
    end if
  end function psb_ssize2d

  function psb_dsize1d(vin)
    integer(psb_ipk_) :: psb_dsize1d
    real(psb_dpk_), allocatable, intent(in) :: vin(:)

    if (.not.allocated(vin)) then 
      psb_dsize1d = 0
    else
      psb_dsize1d = size(vin)
    end if
  end function psb_dsize1d

  function psb_dsize2d(vin,dim)
    integer(psb_ipk_) :: psb_dsize2d
    real(psb_dpk_), allocatable, intent(in) :: vin(:,:)
    integer(psb_ipk_), optional :: dim
    integer(psb_ipk_) :: dim_


    if (.not.allocated(vin)) then 
      psb_dsize2d = 0
    else
      if (present(dim)) then 
        dim_= dim
        psb_dsize2d = size(vin,dim=dim_)
      else
        psb_dsize2d = size(vin)
      end if
    end if
  end function psb_dsize2d


  function psb_csize1d(vin)
    integer(psb_ipk_) :: psb_csize1d
    complex(psb_spk_), allocatable, intent(in) :: vin(:)

    if (.not.allocated(vin)) then 
      psb_csize1d = 0
    else
      psb_csize1d = size(vin)
    end if
  end function psb_csize1d

  function psb_csize2d(vin,dim)
    integer(psb_ipk_) :: psb_csize2d
    complex(psb_spk_), allocatable, intent(in) :: vin(:,:)
    integer(psb_ipk_), optional :: dim
    integer(psb_ipk_) :: dim_

    if (.not.allocated(vin)) then 
      psb_csize2d = 0
    else
      if (present(dim)) then 
        dim_= dim
        psb_csize2d = size(vin,dim=dim_)
      else
        psb_csize2d = size(vin)
      end if
    end if
  end function psb_csize2d

  function psb_zsize1d(vin)
    integer(psb_ipk_) :: psb_zsize1d
    complex(psb_dpk_), allocatable, intent(in) :: vin(:)

    if (.not.allocated(vin)) then 
      psb_zsize1d = 0
    else
      psb_zsize1d = size(vin)
    end if
  end function psb_zsize1d

  function psb_zsize2d(vin,dim)
    integer(psb_ipk_) :: psb_zsize2d
    complex(psb_dpk_), allocatable, intent(in) :: vin(:,:)
    integer(psb_ipk_), optional :: dim
    integer(psb_ipk_) :: dim_

    if (.not.allocated(vin)) then 
      psb_zsize2d = 0
    else
      if (present(dim)) then 
        dim_= dim
        psb_zsize2d = size(vin,dim=dim_)
      else
        psb_zsize2d = size(vin)
      end if
    end if
  end function psb_zsize2d


  Subroutine psb_icksz1d(len,v,info,pad,addsz,newsz)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in)                 :: len
    integer(psb_ipk_),allocatable, intent(inout) :: v(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional, intent(in) :: pad
    integer(psb_ipk_), optional, intent(in) :: addsz,newsz
    ! ...Local Variables
    character(len=20)  :: name
    logical, parameter :: debug=.false.
    integer(psb_ipk_) :: isz, err_act

    name='psb_ensure_size'
    call psb_erractionsave(err_act)
    info=psb_success_

    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    If (len > psb_size(v)) Then
      if (present(newsz)) then 
        isz = (max(len+1,newsz))
      else
        if (present(addsz)) then 
          isz = len+max(1,addsz)
        else
          isz = max(len+10, int(1.25*len))
        endif
      endif
      call psb_realloc(isz,v,info,pad=pad)

      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_realloc')
        goto 9999
      end if
    end If

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return


  End Subroutine psb_icksz1d

#if !defined(LONG_INTEGERS)
  Subroutine psb_i8cksz1d(len,v,info,pad,addsz,newsz)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in)                 :: len
    Integer(psb_long_int_k_),allocatable, intent(inout) :: v(:)
    integer(psb_ipk_) :: info
    integer(psb_long_int_k_), optional, intent(in) :: pad
    integer(psb_ipk_), optional, intent(in) :: addsz,newsz
    ! ...Local Variables
    character(len=20)  :: name
    logical, parameter :: debug=.false.
    integer(psb_ipk_) :: isz, err_act

    name='psb_ensure_size'
    call psb_erractionsave(err_act)
    info=psb_success_

    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    If (len > psb_size(v)) Then
      if (present(newsz)) then 
        isz = (max(len+1,newsz))
      else
        if (present(addsz)) then 
          isz = len+max(1,addsz)
        else
          isz = max(len+10, int(1.25*len))
        endif
      endif
      call psb_realloc(isz,v,info,pad=pad)

      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_realloc')
        goto 9999
      end if
    end If

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  End Subroutine psb_i8cksz1d
#endif

  Subroutine psb_scksz1d(len,v,info,pad,addsz,newsz)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in)                 :: len
    real(psb_spk_),allocatable, intent(inout) :: v(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional, intent(in)          :: addsz,newsz
    real(psb_spk_), optional, intent(in) :: pad
    ! ...Local Variables
    character(len=20)  :: name
    logical, parameter :: debug=.false.
    integer(psb_ipk_) :: isz, err_act

    name='psb_ensure_size'
    call psb_erractionsave(err_act)
    info=psb_success_

    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    If (len > psb_size(v)) Then
      if (present(newsz)) then 
        isz = (max(len+1,newsz))
      else
        if (present(addsz)) then 
          isz = len+max(1,addsz)
        else
          isz = max(len+10, int(1.25*len))
        endif
      endif

      call psb_realloc(isz,v,info,pad=pad)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_realloc')
        goto 9999
      End If
    end If

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return


  End Subroutine psb_scksz1d

  Subroutine psb_dcksz1d(len,v,info,pad,addsz,newsz)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in)                 :: len
    real(psb_dpk_),allocatable, intent(inout) :: v(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional, intent(in)          :: addsz,newsz
    real(psb_dpk_), optional, intent(in) :: pad
    ! ...Local Variables
    character(len=20)  :: name
    logical, parameter :: debug=.false.
    integer(psb_ipk_) :: isz, err_act

    name='psb_ensure_size'
    call psb_erractionsave(err_act)
    info=psb_success_

    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    If (len > psb_size(v)) Then
      if (present(newsz)) then 
        isz = (max(len+1,newsz))
      else
        if (present(addsz)) then 
          isz = len+max(1,addsz)
        else
          isz = max(len+10, int(1.25*len))
        endif
      endif

      call psb_realloc(isz,v,info,pad=pad)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_realloc')
        goto 9999
      End If
    end If

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return


  End Subroutine psb_dcksz1d


  Subroutine psb_ccksz1d(len,v,info,pad,addsz,newsz)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in)                             :: len
    complex(psb_spk_),allocatable, intent(inout) :: v(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional, intent(in)                  :: addsz,newsz
    complex(psb_spk_), optional, intent(in)      :: pad
    ! ...Local Variables
    character(len=20)  :: name
    logical, parameter :: debug=.false.
    integer(psb_ipk_) :: isz, err_act

    name='psb_ensure_size'
    call psb_erractionsave(err_act)
    info=psb_success_

    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    If (len > psb_size(v)) Then
      if (present(newsz)) then 
        isz = (max(len+1,newsz))
      else
        if (present(addsz)) then 
          isz = len+max(1,addsz)
        else
          isz = max(len+10, int(1.25*len))
        endif
      endif
      call psb_realloc(isz,v,info,pad=pad)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_realloc')
        goto 9999
      end if
    end If

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return


  End Subroutine psb_ccksz1d


  Subroutine psb_zcksz1d(len,v,info,pad,addsz,newsz)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in)                             :: len
    complex(psb_dpk_),allocatable, intent(inout) :: v(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional, intent(in)                  :: addsz,newsz
    complex(psb_dpk_), optional, intent(in)      :: pad
    ! ...Local Variables
    character(len=20)  :: name
    logical, parameter :: debug=.false.
    integer(psb_ipk_) :: isz, err_act

    name='psb_ensure_size'
    call psb_erractionsave(err_act)
    info=psb_success_

    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    If (len > psb_size(v)) Then
      if (present(newsz)) then 
        isz = (max(len+1,newsz))
      else
        if (present(addsz)) then 
          isz = len+max(1,addsz)
        else
          isz = max(len+10, int(1.25*len))
        endif
      endif
      call psb_realloc(isz,v,info,pad=pad)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_
        call psb_errpush(info,name,a_err='psb_realloc')
        goto 9999
      end if
    end If

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return


  End Subroutine psb_zcksz1d


  Subroutine psb_reallocate1i(len,rrax,info,pad,lb)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in)                 :: len
    integer(psb_ipk_),allocatable, intent(inout) :: rrax(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional, intent(in) :: pad
    integer(psb_ipk_), optional, intent(in) :: lb
    ! ...Local Variables
    integer(psb_ipk_),allocatable  :: tmp(:)
    integer(psb_ipk_) :: dim, err_act, err,lb_, lbi, ub_
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_reallocate1i' 
    call psb_erractionsave(err_act)
    info=psb_success_

    if (debug) write(psb_err_unit,*) 'reallocate I',len
    if (psb_get_errstatus() /= 0) then 
      if (debug) write(psb_err_unit,*) 'reallocate errstatus /= 0'
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0)) then 
      err=4025
      call psb_errpush(err,name,&
           & i_err=(/len,izero,izero,izero,izero/),a_err='integer')
      goto 9999
    end if
    ub_ = lb_+len-1
    if (debug) write(psb_err_unit,*) 'reallocate : lb ub ',lb_, ub_
    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1) 
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, &
               & i_err=(/len,izero,izero,izero,izero/),a_err='integer')
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        if (debug) write(psb_err_unit,*) 'reallocate : calling move_alloc '
        call psb_move_alloc(tmp,rrax,info)
        if (debug) write(psb_err_unit,*) 'reallocate : from move_alloc ',info
      end if
    else
      dim = 0
      allocate(rrax(lb_:ub_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, &
             & i_err=(/len,izero,izero,izero,izero/),a_err='integer')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb_-1+dim+1:lb_-1+len) = pad
    endif
    if (debug) write(psb_err_unit,*) 'end reallocate : ',info
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err

    call psb_error_handler(err_act)

    return

  End Subroutine psb_reallocate1i

  Subroutine psb_reallocate1s(len,rrax,info,pad,lb)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len
    Real(psb_spk_),allocatable, intent(inout) :: rrax(:)
    integer(psb_ipk_) :: info
    real(psb_spk_), optional, intent(in) :: pad
    integer(psb_ipk_), optional, intent(in) :: lb

    ! ...Local Variables
    Real(psb_spk_),allocatable  :: tmp(:)
    integer(psb_ipk_) :: dim,err_act,err, lb_, lbi,ub_
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_reallocate1s'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (debug) write(psb_err_unit,*) 'reallocate S',len

    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0)) then 
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len,izero,izero,izero,izero/),a_err='real(psb_spk_)')
      goto 9999
    end if
    ub_ = lb_ + len-1

    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1)
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, &
               & i_err=(/len,izero,izero,izero,izero/),a_err='real(psb_spk_)')
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim = 0
      Allocate(rrax(lb_:ub_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, &
             & i_err=(/len,izero,izero,izero,izero/),a_err='real(psb_spk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb_-1+dim+1:lb_-1+len) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_reallocate1s

  Subroutine psb_reallocate1d(len,rrax,info,pad,lb)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len
    Real(psb_dpk_),allocatable, intent(inout) :: rrax(:)
    integer(psb_ipk_) :: info
    real(psb_dpk_), optional, intent(in) :: pad
    integer(psb_ipk_), optional, intent(in) :: lb

    ! ...Local Variables
    Real(psb_dpk_),allocatable  :: tmp(:)
    integer(psb_ipk_) :: dim,err_act,err, lb_, lbi,ub_
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_reallocate1d'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (debug) write(psb_err_unit,*) 'reallocate D',len

    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0)) then 
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len,izero,izero,izero,izero/),a_err='real(psb_dpk_)')
      goto 9999
    end if
    ub_ = lb_ + len-1

    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1)
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, &
               & i_err=(/len,izero,izero,izero,izero/),a_err='real(psb_dpk_)')
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim = 0
      Allocate(rrax(lb_:ub_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, &
             & i_err=(/len,izero,izero,izero,izero/),a_err='real(psb_dpk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb_-1+dim+1:lb_-1+len) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_reallocate1d


  Subroutine psb_reallocate1c(len,rrax,info,pad,lb)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len
    complex(psb_spk_),allocatable, intent(inout):: rrax(:)
    integer(psb_ipk_) :: info
    complex(psb_spk_), optional, intent(in) :: pad
    integer(psb_ipk_), optional, intent(in) :: lb

    ! ...Local Variables
    complex(psb_spk_),allocatable  :: tmp(:)
    integer(psb_ipk_) :: dim,err_act,err,lb_,ub_,lbi
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_reallocate1c'
    call psb_erractionsave(err_act)
    info=psb_success_
    if (debug) write(psb_err_unit,*) 'reallocate C',len    
    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0)) then 
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len,izero,izero,izero,izero/),a_err='complex(psb_spk_)')
      goto 9999
    end if
    ub_ = lb_+len-1

    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1) 
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, &
               & i_err=(/len,izero,izero,izero,izero/),a_err='complex(psb_spk_)')
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        call psb_move_alloc(tmp,rrax,info)
      end if
    else
      dim = 0
      Allocate(rrax(lb_:ub_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, &
             & i_err=(/len,izero,izero,izero,izero/),a_err='complex(psb_spk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb_-1+dim+1:lb_-1+len) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_reallocate1c

  Subroutine psb_reallocate1z(len,rrax,info,pad,lb)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len
    complex(psb_dpk_),allocatable, intent(inout):: rrax(:)
    integer(psb_ipk_) :: info
    complex(psb_dpk_), optional, intent(in) :: pad
    integer(psb_ipk_), optional, intent(in) :: lb

    ! ...Local Variables
    complex(psb_dpk_),allocatable  :: tmp(:)
    integer(psb_ipk_) :: dim,err_act,err,lb_,ub_,lbi
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_reallocate1z'
    call psb_erractionsave(err_act)
    info=psb_success_
    if (debug) write(psb_err_unit,*) 'reallocate Z',len    
    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0)) then 
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len,izero,izero,izero,izero/),a_err='complex(psb_dpk_)')
      goto 9999
    end if
    ub_ = lb_+len-1

    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1) 
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, &
               & i_err=(/len,izero,izero,izero,izero/),a_err='complex(psb_dpk_)')
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        call psb_move_alloc(tmp,rrax,info)
      end if
    else
      dim = 0
      Allocate(rrax(lb_:ub_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, &
             & i_err=(/len,izero,izero,izero,izero/),a_err='complex(psb_dpk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb_-1+dim+1:lb_-1+len) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_reallocate1z



  Subroutine psb_reallocates2(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len1,len2
    Real(psb_spk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    real(psb_spk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    ! ...Local Variables

    Real(psb_spk_),allocatable  :: tmp(:,:)
    integer(psb_ipk_) :: dim,err_act,err, dim2,lb1_, lb2_, ub1_, ub2_,&
         & lbi1, lbi2
    character(len=20)  :: name

    name='psb_reallocates2'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (present(lb1)) then 
      lb1_ = lb1
    else
      lb1_ = 1
    endif
    if (present(lb2)) then 
      lb2_ = lb2
    else
      lb2_ = 1
    endif
    ub1_ = lb1_ + len1 -1
    ub2_ = lb2_ + len2 -1

    if (len1 < 0) then
      err=4025 
      call psb_errpush(err,name, &
           & i_err=(/len1,izero,izero,izero,izero/),a_err='real(psb_spk_)')
      goto 9999
    end if
    if (len2 < 0) then
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len2,izero,izero,izero,izero/),a_err='real(psb_spk_)')
      goto 9999
    end if


    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, &
               & i_err=(/len1*len2,izero,izero,izero,izero/),a_err='real(psb_spk_)')
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, &
             & i_err=(/len1*len2,izero,izero,izero,izero/),a_err='real(psb_spk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb1_-1+dim+1:lb1_-1+len1,:) = pad
      rrax(lb1_:lb1_-1+dim,lb2_-1+dim2+1:lb2_-1+len2) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_reallocates2


  Subroutine psb_reallocated2(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len1,len2
    Real(psb_dpk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    real(psb_dpk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    ! ...Local Variables

    Real(psb_dpk_),allocatable  :: tmp(:,:)
    integer(psb_ipk_) :: dim,err_act,err, dim2,lb1_, lb2_, ub1_, ub2_,&
         & lbi1, lbi2
    character(len=20)  :: name

    name='psb_reallocated2'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (present(lb1)) then 
      lb1_ = lb1
    else
      lb1_ = 1
    endif
    if (present(lb2)) then 
      lb2_ = lb2
    else
      lb2_ = 1
    endif
    ub1_ = lb1_ + len1 -1
    ub2_ = lb2_ + len2 -1

    if (len1 < 0) then
      err=4025 
      call psb_errpush(err,name, &
           & i_err=(/len1,izero,izero,izero,izero/),a_err='real(psb_dpk_)')
      goto 9999
    end if
    if (len2 < 0) then
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len2,izero,izero,izero,izero/),a_err='real(psb_dpk_)')
      goto 9999
    end if


    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, &
               & i_err=(/len1*len2,izero,izero,izero,izero/),a_err='real(psb_dpk_)')
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, &
             & i_err=(/len1*len2,izero,izero,izero,izero/),a_err='real(psb_dpk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb1_-1+dim+1:lb1_-1+len1,:) = pad
      rrax(lb1_:lb1_-1+dim,lb2_-1+dim2+1:lb2_-1+len2) = pad
    endif
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_reallocated2


  Subroutine psb_reallocatec2(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len1,len2
    complex(psb_spk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    complex(psb_spk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    ! ...Local Variables

    complex(psb_spk_),allocatable  :: tmp(:,:)
    integer(psb_ipk_) :: dim,err_act,err,dim2,lb1_, lb2_, ub1_, ub2_,&
         & lbi1, lbi2
    character(len=20)  :: name

    name='psb_reallocatec2'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (present(lb1)) then 
      lb1_ = lb1
    else
      lb1_ = 1
    endif
    if (present(lb2)) then 
      lb2_ = lb2
    else
      lb2_ = 1
    endif
    ub1_ = lb1_ + len1 -1
    ub2_ = lb2_ + len2 -1

    if (len1 < 0) then
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len1,izero,izero,izero,izero/),a_err='complex(psb_spk_)')
      goto 9999
    end if
    if (len2 < 0) then
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len2,izero,izero,izero,izero/),a_err='complex(psb_spk_)')
      goto 9999
    end if


    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, &
               & i_err=(/len1*len2,izero,izero,izero,izero/),a_err='complex(psb_spk_)')
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, &
             & i_err=(/len1*len2,izero,izero,izero,izero/),a_err='complex(psb_spk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb1_-1+dim+1:lb1_-1+len1,:) = pad
      rrax(lb1_:lb1_-1+dim,lb2_-1+dim2+1:lb2_-1+len2) = pad
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_reallocatec2

  Subroutine psb_reallocatez2(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len1,len2
    complex(psb_dpk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    complex(psb_dpk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    ! ...Local Variables

    complex(psb_dpk_),allocatable  :: tmp(:,:)
    integer(psb_ipk_) :: dim,err_act,err,dim2,lb1_, lb2_, ub1_, ub2_,&
         & lbi1, lbi2
    character(len=20)  :: name

    name='psb_reallocatez2'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (present(lb1)) then 
      lb1_ = lb1
    else
      lb1_ = 1
    endif
    if (present(lb2)) then 
      lb2_ = lb2
    else
      lb2_ = 1
    endif
    ub1_ = lb1_ + len1 -1
    ub2_ = lb2_ + len2 -1

    if (len1 < 0) then
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len1,izero,izero,izero,izero/),a_err='complex(psb_dpk_)')
      goto 9999
    end if
    if (len2 < 0) then
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len2,izero,izero,izero,izero/),a_err='complex(psb_dpk_)')
      goto 9999
    end if


    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, &
               & i_err=(/len1*len2,izero,izero,izero,izero/),a_err='complex(psb_dpk_)')
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, &
             & i_err=(/len1*len2,izero,izero,izero,izero/),a_err='complex(psb_dpk_)')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb1_-1+dim+1:lb1_-1+len1,:) = pad
      rrax(lb1_:lb1_-1+dim,lb2_-1+dim2+1:lb2_-1+len2) = pad
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_reallocatez2


  Subroutine psb_reallocatei2(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len1,len2
    integer(psb_ipk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    ! ...Local Variables
    integer(psb_ipk_),allocatable  :: tmp(:,:)
    integer(psb_ipk_) :: dim,err_act,err, dim2,lb1_, lb2_, ub1_, ub2_,&
         & lbi1, lbi2
    character(len=20)  :: name

    name='psb_reallocatei2'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (present(lb1)) then 
      lb1_ = lb1
    else
      lb1_ = 1
    endif
    if (present(lb2)) then 
      lb2_ = lb2
    else
      lb2_ = 1
    endif
    ub1_ = lb1_ + len1 -1
    ub2_ = lb2_ + len2 -1

    if (len1 < 0) then
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len1,izero,izero,izero,izero/),a_err='integer')
      goto 9999
    end if
    if (len2 < 0) then
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len2,izero,izero,izero,izero/),a_err='integer')
      goto 9999
    end if

    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, &
               & i_err=(/len1*len2,izero,izero,izero,izero/),a_err='integer')
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, &
             & i_err=(/len1*len2,izero,izero,izero,izero/),a_err='integer')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb1_-1+dim+1:lb1_-1+len1,:) = pad
      rrax(lb1_:lb1_-1+dim,lb2_-1+dim2+1:lb2_-1+len2) = pad
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_reallocatei2

#if !defined(LONG_INTEGERS)

  Subroutine psb_reallocate1i8(len,rrax,info,pad,lb)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in)                 :: len
    Integer(psb_long_int_k_),allocatable, intent(inout) :: rrax(:)
    integer(psb_ipk_) :: info
    integer(psb_long_int_k_), optional, intent(in) :: pad
    integer(psb_ipk_), optional, intent(in) :: lb
    ! ...Local Variables
    Integer(psb_long_int_k_),allocatable  :: tmp(:)
    integer(psb_ipk_) :: dim, err_act, err,lb_, lbi, ub_
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_reallocate1i' 
    call psb_erractionsave(err_act)
    info=psb_success_

    if (debug) write(psb_err_unit,*) 'reallocate I',len
    if (psb_get_errstatus() /= 0) then 
      if (debug) write(psb_err_unit,*) 'reallocate errstatus /= 0'
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0)) then 
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len,izero,izero,izero,izero/),a_err='integer')
      goto 9999
    end if
    ub_ = lb_+len-1
    if (debug) write(psb_err_unit,*) 'reallocate : lb ub ',lb_, ub_
    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1) 
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, &
               & i_err=(/len,izero,izero,izero,izero/),a_err='integer')
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        if (debug) write(psb_err_unit,*) 'reallocate : calling move_alloc '
        call psb_move_alloc(tmp,rrax,info)
        if (debug) write(psb_err_unit,*) 'reallocate : from move_alloc ',info
      end if
    else
      dim = 0
      allocate(rrax(lb_:ub_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, &
             & i_err=(/len,izero,izero,izero,izero/),a_err='integer')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb_-1+dim+1:lb_-1+len) = pad
    endif
    if (debug) write(psb_err_unit,*) 'end reallocate : ',info
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return


  End Subroutine psb_reallocate1i8

  Subroutine psb_reallocatei8_2(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len1,len2
    integer(psb_long_int_k_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    integer(psb_long_int_k_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    ! ...Local Variables
    integer(psb_long_int_k_),allocatable  :: tmp(:,:)
    integer(psb_ipk_) :: dim,err_act,err, dim2,lb1_, lb2_, ub1_, ub2_,&
         & lbi1, lbi2
    character(len=20)  :: name

    name='psb_reallocatei2'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (present(lb1)) then 
      lb1_ = lb1
    else
      lb1_ = 1
    endif
    if (present(lb2)) then 
      lb2_ = lb2
    else
      lb2_ = 1
    endif
    ub1_ = lb1_ + len1 -1
    ub2_ = lb2_ + len2 -1

    if (len1 < 0) then
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len1,izero,izero,izero,izero/),a_err='integer')
      goto 9999
    end if
    if (len2 < 0) then
      err=4025
      call psb_errpush(err,name, &
           & i_err=(/len2,izero,izero,izero,izero/),a_err='integer')
      goto 9999
    end if

    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= psb_success_) then
          err=4025
          call psb_errpush(err,name, &
               & i_err=(/len1*len2,izero,izero,izero,izero/),a_err='integer')
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= psb_success_) then
        err=4025
        call psb_errpush(err,name, &
             & i_err=(/len1*len2,izero,izero,izero,izero/),a_err='integer')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb1_-1+dim+1:lb1_-1+len1,:) = pad
      rrax(lb1_:lb1_-1+dim,lb2_-1+dim2+1:lb2_-1+len2) = pad
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_reallocatei8_2
#endif

  Subroutine psb_reallocate2i(len,rrax,y,info,pad)
    use psb_error_mod
    ! ...Subroutine Arguments  

    integer(psb_ipk_),Intent(in) :: len
    integer(psb_ipk_),allocatable, intent(inout) :: rrax(:),y(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional, intent(in) :: pad
    character(len=20)  :: name
    integer(psb_ipk_) :: err_act, err

    name='psb_reallocate2i'
    call psb_erractionsave(err_act)
    info=psb_success_

    if(psb_get_errstatus() /= 0) then 
      info=psb_err_from_subroutine_
      goto 9999
    end if

    call psb_reallocate1i(len,rrax,info,pad=pad)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_reallocate1i(len,y,info,pad=pad)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return

  End Subroutine psb_reallocate2i




  Subroutine psb_reallocate2i1s(len,rrax,y,z,info)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len
    integer(psb_ipk_),allocatable, intent(inout)  :: rrax(:),y(:)
    Real(psb_spk_),allocatable, intent(inout) :: z(:)
    integer(psb_ipk_) :: info
    character(len=20)  :: name
    integer(psb_ipk_) :: err_act, err
    logical, parameter :: debug=.false.

    name='psb_reallocate2i1s'
    call psb_erractionsave(err_act)


    info=psb_success_
    call psb_realloc(len,rrax,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,y,info)    
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,z,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  End Subroutine psb_reallocate2i1s


  Subroutine psb_reallocate2i1d(len,rrax,y,z,info)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len
    integer(psb_ipk_),allocatable, intent(inout)  :: rrax(:),y(:)
    Real(psb_dpk_),allocatable, intent(inout) :: z(:)
    integer(psb_ipk_) :: info
    character(len=20)  :: name
    integer(psb_ipk_) :: err_act, err

    name='psb_reallocate2i1d'
    call psb_erractionsave(err_act)

    info=psb_success_

    call psb_realloc(len,rrax,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,y,info)    
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,z,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  End Subroutine psb_reallocate2i1d



  Subroutine psb_reallocate2i1c(len,rrax,y,z,info)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len
    integer(psb_ipk_),allocatable, intent(inout) :: rrax(:),y(:)
    complex(psb_spk_),allocatable, intent(inout) :: z(:)
    integer(psb_ipk_) :: info
    character(len=20)  :: name
    integer(psb_ipk_) :: err_act, err

    name='psb_reallocate2i1c'
    call psb_erractionsave(err_act)


    info=psb_success_
    call psb_realloc(len,rrax,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,y,info)    
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,z,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  End Subroutine psb_reallocate2i1c

  Subroutine psb_reallocate2i1z(len,rrax,y,z,info)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len
    integer(psb_ipk_),allocatable, intent(inout) :: rrax(:),y(:)
    complex(psb_dpk_),allocatable, intent(inout) :: z(:)
    integer(psb_ipk_) :: info
    character(len=20)  :: name
    integer(psb_ipk_) :: err_act, err

    name='psb_reallocate2i1z'
    call psb_erractionsave(err_act)

    info=psb_success_
    call psb_realloc(len,rrax,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,y,info)    
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_realloc(len,z,info)
    if (info /= psb_success_) then
      err=4000
      call psb_errpush(err,name)
      goto 9999
    end if
    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(err_act)

    return
  End Subroutine psb_reallocate2i1z

  Subroutine psb_smove_alloc1d(vin,vout,info)
    use psb_error_mod
    real(psb_spk_), allocatable, intent(inout) :: vin(:),vout(:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC    

    call move_alloc(vin,vout)

#else      
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return
    allocate(vout(lbound(vin,1):ubound(vin,1)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_smove_alloc1d

  Subroutine psb_smove_alloc2d(vin,vout,info)
    use psb_error_mod
    real(psb_spk_), allocatable, intent(inout) :: vin(:,:),vout(:,:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return

    allocate(vout(lbound(vin,1):ubound(vin,1),&
         & lbound(vin,2):ubound(vin,2)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_smove_alloc2d

  Subroutine psb_dmove_alloc1d(vin,vout,info)
    use psb_error_mod
    real(psb_dpk_), allocatable, intent(inout) :: vin(:),vout(:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC    

    call move_alloc(vin,vout)

#else      
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return
    allocate(vout(lbound(vin,1):ubound(vin,1)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_dmove_alloc1d

  Subroutine psb_dmove_alloc2d(vin,vout,info)
    use psb_error_mod
    real(psb_dpk_), allocatable, intent(inout) :: vin(:,:),vout(:,:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return

    allocate(vout(lbound(vin,1):ubound(vin,1),&
         & lbound(vin,2):ubound(vin,2)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_dmove_alloc2d

  Subroutine psb_cmove_alloc1d(vin,vout,info)
    use psb_error_mod
    complex(psb_spk_), allocatable, intent(inout) :: vin(:),vout(:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return
    allocate(vout(lbound(vin,1):ubound(vin,1)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_cmove_alloc1d

  Subroutine psb_cmove_alloc2d(vin,vout,info)
    use psb_error_mod
    complex(psb_spk_), allocatable, intent(inout) :: vin(:,:),vout(:,:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return

    allocate(vout(lbound(vin,1):ubound(vin,1),&
         & lbound(vin,2):ubound(vin,2)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_cmove_alloc2d

  Subroutine psb_zmove_alloc1d(vin,vout,info)
    use psb_error_mod
    complex(psb_dpk_), allocatable, intent(inout) :: vin(:),vout(:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return
    allocate(vout(lbound(vin,1):ubound(vin,1)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_zmove_alloc1d

  Subroutine psb_zmove_alloc2d(vin,vout,info)
    use psb_error_mod
    complex(psb_dpk_), allocatable, intent(inout) :: vin(:,:),vout(:,:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return

    allocate(vout(lbound(vin,1):ubound(vin,1),&
         & lbound(vin,2):ubound(vin,2)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_zmove_alloc2d

  Subroutine psb_imove_alloc1d(vin,vout,info)
    use psb_error_mod
    integer(psb_ipk_), allocatable, intent(inout) :: vin(:),vout(:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return
    allocate(vout(lbound(vin,1):ubound(vin,1)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_imove_alloc1d

  Subroutine psb_imove_alloc2d(vin,vout,info)
    use psb_error_mod
    integer(psb_ipk_), allocatable, intent(inout) :: vin(:,:),vout(:,:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return

    allocate(vout(lbound(vin,1):ubound(vin,1),&
         & lbound(vin,2):ubound(vin,2)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_imove_alloc2d

#if !defined(LONG_INTEGERS)
  Subroutine psb_i8move_alloc1d(vin,vout,info)
    use psb_error_mod
    integer(psb_long_int_k_), allocatable, intent(inout) :: vin(:),vout(:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return
    allocate(vout(lbound(vin,1):ubound(vin,1)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_i8move_alloc1d

  Subroutine psb_i8move_alloc2d(vin,vout,info)
    use psb_error_mod
    integer(psb_long_int_k_), allocatable, intent(inout) :: vin(:,:),vout(:,:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return

    allocate(vout(lbound(vin,1):ubound(vin,1),&
         & lbound(vin,2):ubound(vin,2)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_i8move_alloc2d

#else

  Subroutine psb_i4move_alloc1d(vin,vout,info)
    use psb_error_mod
    integer(psb_mpik_), allocatable, intent(inout) :: vin(:),vout(:)
    integer(psb_mpik_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return
    allocate(vout(lbound(vin,1):ubound(vin,1)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_i4move_alloc1d

  Subroutine psb_i4move_alloc1d_i8(vin,vout,info)
    use psb_error_mod
    integer(psb_mpik_), allocatable, intent(inout) :: vin(:),vout(:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return
    allocate(vout(lbound(vin,1):ubound(vin,1)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_i4move_alloc1d_i8
 
  Subroutine psb_i4move_alloc2d(vin,vout,info)
    use psb_error_mod
    integer(psb_mpik_), allocatable, intent(inout) :: vin(:,:),vout(:,:)
    integer(psb_mpik_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return

    allocate(vout(lbound(vin,1):ubound(vin,1),&
         & lbound(vin,2):ubound(vin,2)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_i4move_alloc2d
 
  Subroutine psb_i4move_alloc2d_i8(vin,vout,info)
    use psb_error_mod
    integer(psb_mpik_), allocatable, intent(inout) :: vin(:,:),vout(:,:)
    integer(psb_ipk_), intent(out) :: info 
    !
    ! 
    info=psb_success_
#ifdef HAVE_MOVE_ALLOC

    call move_alloc(vin,vout)

#else
    if (allocated(vout)) then 
      deallocate(vout,stat=info)
    end if
    if (.not.allocated(vin) ) return

    allocate(vout(lbound(vin,1):ubound(vin,1),&
         & lbound(vin,2):ubound(vin,2)),stat=info)
    if (info /= psb_success_) return
    vout = vin
    deallocate(vin,stat=info)
#endif
  end Subroutine psb_i4move_alloc2d_i8

#endif

#if defined(LONG_INTEGERS)
  Subroutine psb_reallocate1i4(len,rrax,info,pad,lb)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in)                 :: len
    Integer(psb_mpik_),allocatable, intent(inout) :: rrax(:)
    integer(psb_mpik_) :: info
    integer(psb_mpik_), optional, intent(in) :: pad
    integer(psb_mpik_), optional, intent(in) :: lb
    ! ...Local Variables
    Integer(psb_mpik_),allocatable  :: tmp(:)
    integer(psb_mpik_) :: dim, lb_, lbi, ub_
    integer(psb_ipk_)  :: err, err_act, ierr(5)
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_reallocate1i4' 
    call psb_erractionsave(err_act)
    info=psb_success_

    if (debug) write(psb_err_unit,*) 'reallocate I',len
    if (psb_get_errstatus() /= 0) then 
      if (debug) write(psb_err_unit,*) 'reallocate errstatus /= 0'
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0)) then 
      err=4025; ierr(1) = len
      call psb_errpush(err,name,i_err=ierr,a_err='integer')
      goto 9999
    end if
    ub_ = lb_+len-1
    if (debug) write(psb_err_unit,*) 'reallocate : lb ub ',lb_, ub_
    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1) 
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= psb_success_) then
          err=4025; ierr(1) = len
          call psb_errpush(err,name,i_err=ierr,a_err='integer')
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        if (debug) write(psb_err_unit,*) 'reallocate : calling move_alloc '
        call psb_move_alloc(tmp,rrax,info)
        if (debug) write(psb_err_unit,*) 'reallocate : from move_alloc ',info
      end if
    else
      dim = 0
      allocate(rrax(lb_:ub_),stat=info)
      if (info /= psb_success_) then
        err=4025; ierr(1) = len
        call psb_errpush(err,name,i_err=ierr,a_err='integer')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb_-1+dim+1:lb_-1+len) = pad
    endif
    if (debug) write(psb_err_unit,*) 'end reallocate : ',info
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return


  End Subroutine psb_reallocate1i4

  Subroutine psb_reallocate1i4_i8(len,rrax,info,pad,lb)
    use psb_error_mod

    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in)                 :: len
    Integer(psb_mpik_),allocatable, intent(inout) :: rrax(:)
    integer(psb_ipk_) :: info
    integer(psb_mpik_), optional, intent(in) :: pad
    integer(psb_ipk_), optional, intent(in) :: lb
    ! ...Local Variables
    Integer(psb_mpik_),allocatable  :: tmp(:)
    integer(psb_mpik_) :: dim, lb_, lbi, ub_, iinfo
    integer(psb_ipk_)  :: err, err_act, ierr(5)
    character(len=20)  :: name
    logical, parameter :: debug=.false.

    name='psb_reallocate1i4' 
    call psb_erractionsave(err_act)
    info=psb_success_

    if (debug) write(psb_err_unit,*) 'reallocate I',len
    if (psb_get_errstatus() /= 0) then 
      if (debug) write(psb_err_unit,*) 'reallocate errstatus /= 0'
      info=psb_err_from_subroutine_
      goto 9999
    end if

    if (present(lb)) then
      lb_ = lb
    else
      lb_ = 1
    endif
    if ((len<0)) then 
      err=4025; ierr(1) = len
      call psb_errpush(err,name,i_err=ierr,a_err='integer')
      goto 9999
    end if
    ub_ = lb_+len-1
    if (debug) write(psb_err_unit,*) 'reallocate : lb ub ',lb_, ub_
    if (allocated(rrax)) then 
      dim = size(rrax)
      lbi = lbound(rrax,1) 
      If ((dim /= len).or.(lbi /= lb_))  Then
        Allocate(tmp(lb_:ub_),stat=info)
        if (info /= psb_success_) then
          err=4025; ierr(1) = len
          call psb_errpush(err,name,i_err=ierr,a_err='integer')
          goto 9999
        end if
        tmp(lb_:lb_-1+min(len,dim))=rrax(lbi:lbi-1+min(len,dim))
        if (debug) write(psb_err_unit,*) 'reallocate : calling move_alloc '
        call psb_move_alloc(tmp,rrax,iinfo)
        if (debug) write(psb_err_unit,*) 'reallocate : from move_alloc ',iinfo
      end if
    else
      dim = 0
      allocate(rrax(lb_:ub_),stat=info)
      if (info /= psb_success_) then
        err=4025; ierr(1) = len
        call psb_errpush(err,name,i_err=ierr,a_err='integer')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb_-1+dim+1:lb_-1+len) = pad
    endif
    if (debug) write(psb_err_unit,*) 'end reallocate : ',info
    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return


  End Subroutine psb_reallocate1i4_i8

  Subroutine psb_reallocate2i4(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len1,len2
    integer(psb_mpik_),allocatable :: rrax(:,:)
    integer(psb_mpik_) :: info
    integer(psb_mpik_), optional, intent(in) :: pad
    integer(psb_mpik_),Intent(in), optional  :: lb1,lb2

    ! ...Local Variables
    integer(psb_mpik_),allocatable  :: tmp(:,:)
    integer(psb_mpik_) :: dim, dim2,lb1_, lb2_, ub1_, ub2_,&
         & lbi1, lbi2
    integer(psb_ipk_)  :: err,err_act, ierr(5)
    character(len=20)  :: name

    name='psb_reallocatei2'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (present(lb1)) then 
      lb1_ = lb1
    else
      lb1_ = 1
    endif
    if (present(lb2)) then 
      lb2_ = lb2
    else
      lb2_ = 1
    endif
    ub1_ = lb1_ + len1 -1
    ub2_ = lb2_ + len2 -1

    if (len1 < 0) then
      err=4025; ierr(1) = len1
      call psb_errpush(err,name,i_err=ierr,a_err='integer')
      goto 9999
    end if
    if (len2 < 0) then
      err=4025; ierr(1) = len2
      call psb_errpush(err,name,i_err=ierr,a_err='integer')
      goto 9999
    end if

    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= psb_success_) then
          err=4025; ierr(1) = len1*len2
          call psb_errpush(err,name,i_err=ierr,a_err='integer')
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= psb_success_) then
        err=4025; ierr(1) = len1*len2
        call psb_errpush(err,name,i_err=ierr,a_err='integer')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb1_-1+dim+1:lb1_-1+len1,:) = pad
      rrax(lb1_:lb1_-1+dim,lb2_-1+dim2+1:lb2_-1+len2) = pad
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_reallocate2i4

  Subroutine psb_reallocate2i4_i8(len1,len2,rrax,info,pad,lb1,lb2)
    use psb_error_mod
    ! ...Subroutine Arguments  
    integer(psb_ipk_),Intent(in) :: len1,len2
    integer(psb_mpik_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    integer(psb_mpik_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    ! ...Local Variables
    integer(psb_mpik_),allocatable  :: tmp(:,:)
    integer(psb_mpik_) :: dim, dim2,lb1_, lb2_, ub1_, ub2_,&
         & lbi1, lbi2
    integer(psb_ipk_)  :: err,err_act, ierr(5)
    character(len=20)  :: name

    name='psb_reallocatei2'
    call psb_erractionsave(err_act)
    info=psb_success_ 
    if (present(lb1)) then 
      lb1_ = lb1
    else
      lb1_ = 1
    endif
    if (present(lb2)) then 
      lb2_ = lb2
    else
      lb2_ = 1
    endif
    ub1_ = lb1_ + len1 -1
    ub2_ = lb2_ + len2 -1

    if (len1 < 0) then
      err=4025; ierr(1) = len1
      call psb_errpush(err,name,i_err=ierr,a_err='integer')
      goto 9999
    end if
    if (len2 < 0) then
      err=4025; ierr(1) = len2
      call psb_errpush(err,name,i_err=ierr,a_err='integer')
      goto 9999
    end if

    if (allocated(rrax)) then 
      dim  = size(rrax,1)
      lbi1 = lbound(rrax,1)
      dim2 = size(rrax,2)
      lbi2 = lbound(rrax,2)
      If ((dim /= len1).or.(dim2 /= len2).or.(lbi1 /= lb1_)&
           &  .or.(lbi2 /= lb2_)) Then
        Allocate(tmp(lb1_:ub1_,lb2_:ub2_),stat=info)
        if (info /= psb_success_) then
          err=4025; ierr(1) = len1*len2
          call psb_errpush(err,name,i_err=ierr,a_err='integer')
          goto 9999
        end if
        tmp(lb1_:lb1_-1+min(len1,dim),lb2_:lb2_-1+min(len2,dim2)) = &
             & rrax(lbi1:lbi1-1+min(len1,dim),lbi2:lbi2-1+min(len2,dim2))
        call psb_move_alloc(tmp,rrax,info)
      End If
    else
      dim  = 0
      dim2 = 0
      Allocate(rrax(lb1_:ub1_,lb2_:ub2_),stat=info)
      if (info /= psb_success_) then
        err=4025; ierr(1) = len1*len2
        call psb_errpush(err,name,i_err=ierr,a_err='integer')
        goto 9999
      end if
    endif
    if (present(pad)) then 
      rrax(lb1_-1+dim+1:lb1_-1+len1,:) = pad
      rrax(lb1_:lb1_-1+dim,lb2_-1+dim2+1:lb2_-1+len2) = pad
    endif

    call psb_erractionrestore(err_act)
    return

9999 continue
    info = err
    call psb_error_handler(err_act)
    return

  End Subroutine psb_reallocate2i4_i8


  Subroutine psb_rp1i1(len,rrax,info,pad,lb)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len
    integer(psb_ipk_),allocatable, intent(inout) :: rrax(:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional, intent(in) :: pad
    integer(psb_mpik_), optional, intent(in) :: lb

    integer(psb_ipk_) :: ilen, ilb

    ilen=len
    if (present(lb)) then 
      ilb=lb
    else
      ilb = 1
    end if
    call psb_realloc(ilen,rrax,info,lb=ilb,pad=pad)

  end Subroutine psb_rp1i1


  subroutine  psb_rp1i2i2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len1
    integer(psb_ipk_),Intent(in) :: len2
    integer(psb_ipk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len1_ = len1
    call psb_realloc(len1_,len2,rrax,info,pad,lb1,lb2)
  end subroutine psb_rp1i2i2

  subroutine  psb_ri1p2i2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len2
    integer(psb_ipk_),Intent(in) :: len1
    integer(psb_ipk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len2_ = len2
    call psb_realloc(len1,len2_,rrax,info,pad,lb1,lb2)
  end subroutine psb_ri1p2i2

  subroutine  psb_rp1p2i2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len1
    integer(psb_mpik_),Intent(in) :: len2
    integer(psb_ipk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    integer(psb_ipk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len1_ = len1
    len2_ = len2
    call psb_realloc(len1_,len2_,rrax,info,pad,lb1,lb2)
  end subroutine psb_rp1p2i2

  Subroutine psb_rp1s1(len,rrax,info,pad,lb)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len
    real(psb_spk_),allocatable, intent(inout) :: rrax(:)
    integer(psb_ipk_) :: info
    real(psb_spk_), optional, intent(in) :: pad
    integer(psb_mpik_), optional, intent(in) :: lb

    integer(psb_ipk_) :: ilen, ilb

    ilen=len
    if (present(lb)) then 
      ilb=lb
    else
      ilb = 1
    end if
    call psb_realloc(ilen,rrax,info,lb=ilb,pad=pad)

  end Subroutine psb_rp1s1

  subroutine  psb_rp1i2s2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len1
    integer(psb_ipk_),Intent(in) :: len2
    real(psb_spk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    real(psb_spk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len1_ = len1
    call psb_realloc(len1_,len2,rrax,info,pad,lb1,lb2)
  end subroutine psb_rp1i2s2

  subroutine  psb_ri1p2s2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len2
    integer(psb_ipk_),Intent(in) :: len1
    real(psb_spk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    real(psb_spk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len2_ = len2
    call psb_realloc(len1,len2_,rrax,info,pad,lb1,lb2)
  end subroutine psb_ri1p2s2

  subroutine  psb_rp1p2s2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len1
    integer(psb_mpik_),Intent(in) :: len2
    real(psb_spk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    real(psb_spk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len1_ = len1
    len2_ = len2
    call psb_realloc(len1_,len2_,rrax,info,pad,lb1,lb2)
  end subroutine psb_rp1p2s2



  Subroutine psb_rp1d1(len,rrax,info,pad,lb)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len
    Real(psb_dpk_),allocatable, intent(inout) :: rrax(:)
    integer(psb_ipk_) :: info
    real(psb_dpk_), optional, intent(in) :: pad
    integer(psb_mpik_), optional, intent(in) :: lb

    integer(psb_ipk_) :: ilen, ilb

    ilen=len
    if (present(lb)) then 
      ilb=lb
    else
      ilb = 1
    end if
    call psb_realloc(ilen,rrax,info,lb=ilb,pad=pad)

  end Subroutine psb_rp1d1


  subroutine  psb_rp1i2d2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len1
    integer(psb_ipk_),Intent(in) :: len2
    Real(psb_dpk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    real(psb_dpk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len1_ = len1
    call psb_realloc(len1_,len2,rrax,info,pad,lb1,lb2)
  end subroutine psb_rp1i2d2

  subroutine  psb_ri1p2d2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len2
    integer(psb_ipk_),Intent(in) :: len1
    Real(psb_dpk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    real(psb_dpk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len2_ = len2
    call psb_realloc(len1,len2_,rrax,info,pad,lb1,lb2)
  end subroutine psb_ri1p2d2

  subroutine  psb_rp1p2d2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len1
    integer(psb_mpik_),Intent(in) :: len2
    Real(psb_dpk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    real(psb_dpk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len1_ = len1
    len2_ = len2
    call psb_realloc(len1_,len2_,rrax,info,pad,lb1,lb2)
  end subroutine psb_rp1p2d2



  Subroutine psb_rp1c1(len,rrax,info,pad,lb)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len
    complex(psb_spk_),allocatable, intent(inout) :: rrax(:)
    integer(psb_ipk_) :: info
    complex(psb_spk_), optional, intent(in) :: pad
    integer(psb_mpik_), optional, intent(in) :: lb

    integer(psb_ipk_) :: ilen, ilb

    ilen=len
    if (present(lb)) then 
      ilb=lb
    else
      ilb = 1
    end if
    call psb_realloc(ilen,rrax,info,lb=ilb,pad=pad)

  end Subroutine psb_rp1c1

  subroutine  psb_rp1i2c2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len1
    integer(psb_ipk_),Intent(in) :: len2
    complex(psb_spk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    complex(psb_spk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len1_ = len1
    call psb_realloc(len1_,len2,rrax,info,pad,lb1,lb2)
  end subroutine psb_rp1i2c2

  subroutine  psb_ri1p2c2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len2
    integer(psb_ipk_),Intent(in) :: len1
    complex(psb_spk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    complex(psb_spk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len2_ = len2
    call psb_realloc(len1,len2_,rrax,info,pad,lb1,lb2)
  end subroutine psb_ri1p2c2

  subroutine  psb_rp1p2c2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len1
    integer(psb_mpik_),Intent(in) :: len2
    complex(psb_spk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    complex(psb_spk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len1_ = len1
    len2_ = len2
    call psb_realloc(len1_,len2_,rrax,info,pad,lb1,lb2)
  end subroutine psb_rp1p2c2


  Subroutine psb_rp1z1(len,rrax,info,pad,lb)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len
    Complex(psb_dpk_),allocatable, intent(inout) :: rrax(:)
    integer(psb_ipk_) :: info
    complex(psb_dpk_), optional, intent(in) :: pad
    integer(psb_mpik_), optional, intent(in) :: lb

    integer(psb_ipk_) :: ilen, ilb

    ilen=len
    if (present(lb)) then 
      ilb=lb
    else
      ilb = 1
    end if
    call psb_realloc(ilen,rrax,info,lb=ilb,pad=pad)

  end Subroutine psb_rp1z1

  subroutine  psb_rp1i2z2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len1
    integer(psb_ipk_),Intent(in) :: len2
    Complex(psb_dpk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    complex(psb_dpk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len1_ = len1
    call psb_realloc(len1_,len2,rrax,info,pad,lb1,lb2)
  end subroutine psb_rp1i2z2

  subroutine  psb_ri1p2z2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len2
    integer(psb_ipk_),Intent(in) :: len1
    Complex(psb_dpk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    complex(psb_dpk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len2_ = len2
    call psb_realloc(len1,len2_,rrax,info,pad,lb1,lb2)
  end subroutine psb_ri1p2z2

  subroutine  psb_rp1p2z2(len1,len2,rrax,info,pad,lb1,lb2)
    ! ...Subroutine Arguments  
    integer(psb_mpik_),Intent(in) :: len1
    integer(psb_mpik_),Intent(in) :: len2
    Complex(psb_dpk_),allocatable :: rrax(:,:)
    integer(psb_ipk_) :: info
    complex(psb_dpk_), optional, intent(in) :: pad
    integer(psb_ipk_),Intent(in), optional  :: lb1,lb2

    integer(psb_ipk_) :: len1_, len2_
    len1_ = len1
    len2_ = len2
    call psb_realloc(len1_,len2_,rrax,info,pad,lb1,lb2)
  end subroutine psb_rp1p2z2

#endif


  subroutine i_trans(a,at)
    implicit none 
    integer(psb_ipk_) :: nr,nc
    integer(psb_ipk_) :: a(:,:)
    integer(psb_ipk_), allocatable, intent(out) :: at(:,:)
    integer(psb_ipk_) :: i,j,ib, ii
    integer(psb_ipk_), parameter :: nb=32

    nr = size(a,1)
    nc = size(a,2)
    allocate(at(nc,nr))
    do i=1,nr,nb
      ib=min(nb,nr-i+1)
      do ii=i,i+ib-1
        do j=1,nc
          at(j,ii) = a(ii,j)
        end do
      end do
    end do
  end subroutine i_trans

  subroutine d_trans(a,at)
    implicit none 
    integer(psb_ipk_) :: nr,nc
    real(psb_dpk_) :: a(:,:)
    real(psb_dpk_), allocatable, intent(out) :: at(:,:)
    integer(psb_ipk_) :: i,j,ib, ii
    integer(psb_ipk_), parameter :: nb=32

    nr = size(a,1)
    nc = size(a,2)
    allocate(at(nc,nr))
    do i=1,nr,nb
      ib=min(nb,nr-i+1)
      do ii=i,i+ib-1
        do j=1,nc
          at(j,ii) = a(ii,j)
        end do
      end do
    end do
  end subroutine d_trans

end module psb_realloc_mod
