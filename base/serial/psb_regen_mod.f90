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
module psb_regen_mod

  interface csr_regen
    module procedure csr_dsp_regen, csr_zsp_regen
  end interface
  interface coo_regen
    module procedure coo_dsp_regen, coo_zsp_regen
  end interface
  interface jad_regen
    module procedure jad_dsp_regen, jad_zsp_regen
  end interface

contains
  
  subroutine csr_dsp_regen(a,info)

    use psb_spmat_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    type(psb_dspmat_type), intent(inout)    :: a
    integer                              :: info

    integer :: i, ip1,ip2,nnz,iflag,ichk,nnzt
    real(psb_dpk_), allocatable  :: work(:)
    integer :: err_act
    character(len=20)  :: name
    integer            :: debug_level, debug_unit


    name='psb_spcnv'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()


    !
    !   dupl_ and upd_ fields should not be changed. 
    !
    select case(psb_sp_getifld(psb_upd_,a,info))

    case(psb_upd_perm_)

      allocate(work(size(a%aspk)+1000),stat=info)
      if (info /= 0) then
        info=2040
        call psb_errpush(info,name)
        goto 9999
      end if

      if (debug_level >= psb_debug_serial_)&
           & write(debug_unit,*)  trim(name),'Regeneration with psb_upd_perm_'
      ip1   = psb_sp_getifld(psb_upd_pnt_,a,info) 
      ip2   = a%ia2(ip1+psb_ip2_)
      nnz   = a%ia2(ip1+psb_nnz_)
      iflag = a%ia2(ip1+psb_iflag_)
      ichk  = a%ia2(ip1+psb_ichk_)
      nnzt  = a%ia2(ip1+psb_nnzt_)
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*)  trim(name),'Regeneration start: ',&
           &   a%infoa(psb_upd_),psb_upd_perm_,nnz,nnzt ,iflag,info

      if ((ichk/=nnzt+iflag).or.(nnz/=nnzt)) then               
        info = 8889
        write(debug_unit,*) trim(name),'Regeneration start error: ',&
             &   a%infoa(psb_upd_),psb_upd_perm_,nnz,nnzt ,iflag,ichk    
        call psb_errpush(info,name)
        goto 9999
      endif
      do i= 1, nnz
        work(i) = dzero
      enddo
      select case(iflag) 
      case(psb_dupl_ovwrt_,psb_dupl_err_) 
        do i=1, nnz 
          work(a%ia2(ip2+i-1)) = a%aspk(i) 
        enddo
      case(psb_dupl_add_) 
        do i=1, nnz 
          work(a%ia2(ip2+i-1)) = a%aspk(i) + work(a%ia2(ip2+i-1)) 
        enddo
      case default
        info = 8887
        call psb_errpush(info,name)
        goto 9999          
      end select

      do i=1, nnz
        a%aspk(i) = work(i)
      enddo


    case(psb_upd_srch_)
      ! Nothing to be done  here. 
      if (debug_level >= psb_debug_serial_)&
           & write(debug_unit,*)  trim(name),&
           & 'Going through on regeneration with psb_upd_srch_'
    case default
      ! Wrong value
      info = 8888
      call psb_errpush(info,name)
      goto 9999

    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine csr_dsp_regen

  subroutine coo_dsp_regen(a,info)

    use psb_spmat_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    type(psb_dspmat_type), intent(inout)    :: a
    integer                              :: info

    integer :: i, ip1,ip2,nnz,iflag,ichk,nnzt
    real(psb_dpk_), allocatable  :: work(:)
    integer :: err_act
    character(len=20)  :: name
    integer            :: debug_level, debug_unit


    name='psb_spcnv'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()


    !
    !   dupl_ and upd_ fields should not be changed. 
    !
    select case(psb_sp_getifld(psb_upd_,a,info))

    case(psb_upd_perm_)

      allocate(work(size(a%aspk)+1000),stat=info)
      if (info /= 0) then
        info=2040
        call psb_errpush(info,name)
        goto 9999
      end if

      if (debug_level >= psb_debug_serial_)&
           & write(debug_unit,*)  trim(name),'Regeneration with psb_upd_perm_'
      ip1   = psb_sp_getifld(psb_upd_pnt_,a,info) 
      ip2   = a%ia2(ip1+psb_ip2_)
      nnz   = a%ia2(ip1+psb_nnz_)
      iflag = a%ia2(ip1+psb_iflag_)
      ichk  = a%ia2(ip1+psb_ichk_)
      nnzt  = a%ia2(ip1+psb_nnzt_)
      if (debug_level >= psb_debug_serial_)&
           & write(debug_unit,*) trim(name),'Regeneration start: ',&
           &   a%infoa(psb_upd_),psb_upd_perm_,nnz,nnzt ,iflag,info

      if ((ichk/=nnzt+iflag).or.(nnz/=nnzt)) then               
        info = 8889
        write(debug_unit,*) trim(name),'Regeneration start error: ',&
             &   a%infoa(psb_upd_),psb_upd_perm_,nnz,nnzt ,iflag,ichk    
        call psb_errpush(info,name)
        goto 9999
      endif
      do i= 1, nnz
        work(i) = dzero
      enddo
      select case(iflag) 
      case(psb_dupl_ovwrt_,psb_dupl_err_) 
        do i=1, nnz 
          work(a%ia2(ip2+i-1)) = a%aspk(i) 
        enddo
      case(psb_dupl_add_) 
        do i=1, nnz 
          work(a%ia2(ip2+i-1)) = a%aspk(i) + work(a%ia2(ip2+i-1)) 
        enddo
      case default
        info = 8887
        call psb_errpush(info,name)
        goto 9999          
      end select

      do i=1, nnz
        a%aspk(i) = work(i)
      enddo


    case(psb_upd_srch_)
      ! Nothing to be done  here. 
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*) trim(name),&
           & 'Going through on regeneration with psb_upd_srch_'
    case default
      ! Wrong value
      info = 8888
      call psb_errpush(info,name)
      goto 9999

    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine coo_dsp_regen

  subroutine jad_dsp_regen(a,info)

    use psb_spmat_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    type(psb_dspmat_type), intent(inout)    :: a
    integer                              :: info

    integer :: i, ip1,ip2,nnz,iflag,ichk,nnzt
    real(psb_dpk_), allocatable  :: work(:)
    integer :: err_act
    character(len=20)  :: name
    integer            :: debug_level, debug_unit

    name='psb_spcnv'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()


    !
    !   dupl_ and upd_ fields should not be changed. 
    !
    select case(psb_sp_getifld(psb_upd_,a,info))

    case(psb_upd_perm_)

      allocate(work(size(a%aspk)+1000),stat=info)
      if (info /= 0) then
        info=2040
        call psb_errpush(info,name)
        goto 9999
      end if

      if (debug_level >= psb_debug_serial_)&
           & write(debug_unit,*) trim(name),'Regeneration with psb_upd_perm_'
      ip1   = psb_sp_getifld(psb_upd_pnt_,a,info) 
      ip2   = a%ia1(ip1+psb_ip2_)
      nnz   = a%ia1(ip1+psb_nnz_)
      iflag = a%ia1(ip1+psb_iflag_)
      ichk  = a%ia1(ip1+psb_ichk_)
      nnzt  = a%ia1(ip1+psb_nnzt_)
      if (debug_level >= psb_debug_serial_)&
           & write(debug_unit,*) trim(name),'Regeneration start: ',&
           &   a%infoa(psb_upd_),psb_upd_perm_,nnz,nnzt ,iflag,info

      if ((ichk/=nnzt+iflag).or.(nnz/=nnzt)) then               
        info = 8889
        write(debug_unit,*) trim(name),'Regeneration start error: ',&
             &   a%infoa(psb_upd_),psb_upd_perm_,nnz,nnzt ,iflag,ichk    
        call psb_errpush(info,name)
        goto 9999
      endif
      do i= 1, nnz
        work(i) = dzero
      enddo
      select case(iflag) 
      case(psb_dupl_ovwrt_,psb_dupl_err_) 
        do i=1, nnz 
          work(a%ia1(ip2+i-1)) = a%aspk(i) 
        enddo
      case(psb_dupl_add_) 
        do i=1, nnz 
          work(a%ia1(ip2+i-1)) = a%aspk(i) + work(a%ia1(ip2+i-1)) 
        enddo
      case default
        info = 8887
        call psb_errpush(info,name)
        goto 9999          
      end select

      do i=1, nnz
        a%aspk(i) = work(i)
      enddo


    case(psb_upd_srch_)
      ! Nothing to be done  here. 
      if (debug_level >= psb_debug_serial_)&
           & write(debug_unit,*) trim(name),&
           & 'Going through on regeneration with psb_upd_srch_'
    case default
      ! Wrong value
      info = 8888
      call psb_errpush(info,name)
      goto 9999

    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine jad_dsp_regen


  
  subroutine csr_zsp_regen(a,info)

    use psb_spmat_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    type(psb_zspmat_type), intent(inout)    :: a
    integer                              :: info

    integer :: i, ip1,ip2,nnz,iflag,ichk,nnzt
    complex(psb_dpk_), allocatable  :: work(:)
    integer :: err_act
    character(len=20)  :: name
    integer            :: debug_level, debug_unit
  
    name='psb_spcnv'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()


    !
    !   dupl_ and upd_ fields should not be changed. 
    !
    select case(psb_sp_getifld(psb_upd_,a,info))

    case(psb_upd_perm_)

      allocate(work(size(a%aspk)+1000),stat=info)
      if (info /= 0) then
        info=2040
        call psb_errpush(info,name)
        goto 9999
      end if

      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*)  trim(name),'Regeneration with psb_upd_perm_'
      ip1   = psb_sp_getifld(psb_upd_pnt_,a,info) 
      ip2   = a%ia2(ip1+psb_ip2_)
      nnz   = a%ia2(ip1+psb_nnz_)
      iflag = a%ia2(ip1+psb_iflag_)
      ichk  = a%ia2(ip1+psb_ichk_)
      nnzt  = a%ia2(ip1+psb_nnzt_)
      if (debug_level >= psb_debug_serial_)&
           &  write(debug_unit,*) trim(name),'Regeneration start: ',&
           &   a%infoa(psb_upd_),psb_upd_perm_,nnz,nnzt ,iflag,info

      if ((ichk/=nnzt+iflag).or.(nnz/=nnzt)) then               
        info = 8889
        write(debug_unit,*) trim(name),'Regeneration start error: ',&
             &   a%infoa(psb_upd_),psb_upd_perm_,nnz,nnzt ,iflag,ichk    
        call psb_errpush(info,name)
        goto 9999
      endif
      do i= 1, nnz
        work(i) = zzero
      enddo
      select case(iflag) 
      case(psb_dupl_ovwrt_,psb_dupl_err_) 
        do i=1, nnz 
          work(a%ia2(ip2+i-1)) = a%aspk(i) 
        enddo
      case(psb_dupl_add_) 
        do i=1, nnz 
          work(a%ia2(ip2+i-1)) = a%aspk(i) + work(a%ia2(ip2+i-1)) 
        enddo
      case default
        info = 8887
        call psb_errpush(info,name)
        goto 9999          
      end select

      do i=1, nnz
        a%aspk(i) = work(i)
      enddo


    case(psb_upd_srch_)
      ! Nothing to be done  here. 
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*) trim(name),&
           & 'Going through on regeneration with psb_upd_srch_'
    case default
      ! Wrong value
      info = 8888
      call psb_errpush(info,name)
      goto 9999

    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine csr_zsp_regen

  subroutine coo_zsp_regen(a,info)

    use psb_spmat_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    type(psb_zspmat_type), intent(inout)    :: a
    integer                              :: info

    integer :: i, ip1,ip2,nnz,iflag,ichk,nnzt
    complex(psb_dpk_), allocatable  :: work(:)
    integer :: err_act
    character(len=20)  :: name
    integer            :: debug_level, debug_unit

    name='psb_spcnv'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    !
    !   dupl_ and upd_ fields should not be changed. 
    !
    select case(psb_sp_getifld(psb_upd_,a,info))

    case(psb_upd_perm_)

      allocate(work(size(a%aspk)+1000),stat=info)
      if (info /= 0) then
        info=2040
        call psb_errpush(info,name)
        goto 9999
      end if

      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*) trim(name),'Regeneration with psb_upd_perm_'
      ip1   = psb_sp_getifld(psb_upd_pnt_,a,info) 
      ip2   = a%ia2(ip1+psb_ip2_)
      nnz   = a%ia2(ip1+psb_nnz_)
      iflag = a%ia2(ip1+psb_iflag_)
      ichk  = a%ia2(ip1+psb_ichk_)
      nnzt  = a%ia2(ip1+psb_nnzt_)
      if (debug_level >= psb_debug_serial_)&
           & write(debug_unit,*) trim(name),'Regeneration start: ',&
           &   a%infoa(psb_upd_),psb_upd_perm_,nnz,nnzt ,iflag,info

      if ((ichk/=nnzt+iflag).or.(nnz/=nnzt)) then               
        info = 8889
        write(debug_unit,*) trim(name),'Regeneration start error: ',&
             &   a%infoa(psb_upd_),psb_upd_perm_,nnz,nnzt ,iflag,ichk    
        call psb_errpush(info,name)
        goto 9999
      endif
      do i= 1, nnz
        work(i) = zzero
      enddo
      select case(iflag) 
      case(psb_dupl_ovwrt_,psb_dupl_err_) 
        do i=1, nnz 
          work(a%ia2(ip2+i-1)) = a%aspk(i) 
        enddo
      case(psb_dupl_add_) 
        do i=1, nnz 
          work(a%ia2(ip2+i-1)) = a%aspk(i) + work(a%ia2(ip2+i-1)) 
        enddo
      case default
        info = 8887
        call psb_errpush(info,name)
        goto 9999          
      end select

      do i=1, nnz
        a%aspk(i) = work(i)
      enddo


    case(psb_upd_srch_)
      ! Nothing to be done  here. 
      if (debug_level >= psb_debug_serial_)&
           &  write(debug_unit,*) trim(name),&
           & 'Going through on regeneration with psb_upd_srch_'
    case default
      ! Wrong value
      info = 8888
      call psb_errpush(info,name)
      goto 9999

    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine coo_zsp_regen

  subroutine jad_zsp_regen(a,info)

    use psb_spmat_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    type(psb_zspmat_type), intent(inout)    :: a
    integer                              :: info

    integer :: i, ip1,ip2,nnz,iflag,ichk,nnzt
    complex(psb_dpk_), allocatable  :: work(:)
    integer :: err_act
    character(len=20)  :: name
    integer            :: debug_level, debug_unit

    name='psb_spcnv'
    info = 0
    call psb_erractionsave(err_act)
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()


    !
    !   dupl_ and upd_ fields should not be changed. 
    !
    select case(psb_sp_getifld(psb_upd_,a,info))

    case(psb_upd_perm_)

      allocate(work(size(a%aspk)+1000),stat=info)
      if (info /= 0) then
        info=2040
        call psb_errpush(info,name)
        goto 9999
      end if

      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*) trim(name),'Regeneration with psb_upd_perm_'
      ip1   = psb_sp_getifld(psb_upd_pnt_,a,info) 
      ip2   = a%ia1(ip1+psb_ip2_)
      nnz   = a%ia1(ip1+psb_nnz_)
      iflag = a%ia1(ip1+psb_iflag_)
      ichk  = a%ia1(ip1+psb_ichk_)
      nnzt  = a%ia1(ip1+psb_nnzt_)
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*) trim(name),'Regeneration start: ',&
           & a%infoa(psb_upd_),psb_upd_perm_,nnz,nnzt ,iflag,info

      if ((ichk/=nnzt+iflag).or.(nnz/=nnzt)) then               
        info = 8889
        write(debug_unit,*) trim(name),'Regeneration start error: ',&
             & a%infoa(psb_upd_),psb_upd_perm_,nnz,nnzt ,iflag,ichk    
        call psb_errpush(info,name)
        goto 9999
      endif
      do i= 1, nnz
        work(i) = zzero
      enddo
      select case(iflag) 
      case(psb_dupl_ovwrt_,psb_dupl_err_) 
        do i=1, nnz 
          work(a%ia1(ip2+i-1)) = a%aspk(i) 
        enddo
      case(psb_dupl_add_) 
        do i=1, nnz 
          work(a%ia1(ip2+i-1)) = a%aspk(i) + work(a%ia1(ip2+i-1)) 
        enddo
      case default
        info = 8887
        call psb_errpush(info,name)
        goto 9999          
      end select

      do i=1, nnz
        a%aspk(i) = work(i)
      enddo


    case(psb_upd_srch_)
      ! Nothing to be done  here. 
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*) trim(name),&
           & 'Going through on regeneration with psb_upd_srch_'
    case default
      ! Wrong value
      info = 8888
      call psb_errpush(info,name)
      goto 9999

    end select

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act == psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine jad_zsp_regen

  

end module psb_regen_mod
