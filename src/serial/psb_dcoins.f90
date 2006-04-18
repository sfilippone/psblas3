!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
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
! File:  psbdcoins.f90 
 ! Subroutine: 
 ! Parameters:
subroutine psb_dcoins(nz,ia,ja,val,a,gtl,imin,imax,jmin,jmax,info)

  use psb_spmat_type
  use psb_const_mod
  use psb_realloc_mod
  use psb_string_mod
  use psb_error_mod
  use psb_serial_mod, only : psb_spinfo
  implicit none 

  integer, intent(in)                  :: nz, imin,imax,jmin,jmax
  integer, intent(in)                  :: ia(:),ja(:),gtl(:)
  real(kind(1.d0)), intent(in)         :: val(:)
  type(psb_dspmat_type), intent(inout) :: a
  integer, intent(out)                 :: info

  character(len=5)     :: ufida
  integer              :: i,j,ir,ic,nr,nc, ng, nza, isza,spstate, nnz,&
       & ip1, nzl, err_act, int_err(5), iupd
  logical, parameter   :: debug=.false.
  character(len=20)    :: name, ch_err

  name='psb_dcoins'
  info  = 0
  call psb_erractionsave(err_act)

  info = 0
  if (nz <= 0) then 
    info = 10
    int_err(1)=1
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if
  if (size(ia) < nz) then 
    info = 35
    int_err(1)=2
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  if (size(ja) < nz) then 
    info = 35
    int_err(1)=3
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if
  if (size(val) < nz) then 
    info = 35
    int_err(1)=4
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if


!!$  ufida = toupper(a%fida)
  call touppers(a%fida,ufida)
  ng = size(gtl)
  spstate = psb_sp_getifld(psb_state_,a,info) 

  select case(spstate) 

  case(psb_spmat_bld_) 
    if ((ufida /= 'COO').and.(ufida/='COI')) then 
      info = 134
      ch_err(1:3)=ufida(1:3)
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    call psb_spinfo(psb_nztotreq_,a,nza,info)
    call psb_spinfo(psb_nzsizereq_,a,isza,info)
    if(info.ne.izero) then
      info=4010
      ch_err='psb_spinfo'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    endif

    if ((nza+nz)>isza) then 
      call psb_sp_reall(a,max(nza+nz,int(1.5*isza)),info)
      if(info.ne.izero) then
        info=4010
        ch_err='psb_sp_reall'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      endif
    endif
    call psb_inner_ins(nz,ia,ja,val,nza,a%ia1,a%ia2,a%aspk,gtl,ng,&
         & imin,imax,jmin,jmax,info)
    if(info.ne.izero) then
      info=4010
      ch_err='psb_inner_ins'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    endif
    if (debug) then 
      if ((nza - a%infoa(psb_nnz_)) /= nz) then 
        write(0,*) 'PSB_COINS: insert discarded items '
      end if
    end if
    if ((nza - psb_sp_getifld(psb_nnz_,a,info)) /= nz) then
      call psb_sp_setifld(nza,psb_del_bnd_,a,info)
    endif
    call psb_sp_setifld(nza,psb_nnz_,a,info)

  case(psb_spmat_upd_)

    iupd = psb_sp_getifld(psb_upd_,a,info)
    select case (iupd)
    case (psb_upd_perm_)
      ip1 = psb_sp_getifld(psb_upd_pnt_,a,info)      
      nzl = psb_sp_getifld(psb_del_bnd_,a,info)      
      nza = a%ia2(ip1+psb_nnz_)

      nza = a%ia2(ip1+psb_nnz_)
      nzl = a%infoa(psb_del_bnd_)

      call psb_inner_upd(nz,ia,ja,val,nza,a%aspk,gtl,ng,&
           & imin,imax,jmin,jmax,nzl,info)

      if(info.ne.izero) then
        info=4010
        ch_err='psb_inner_upd'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      endif
      if (debug) then 
        if ((nza - a%ia2(ip1+psb_nnz_)) /= nz) then 
          write(0,*) 'PSB_COINS: update discarded items '
        end if
      end if
      a%ia2(ip1+psb_nnz_) = nza

      if (debug) write(0,*) 'From COINS(UPD) : NZA:',nza

    case (psb_upd_dflt_, psb_upd_srch_)

      select case(toupper(a%fida))
      case ('CSR') 
!!$      write(0,*) 'Calling csr_inner_upd'
        call  csr_inner_upd(nz,ia,ja,val,nza,a,gtl,ng,&
             & imin,imax,jmin,jmax,nzl,info)
!!$      write(0,*) 'From csr_inner_upd:',info
      case ('COO') 
        call  coo_inner_upd(nz,ia,ja,val,nza,a,gtl,ng,&
             & imin,imax,jmin,jmax,nzl,info)

      case default

        info = 2230
        call psb_errpush(info,name)
        goto 9999 
      end select

    case default  

      info = 2231
      call psb_errpush(info,name)
      goto 9999 
    end select


  case default
    info = 2232
    call psb_errpush(info,name)
    goto 9999
  end select
  return

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error()
    return
  end if
  return

contains

  subroutine psb_inner_upd(nz,ia,ja,val,nza,aspk,gtl,ng,&
       & imin,imax,jmin,jmax,nzl,info)
    implicit none 

    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl,ng
    integer, intent(in) :: ia(*),ja(*),gtl(*)
    integer, intent(inout) :: nza
    real(kind(1.d0)), intent(in) :: val(*)
    real(kind(1.d0)), intent(inout) :: aspk(*)
    integer, intent(out) :: info
    integer  :: i,ir,ic

    if (nza >= nzl) then 
      do i=1, nz 
        nza = nza + 1 
        aspk(nza) = val(i)
      end do
    else
      do i=1, nz 
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          if ((ir >=imin).and.(ir<=imax).and.(ic>=jmin).and.(ic<=jmax)) then 
            nza = nza + 1 
            aspk(nza) = val(i)
          end if
        end if
      end do
    end if

  end subroutine psb_inner_upd

  subroutine psb_inner_ins(nz,ia,ja,val,nza,ia1,ia2,aspk,gtl,ng,&
       & imin,imax,jmin,jmax,info)
    implicit none 

    integer, intent(in) :: nz, imin,imax,jmin,jmax,ng
    integer, intent(in) :: ia(*),ja(*),gtl(*)
    integer, intent(inout) :: nza,ia1(*),ia2(*)
    real(kind(1.d0)), intent(in) :: val(*)
    real(kind(1.d0)), intent(inout) :: aspk(*)
    integer, intent(out) :: info

    integer :: i,ir,ic

    info = 0
    do i=1, nz 
      ir = ia(i)
      ic = ja(i) 

      if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
        ir = gtl(ir)
        ic = gtl(ic) 
        if ((ir >=imin).and.(ir<=imax).and.(ic>=jmin).and.(ic<=jmax)) then 
          nza = nza + 1 
          ia1(nza) = ir
          ia2(nza) = ic
          aspk(nza) = val(i)
        end if
      end if
    end do

  end subroutine psb_inner_ins


  subroutine csr_inner_upd(nz,ia,ja,val,nza,a,gtl,ng,&
       & imin,imax,jmin,jmax,nzl,info)
    implicit none 

    type(psb_dspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl,ng
    integer, intent(in) :: ia(*),ja(*),gtl(*)
    integer, intent(inout) :: nza
    real(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer  :: i,ir,ic,check_flag, ilr, ilc, ip, &
         & i1,i2,nc,lb,ub,m,nnz,dupl

    info = 0

    dupl = psb_sp_getifld(psb_dupl_,a,info)

    select case(dupl)
    case(psb_dupl_ovwrt_,psb_dupl_err_)
      ! Overwrite.
      ! Cannot test for error, should have been caught earlier.

      ilr = -1 
      ilc = -1 
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          i1 = a%ia2(ir)
          i2 = a%ia2(ir+1)
          nc=i2-i1


          if (.true.) then 
            call issrch(ip,ic,nc,a%ia1(i1:i2-1))    
            if (ip>0) then 
              a%aspk(i1+ip-1) = val(i)
            else
              write(0,*)'Was searching ',ic,' in: ',i1,i2,' : ',a%ia1(i1:i2-1)
              info = -1
              return
            end if

          else
!!$ 
            ip = -1
            lb = i1
            ub = i2-1
            do
              if (lb > ub) exit
              m = (lb+ub)/2
!!$              write(0,*) 'Debug: ',lb,m,ub
              if (ic == a%ia1(m))  then
                ip   = m 
                lb   = ub + 1
              else if (ic < a%ia1(m))  then
                ub = m-1
              else 
                lb = m + 1
              end if
            enddo

            if (ip>0) then 
              a%aspk(ip) = val(i)                
            else
              write(0,*)'Was searching ',ic,' in: ',i1,i2,' : ',a%ia1(i1:i2-1)
              info = -1
              return
            end if

          end if
        end if
      end do

    case(psb_dupl_add_)
      ! Add
      ilr = -1 
      ilc = -1 
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          i1 = a%ia2(ir)
          i2 = a%ia2(ir+1)
          nc = i2-i1
          call issrch(ip,ic,nc,a%ia1(i1:i2-1))
          if (ip>0) then 
            a%aspk(i1+ip-1) = a%aspk(i1+ip-1) + val(i)
          else
            info = -2 
            return
          end if
        end if
      end do

    case default
      info = -3
      write(0,*) 'Duplicate handling: ',dupl
    end select
  end subroutine csr_inner_upd

  subroutine coo_inner_upd(nz,ia,ja,val,nza,a,gtl,ng,&
       & imin,imax,jmin,jmax,nzl,info)
    implicit none 

    type(psb_dspmat_type), intent(inout) :: a
    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl,ng
    integer, intent(in) :: ia(*),ja(*),gtl(*)
    integer, intent(inout) :: nza
    real(kind(1.d0)), intent(in) :: val(*)
    integer, intent(out) :: info
    integer  :: i,ir,ic,check_flag, ilr, ilc, ip, &
         & i1,i2,nc,lb,ub,m,nnz,dupl,isrt

    info = 0

    dupl = psb_sp_getifld(psb_dupl_,a,info)

    if (psb_sp_getifld(psb_srtd_,a,info) /= psb_isrtdcoo_) then 
      info = -4
      return
    end if

    ilr = -1 
    ilc = -1 
    nnz = psb_sp_getifld(psb_nnz_,a,info)


    select case(dupl)
    case(psb_dupl_ovwrt_,psb_dupl_err_)
      ! Overwrite.
      ! Cannot test for error, should have been caught earlier.
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          if (ir /= ilr) then 
            call ibsrch(i1,ir,nnz,a%ia1)
            i2 = i1
            do 
              if (i2+1 > nnz) exit
              if (a%ia1(i2+1) /= a%ia1(i2)) exit
              i2 = i2 + 1
            end do
            do 
              if (i1-1 < 1) exit
              if (a%ia1(i1-1) /= a%ia1(i1)) exit
              i1 = i1 - 1
            end do
            ilr = ir
          end if
          nc = i2-i1+1
          call issrch(ip,ic,nc,a%ia2(i1:i2))
          if (ip>0) then 
            a%aspk(i1+ip-1) = val(i)
          else
            info = -2 
            return
          end if
        end if
      end do
    case(psb_dupl_add_)
      ! Add
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          if (ir /= ilr) then 
            call ibsrch(i1,ir,nnz,a%ia1)
            i2 = i1
            do 
              if (i2+1 > nnz) exit
              if (a%ia1(i2+1) /= a%ia1(i2)) exit
              i2 = i2 + 1
            end do
            do 
              if (i1-1 < 1) exit
              if (a%ia1(i1-1) /= a%ia1(i1)) exit
              i1 = i1 - 1
            end do
            ilr = ir
          end if
          nc = i2-i1+1
          call issrch(ip,ic,nc,a%ia2(i1:i2))
          if (ip>0) then 
            a%aspk(i1+ip-1) = a%aspk(i1+ip-1) + val(i)
          else
            info = -2 
            return
          end if
        end if
      end do

    case default
      info = -3
      write(0,*) 'Duplicate handling: ',dupl
    end select

  end subroutine coo_inner_upd

end subroutine psb_dcoins

