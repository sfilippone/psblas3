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
subroutine psb_dcoins(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl,rebuild)

  use psb_spmat_type
  use psb_const_mod
  use psb_realloc_mod
  use psb_string_mod
  use psb_error_mod
  use psb_serial_mod, only : psb_csdp
  use psb_update_mod
  implicit none 

  integer, intent(in)                  :: nz, imin,imax,jmin,jmax
  integer, intent(in)                  :: ia(:),ja(:)
  real(kind(1.d0)), intent(in)         :: val(:)
  type(psb_dspmat_type), intent(inout) :: a
  integer, intent(out)                 :: info
  integer, intent(in), optional        :: gtl(:)
  logical, intent(in), optional        :: rebuild

  character(len=5)     :: ufida
  integer              :: i,j,ir,ic,nr,nc, ng, nza, isza,spstate, nnz,&
       & ip1, nzl, err_act, int_err(5), iupd, irst
  logical, parameter   :: debug=.false.
  logical              :: rebuild_
  character(len=20)    :: name, ch_err
  type(psb_dspmat_type) :: tmp

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

  if (present(rebuild)) then 
    rebuild_ = rebuild
  else
    rebuild_ = .false.
  end if

  call touppers(a%fida,ufida)
  spstate = psb_sp_getifld(psb_state_,a,info) 


  if (present(gtl)) then 

    ng = size(gtl)
    select case(spstate) 
    case(psb_spmat_bld_) 
      if ((ufida /= 'COO').and.(ufida/='COI')) then 
        info = 134
        ch_err(1:3)=ufida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      call psb_sp_info(psb_nztotreq_,a,nza,info)
      call psb_sp_info(psb_nzsizereq_,a,isza,info)
      if(info /= izero) then
        info=4010
        ch_err='psb_spinfo'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      endif

      if ((nza+nz)>isza) then 
        call psb_sp_reall(a,max(nza+nz,int(1.5*isza)),info)
        if(info /= izero) then
          info=4010
          ch_err='psb_sp_reall'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        endif
      endif

      call psb_inner_ins(nz,ia,ja,val,nza,a%ia1,a%ia2,a%aspk,&
           & min(size(a%ia1),size(a%ia2),size(a%aspk)),&
           & imin,imax,jmin,jmax,info,gtl,ng)

      if(info /= izero) then
        
        ch_err='psb_inner_ins'
        call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
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


        call psb_inner_upd(nz,ia,ja,val,nza,a%aspk,size(a%aspk),&
             & imin,imax,jmin,jmax,nzl,info,gtl,ng)

        if (info /= izero) then
          ch_err='psb_inner_upd'
          call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
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

        call  psb_srch_upd(nz,ia,ja,val,nza,a,&
             & imin,imax,jmin,jmax,nzl,info,gtl,ng)

        if (info > 0) then 
          if (rebuild_) then 
            if (debug)  write(0,*)&
                 &  'COINS: Going through rebuild_ fingers crossed!'
            irst = info
            call psb_nullify_sp(tmp)
            tmp%fida='COO'
            call psb_csdp(a,tmp,info)
            if(info /= izero) then
              info=4010
              ch_err='psb_csdp'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            endif
            call psb_sp_setifld(psb_spmat_bld_,psb_state_,tmp,info)
            if (debug) then
              write(0,*) 'COINS  Rebuild: size',tmp%infoa(psb_nnz_) ,irst            
            endif
            call psb_sp_transfer(tmp,a,info)
            if(info /= izero) then
              info=4010
              ch_err='psb_sp_transfer'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            endif
            call psb_sp_info(psb_nztotreq_,a,nza,info)
            call psb_sp_info(psb_nzsizereq_,a,isza,info)
            if(info /= izero) then
              info=4010
              ch_err='psb_spinfo'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            endif

            if (debug)  write(0,*)&
                 &  'COINS: Reinserting',a%fida,nza,isza,irst,nz
            if ((nza+nz)>isza) then 
              call psb_sp_reall(a,max(nza+nz,int(1.5*isza)),info)
              if(info /= izero) then
                info=4010
                ch_err='psb_sp_reall'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              endif
            endif
            if (irst <= nz) then 
              call psb_inner_ins((nz-irst+1),ia(irst:nz),ja(irst:nz),val(irst:nz),&
                   & nza,a%ia1,a%ia2,a%aspk,&
                   & min(size(a%ia1),size(a%ia2),size(a%aspk)),&
                   &imin,imax,jmin,jmax,info,gtl,ng)
              if (info /= izero) then
                ch_err='psb_inner_ins'
                call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
              endif
        
              call psb_sp_setifld(nza,psb_del_bnd_,a,info)
              call psb_sp_setifld(nza,psb_nnz_,a,info)
            end if

          else
            info = 2231
            call psb_errpush(info,name)
            goto 9999 
          end if
        else if (info < 0) then
          info = 2230
          call psb_errpush(info,name)
          goto 9999 

        end if

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





  else

    ng = -1 

    select case(spstate) 
    case(psb_spmat_bld_) 
      if ((ufida /= 'COO').and.(ufida/='COI')) then 
        info = 134
        ch_err(1:3)=ufida(1:3)
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      call psb_sp_info(psb_nztotreq_,a,nza,info)
      call psb_sp_info(psb_nzsizereq_,a,isza,info)
      if(info /= izero) then
        info=4010
        ch_err='psb_spinfo'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      endif

      if ((nza+nz)>isza) then 
        call psb_sp_reall(a,max(nza+nz,int(1.5*isza)),info)
        if(info /= izero) then
          info=4010
          ch_err='psb_sp_reall'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        endif
      endif

      call psb_inner_ins(nz,ia,ja,val,nza,a%ia1,a%ia2,a%aspk,&
           & min(size(a%ia1),size(a%ia2),size(a%aspk)),&
           & imin,imax,jmin,jmax,info)

      if(info /= izero) then
        ch_err='psb_inner_ins'
        call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))

        goto 9999
      endif
      if (debug) then 
        if ((nza - a%infoa(psb_nnz_)) /= nz) then 
          write(0,*) 'PSB_COINS: insert discarded items '
        end if
      end if
      if ((nza - psb_sp_getifld(psb_nnz_,a,info)) /= nz) then
        call psb_sp_setifld(nza,psb_del_bnd_,a,info)
!!$        write(0,*) 'Settind del_bnd_ 2: ',nza
      endif
      call psb_sp_setifld(nza,psb_nnz_,a,info)

    case(psb_spmat_upd_)

      iupd = psb_sp_getifld(psb_upd_,a,info)
      select case (iupd)
      case (psb_upd_perm_)
        ip1 = psb_sp_getifld(psb_upd_pnt_,a,info)      
        nzl = psb_sp_getifld(psb_del_bnd_,a,info)      
        nza = a%ia2(ip1+psb_nnz_)


        call psb_inner_upd(nz,ia,ja,val,nza,a%aspk,size(a%aspk),&
             & imin,imax,jmin,jmax,nzl,info)
        if (info /= izero) then
          ch_err='psb_inner_upd'
          call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
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

        call  psb_srch_upd(nz,ia,ja,val,nza,a,&
             & imin,imax,jmin,jmax,nzl,info)

        if (info > 0) then 
          if (rebuild_) then 
            if (debug)  write(0,*)&
                 &  'COINS: Going through rebuild_ fingers crossed!'
            irst = info
            call psb_nullify_sp(tmp)
            tmp%fida='COO'
            call psb_csdp(a,tmp,info)
            call psb_sp_setifld(psb_spmat_bld_,psb_state_,tmp,info)
            if (debug) then
              write(0,*) 'COINS  Rebuild: size',tmp%infoa(psb_nnz_) ,irst            
            endif
            call psb_sp_transfer(tmp,a,info)
            call psb_sp_info(psb_nztotreq_,a,nza,info)
            call psb_sp_info(psb_nzsizereq_,a,isza,info)
            if(info /= izero) then
              info=4010
              ch_err='psb_spinfo'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            endif

            if (debug)  write(0,*)&
                 &  'COINS: Reinserting',a%fida,nza,isza
            if ((nza+nz)>isza) then 
              call psb_sp_reall(a,max(nza+nz,int(1.5*isza)),info)
              if(info /= izero) then
                info=4010
                ch_err='psb_sp_reall'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              endif
            endif
            if (irst <= nz) then 
              call psb_inner_ins((nz-irst+1),ia(irst:nz),ja(irst:nz),val(irst:nz),&
                   & nza,a%ia1,a%ia2,a%aspk,&
                   & min(size(a%ia1),size(a%ia2),size(a%aspk)),&
                   & imin,imax,jmin,jmax,info)
              if (info /= izero) then 
                ch_err='psb_inner_ins'
                call psb_errpush(4013,name,a_err=ch_err,i_err=(/info,0,0,0,0/))
              endif
        
              call psb_sp_setifld(nza,psb_del_bnd_,a,info)
              call psb_sp_setifld(nza,psb_nnz_,a,info)
            end if

          else
            info = 2231
            call psb_errpush(info,name)
            goto 9999 
          end if
        else if (info < 0) then
          info = 2230
          call psb_errpush(info,name)
          goto 9999 

        end if

      case default  

        info = 2233
        call psb_errpush(info,name)
        goto 9999 
      end select


    case default
      info = 2232
      call psb_errpush(info,name)
      goto 9999
    end select
  endif


  return

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

contains

  subroutine psb_inner_upd(nz,ia,ja,val,nza,aspk,maxsz,&
       & imin,imax,jmin,jmax,nzl,info,gtl,ng)
    implicit none 

    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl,maxsz
    integer, intent(in) :: ia(*),ja(*)
    integer, intent(inout) :: nza
    real(kind(1.d0)), intent(in) :: val(*)
    real(kind(1.d0)), intent(inout) :: aspk(*)
    integer, intent(out) :: info
    integer, intent(in), optional  :: ng,gtl(*)
    integer  :: i,ir,ic

    if (present(gtl)) then 
      if (.not.present(ng)) then 
        info = -1
        return
      endif
      if ((nza > nzl)) then 
        do i=1, nz 
          nza = nza + 1 
          if (nza>maxsz) then 
            write(0,*) 'Out of bounds in INNER_UPD ',nza,maxsz
            info = -71
            return
          endif
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
              if (nza>maxsz) then 
                info = -72
                return
              endif
              aspk(nza) = val(i)
            end if
          end if
        end do
      end if
    else
      if ((nza >= nzl)) then 
        do i=1, nz 
          nza = nza + 1 
          if (nza>maxsz) then 
            info = -73
            return
          endif
          aspk(nza) = val(i)
        end do
      else
        do i=1, nz 
          ir = ia(i)
          ic = ja(i) 
          if ((ir >=imin).and.(ir<=imax).and.(ic>=jmin).and.(ic<=jmax)) then 
            nza = nza + 1 
            if (nza>maxsz) then 
              info = -74
              return
            endif
            aspk(nza) = val(i)
          end if
        end do
      end if
    end if
  end subroutine psb_inner_upd

  subroutine psb_inner_ins(nz,ia,ja,val,nza,ia1,ia2,aspk,maxsz,&
       & imin,imax,jmin,jmax,info,gtl,ng)
    implicit none 

    integer, intent(in) :: nz, imin,imax,jmin,jmax,maxsz
    integer, intent(in) :: ia(*),ja(*)
    integer, intent(inout) :: nza,ia1(*),ia2(*)
    real(kind(1.d0)), intent(in) :: val(*)
    real(kind(1.d0)), intent(inout) :: aspk(*)
    integer, intent(out) :: info
    integer, intent(in), optional  :: ng,gtl(*)
    integer :: i,ir,ic

    info = 0
    if (present(gtl)) then 
      if (.not.present(ng)) then 
        info = -1
        return
      endif
      do i=1, nz 
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=1).and.(ir<=ng).and.(ic>=1).and.(ic<=ng)) then 
          ir = gtl(ir)
          ic = gtl(ic) 
          if ((ir >=imin).and.(ir<=imax).and.(ic>=jmin).and.(ic<=jmax)) then 
            nza = nza + 1 
            if (nza > maxsz) then 
              info = -91
              return
            endif
            ia1(nza) = ir
            ia2(nza) = ic
            aspk(nza) = val(i)
          end if
        end if
      end do
    else

      do i=1, nz 
        ir = ia(i)
        ic = ja(i) 
        if ((ir >=imin).and.(ir<=imax).and.(ic>=jmin).and.(ic<=jmax)) then 
          nza = nza + 1 
          if (nza > maxsz) then 
            info = -92
            return
          endif
          ia1(nza) = ir
          ia2(nza) = ic
          aspk(nza) = val(i)
        end if
      end do
    end if

  end subroutine psb_inner_ins



end subroutine psb_dcoins

