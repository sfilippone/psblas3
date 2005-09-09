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
       & ip1, nzl, err_act, int_err(5)
  logical, parameter   :: debug=.true.
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
  spstate = a%infoa(psb_state_) 

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
        call psb_spreall(a,nza+nz,info)
        if(info.ne.izero) then
           info=4010
           ch_err='psb_spreall'
           call psb_errpush(info,name,a_err=ch_err)
           goto 9999
        endif
     endif
     call psb_inner_ins(nz,ia,ja,val,nza,a%ia1,a%ia2,a%aspk,gtl,&
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
     if ((nza - a%infoa(psb_nnz_)) /= nz) then
        a%infoa(psb_del_bnd_) = nza
     endif
     a%infoa(psb_nnz_) = nza

  case(psb_spmat_upd_)

     if (ibits(a%infoa(psb_upd_),2,1).eq.1) then 
        ip1 = a%infoa(psb_upd_pnt_)      
        nza = a%ia2(ip1+psb_nnz_)
        nzl = a%infoa(psb_del_bnd_)

        call psb_inner_upd(nz,ia,ja,val,nza,a%aspk,gtl,&
             & imin,imax,jmin,jmax,nzl,info)
        if(info.ne.izero) then
           info=4010
           ch_err='psb_inner_upd'
           call psb_errpush(info,name,a_err=ch_err)
           goto 9999
        endif
!!$      if (debug) then 
!!$        if ((nza - a%ia2(ip1+nnz_)) /= nz) then 
!!$          write(0,*) 'PSB_COINS: update discarded items '
!!$        end if
!!$      end if

        a%ia2(ip1+psb_nnz_) = nza
     else 
        info = 2231
        call psb_errpush(info,name)
        goto 9999 
     endif

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
  subroutine psb_inner_upd(nz,ia,ja,val,nza,aspk,gtl,imin,imax,jmin,jmax,nzl,info)
    implicit none 

    integer, intent(in) :: nz, imin,imax,jmin,jmax,nzl
    integer, intent(in) :: ia(*),ja(*),gtl(*)
    integer, intent(inout) :: nza
    real(kind(1.d0)), intent(in) :: val(*)
    real(kind(1.d0)), intent(inout) :: aspk(*)
    integer, intent(out) :: info
    integer  :: i,ir,ic
    
    info = 0

    if (nza >= nzl) then 
      do i=1, nz 
        nza = nza + 1 
        a%aspk(nza) = val(i)
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
            a%aspk(nza) = val(i)
          end if
        end if
      end do
    end if
    
  end subroutine psb_inner_upd

  subroutine psb_inner_ins(nz,ia,ja,val,nza,ia1,ia2,aspk,gtl,&
       & imin,imax,jmin,jmax,info)
    implicit none 

    integer, intent(in) :: nz, imin,imax,jmin,jmax
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
          a%ia1(nza) = ir
          a%ia2(nza) = ic
          a%aspk(nza) = val(i)
        end if
      end if
    end do

  end subroutine psb_inner_ins
end subroutine psb_dcoins

