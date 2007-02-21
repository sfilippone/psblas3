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
module psb_serial_mod
  use psb_spmat_type
  use psb_string_mod

  interface psb_csdp
     subroutine psb_dcsdp(a, b,info,ifc,check,trans,unitd,upd,dupl)
       use psb_spmat_type
       type(psb_dspmat_type), intent(in)   :: a
       type(psb_dspmat_type), intent(inout)  :: b
       integer, intent(out)        :: info
       integer, intent(in), optional :: ifc,upd,dupl
       character, intent(in), optional :: check,trans,unitd
     end subroutine psb_dcsdp
     subroutine psb_zcsdp(a, b,info,ifc,check,trans,unitd,upd,dupl)
       use psb_spmat_type
       type(psb_zspmat_type), intent(in)   :: a
       type(psb_zspmat_type), intent(inout)  :: b
       integer, intent(out)        :: info
       integer, intent(in), optional :: ifc,upd,dupl
       character, intent(in), optional :: check,trans,unitd
     end subroutine psb_zcsdp
  end interface

  interface psb_csrws
     subroutine psb_dcsrws(rw,a,info,trans)
       use psb_spmat_type
       type(psb_dspmat_type) :: a
       real(kind(1.d0)), allocatable   :: rw(:) 
       integer :: info
       character, optional :: trans
     end subroutine psb_dcsrws
     subroutine psb_zcsrws(rw,a,info,trans)
       use psb_spmat_type
       type(psb_zspmat_type) :: a
       complex(kind(1.d0)), allocatable :: rw(:) 
       integer :: info
       character, optional :: trans
     end subroutine psb_zcsrws
  end interface



  interface psb_cssm
     subroutine psb_dcssm(alpha,t,b,beta,c,info,trans,unitd,d)
       use psb_spmat_type
       type(psb_dspmat_type) :: t
       real(kind(1.d0)) :: alpha, beta, b(:,:), c(:,:)
       integer :: info
       character, optional :: trans, unitd
       real(kind(1.d0)), optional, target :: d(:)
     end subroutine psb_dcssm
     subroutine psb_dcssv(alpha,t,b,beta,c,info,trans,unitd,d)
       use psb_spmat_type
       type(psb_dspmat_type) :: t
       real(kind(1.d0)) :: alpha, beta, b(:), c(:)
       integer :: info
       character, optional :: trans, unitd
       real(kind(1.d0)), optional, target :: d(:)
     end subroutine psb_dcssv
     subroutine psb_zcssm(alpha,t,b,beta,c,info,trans,unitd,d)
       use psb_spmat_type
       type(psb_zspmat_type) :: t
       complex(kind(1.d0)) :: alpha, beta, b(:,:), c(:,:)
       integer :: info
       character, optional :: trans, unitd
       complex(kind(1.d0)), optional, target :: d(:)
     end subroutine psb_zcssm
     subroutine psb_zcssv(alpha,t,b,beta,c,info,trans,unitd,d)
       use psb_spmat_type
       type(psb_zspmat_type) :: t
       complex(kind(1.d0)) :: alpha, beta, b(:), c(:)
       integer :: info
       character, optional :: trans, unitd
       complex(kind(1.d0)), optional, target :: d(:)
     end subroutine psb_zcssv
  end interface

  interface psb_csmm
     subroutine psb_dcsmv(alpha,a,b,beta,c,info,trans)
       use psb_spmat_type
       type(psb_dspmat_type) :: a
       real(kind(1.d0)) :: alpha, beta, b(:), c(:)
       integer :: info
       character, optional :: trans
     end subroutine psb_dcsmv
     subroutine psb_dcsmm(alpha,a,b,beta,c,info,trans)
       use psb_spmat_type
       type(psb_dspmat_type) :: a
       real(kind(1.d0)) :: alpha, beta, b(:,:), c(:,:)
       integer :: info
       character, optional :: trans
     end subroutine psb_dcsmm
     subroutine psb_zcsmv(alpha,a,b,beta,c,info,trans)
       use psb_spmat_type
       type(psb_zspmat_type) :: a
       complex(kind(1.d0)) :: alpha, beta, b(:), c(:)
       integer :: info
       character, optional :: trans
     end subroutine psb_zcsmv
     subroutine psb_zcsmm(alpha,a,b,beta,c,info,trans)
       use psb_spmat_type
       type(psb_zspmat_type) :: a
       complex(kind(1.d0)) :: alpha, beta, b(:,:), c(:,:)
       integer :: info
       character, optional :: trans
     end subroutine psb_zcsmm
  end interface

  interface psb_fixcoo
     subroutine psb_dfixcoo(a,info,idir)
       use psb_spmat_type
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
       integer, intent(in), optional :: idir
     end subroutine psb_dfixcoo
     subroutine psb_zfixcoo(a,info,idir)
       use psb_spmat_type
       type(psb_zspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
       integer, intent(in), optional :: idir
     end subroutine psb_zfixcoo
  end interface

  interface psb_ipcoo2csr
     subroutine psb_dipcoo2csr(a,info,rwshr)
       use psb_spmat_type
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
       logical, optional :: rwshr
     end subroutine psb_dipcoo2csr
     subroutine psb_zipcoo2csr(a,info,rwshr)
       use psb_spmat_type
       type(psb_zspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
       logical, optional :: rwshr
     end subroutine psb_zipcoo2csr
  end interface

  interface psb_ipcoo2csc
     subroutine psb_dipcoo2csc(a,info,clshr)
       use psb_spmat_type
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
       logical, optional :: clshr
     end subroutine psb_dipcoo2csc
     subroutine psb_zipcoo2csc(a,info,clshr)
       use psb_spmat_type
       type(psb_zspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
       logical, optional :: clshr
     end subroutine psb_zipcoo2csc
  end interface

  interface psb_ipcsr2coo
     subroutine psb_dipcsr2coo(a,info)
       use psb_spmat_type
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
     end subroutine psb_dipcsr2coo
     subroutine psb_zipcsr2coo(a,info)
       use psb_spmat_type
       type(psb_zspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
     end subroutine psb_zipcsr2coo
  end interface

  interface psb_csprt
     subroutine psb_dcsprt(iout,a,iv,irs,ics,head,ivr,ivc)
       use psb_spmat_type
       integer, intent(in)       :: iout
       type(psb_dspmat_type), intent(in) :: a
       integer, intent(in), optional :: iv(:)
       integer, intent(in), optional :: irs,ics
       character(len=*), optional    :: head
       integer, intent(in), optional :: ivr(:),ivc(:)
     end subroutine psb_dcsprt
     subroutine psb_zcsprt(iout,a,iv,irs,ics,head,ivr,ivc)
       use psb_spmat_type
       integer, intent(in)       :: iout
       type(psb_zspmat_type), intent(in) :: a
       integer, intent(in), optional :: iv(:)
       integer, intent(in), optional :: irs,ics
       character(len=*), optional    :: head
       integer, intent(in), optional :: ivr(:),ivc(:)
     end subroutine psb_zcsprt
  end interface

  interface psb_neigh
     subroutine psb_dneigh(a,idx,neigh,n,info,lev)
       use psb_spmat_type
       type(psb_dspmat_type), intent(in) :: a   
       integer, intent(in)       :: idx 
       integer, intent(out)      :: n   
       integer, allocatable          :: neigh(:)
       integer, intent(out)  :: info
       integer, optional, intent(in) :: lev 
     end subroutine psb_dneigh
     subroutine psb_zneigh(a,idx,neigh,n,info,lev)
       use psb_spmat_type
       type(psb_zspmat_type), intent(in) :: a   
       integer, intent(in)       :: idx 
       integer, intent(out)      :: n   
       integer, allocatable      :: neigh(:)
       integer, intent(out)  :: info
       integer, optional, intent(in) :: lev 
     end subroutine psb_zneigh
  end interface

  interface psb_coins
     subroutine psb_dcoins(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl,rebuild)
       use psb_spmat_type
       integer, intent(in) :: nz, imin,imax,jmin,jmax
       integer, intent(in) :: ia(:),ja(:)
       real(kind(1.d0)), intent(in) :: val(:)
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out) :: info
       integer, intent(in), optional :: gtl(:)
       logical, optional, intent(in) :: rebuild
     end subroutine psb_dcoins
     subroutine psb_zcoins(nz,ia,ja,val,a,imin,imax,jmin,jmax,info,gtl,rebuild)
       use psb_spmat_type
       integer, intent(in) :: nz, imin,imax,jmin,jmax
       integer, intent(in) :: ia(:),ja(:)
       complex(kind(1.d0)), intent(in) :: val(:)
       type(psb_zspmat_type), intent(inout) :: a
       integer, intent(out) :: info
       integer, intent(in), optional :: gtl(:)
       logical, optional, intent(in) :: rebuild
     end subroutine psb_zcoins
  end interface


  interface psb_symbmm
     subroutine psb_dsymbmm(a,b,c,info)
       use psb_spmat_type
       type(psb_dspmat_type) :: a,b,c
       integer               :: info
     end subroutine psb_dsymbmm
     subroutine psb_zsymbmm(a,b,c,info)
       use psb_spmat_type
       type(psb_zspmat_type) :: a,b,c
       integer               :: info
     end subroutine psb_zsymbmm
  end interface

  interface psb_numbmm
     subroutine psb_dnumbmm(a,b,c)
       use psb_spmat_type
       type(psb_dspmat_type) :: a,b,c
     end subroutine psb_dnumbmm
     subroutine psb_znumbmm(a,b,c)
       use psb_spmat_type
       type(psb_zspmat_type) :: a,b,c
     end subroutine psb_znumbmm
  end interface

  interface psb_transp
     subroutine psb_dtransp(a,b,c,fmt)
       use psb_spmat_type
       type(psb_dspmat_type) :: a,b
       integer, optional :: c
       character(len=*), optional :: fmt
     end subroutine psb_dtransp
     subroutine psb_ztransp(a,b,c,fmt)
       use psb_spmat_type
       type(psb_zspmat_type) :: a,b
       integer, optional :: c
       character(len=*), optional :: fmt
     end subroutine psb_ztransp
  end interface

  interface psb_transc
     subroutine psb_ztransc(a,b,c,fmt)
       use psb_spmat_type
       type(psb_zspmat_type) :: a,b
       integer, optional :: c
       character(len=*), optional :: fmt
     end subroutine psb_ztransc
  end interface

  interface psb_rwextd
     subroutine psb_drwextd(nr,a,info,b,rowscale)
       use psb_spmat_type
       integer, intent(in) :: nr
       type(psb_dspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
       type(psb_dspmat_type), intent(in), optional  :: b
       logical, intent(in), optional :: rowscale
     end subroutine psb_drwextd
     subroutine psb_zrwextd(nr,a,info,b,rowscale)
       use psb_spmat_type
       integer, intent(in) :: nr
       type(psb_zspmat_type), intent(inout) :: a
       integer, intent(out)         :: info
       type(psb_zspmat_type), intent(in), optional  :: b
       logical, intent(in), optional :: rowscale
     end subroutine psb_zrwextd
  end interface

  interface psb_csnmi
     real(kind(1.d0)) function psb_dcsnmi(a,info,trans)
       use psb_spmat_type
       type(psb_dspmat_type), intent(in)  :: a
       integer, intent(out)       :: info
       character, optional        :: trans
     end function psb_dcsnmi
     real(kind(1.d0)) function psb_zcsnmi(a,info,trans)
       use psb_spmat_type
       type(psb_zspmat_type), intent(in)  :: a
       integer, intent(out)       :: info
       character, optional        :: trans
     end function psb_zcsnmi
  end interface


  interface psb_sp_getdiag
     subroutine psb_dspgtdiag(a,d,info)
       use psb_spmat_type
       type(psb_dspmat_type), intent(in)     :: a
       real(kind(1.d0)), intent(inout) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_dspgtdiag
     subroutine psb_zspgtdiag(a,d,info)
       use psb_spmat_type
       type(psb_zspmat_type), intent(in)     :: a
       complex(kind(1.d0)), intent(inout) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_zspgtdiag
  end interface

  interface psb_sp_scal
     subroutine psb_dspscal(a,d,info)
       use psb_spmat_type
       type(psb_dspmat_type), intent(inout) :: a
       real(kind(1.d0)), intent(in) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_dspscal
     subroutine psb_zspscal(a,d,info)
       use psb_spmat_type
       type(psb_zspmat_type), intent(inout) :: a
       complex(kind(1.d0)), intent(in) :: d(:) 
       integer, intent(out)  :: info
     end subroutine psb_zspscal
  end interface

  interface psb_sp_getblk
     subroutine psb_dspgtblk(irw,a,b,info,append,iren,lrw)
       use psb_spmat_type
       type(psb_dspmat_type), intent(in) :: a
       integer, intent(in)       :: irw
       type(psb_dspmat_type), intent(inout)    :: b
       logical, intent(in), optional :: append
       integer, intent(in), target, optional :: iren(:)
       integer, intent(in), optional :: lrw
       integer, intent(out)  :: info
     end subroutine psb_dspgtblk
     subroutine psb_zspgtblk(irw,a,b,info,append,iren,lrw)
       use psb_spmat_type
       type(psb_zspmat_type), intent(in) :: a
       integer, intent(in)       :: irw
       type(psb_zspmat_type), intent(inout)    :: b
       logical, intent(in), optional :: append
       integer, intent(in), target, optional :: iren(:)
       integer, intent(in), optional :: lrw
       integer, intent(out)  :: info
     end subroutine psb_zspgtblk
  end interface

  interface psb_sp_getrow
     subroutine psb_dspgetrow(irw,a,nz,ia,ja,val,info,iren,lrw)
       use psb_spmat_type
       type(psb_dspmat_type), intent(in) :: a
       integer, intent(in)       :: irw
       integer, intent(out)      :: nz
       integer, intent(inout)    :: ia(:), ja(:)
       real(kind(1.d0)),  intent(inout)    :: val(:)
       integer, intent(in), target, optional :: iren(:)
       integer, intent(in), optional :: lrw
       integer, intent(out)  :: info
     end subroutine psb_dspgetrow
     subroutine psb_zspgetrow(irw,a,nz,ia,ja,val,info,iren,lrw)
       use psb_spmat_type
       type(psb_zspmat_type), intent(in) :: a
       integer, intent(in)       :: irw
       integer, intent(out)      :: nz
       integer, intent(inout)    :: ia(:), ja(:)
       complex(kind(1.d0)),  intent(inout)    :: val(:)
       integer, intent(in), target, optional :: iren(:)
       integer, intent(in), optional :: lrw
       integer, intent(out)  :: info
     end subroutine psb_zspgetrow
  end interface
  
  interface psb_msort
    module procedure imsort
  end interface
  interface psb_qsort
    module procedure iqsort
  end interface

  integer, parameter :: psb_sort_up_=1, psb_sort_down_=-1
  integer, parameter :: psb_sort_ovw_idx_=0, psb_sort_keep_idx_=1

contains

  subroutine imsort(x,ix,dir,flag)
    use psb_error_mod
    implicit none 
    integer, intent(inout)           :: x(:) 
    integer, optional, intent(in)    :: dir, flag
    integer, optional, intent(inout) :: ix(:)
    
    integer  :: dir_, flag_, n, err_act
    
    character(len=20)  :: name

    name='psb_msort'

    if (present(dir)) then 
      dir_ = dir
    else
      dir_= psb_sort_up_
    end if
    select case(dir_) 
    case( psb_sort_up_, psb_sort_down_)
      ! OK keep going
    case default
      call psb_errpush(30,name,i_err=(/3,dir_,0,0,0/))
      goto 9999
    end select
      
    n = size(x)
 
    if (present(ix)) then 
      if (size(ix) < n) then 
        call psb_errpush(35,name,i_err=(/2,size(ix),0,0,0/))
        goto 9999
      end if
      if (present(flag)) then 
        flag_ = flag
      else 
        flag_ = psb_sort_ovw_idx_
      end if
      select case(flag_) 
      case( psb_sort_ovw_idx_, psb_sort_keep_idx_)
        ! OK keep going
      case default
        call psb_errpush(30,name,i_err=(/4,flag_,0,0,0/))
        goto 9999
      end select

      call imsrx(n,x,ix,dir_,flag_)
    else
      call imsr(n,x,dir_)
    end if

9999 continue 
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
  end subroutine imsort


  subroutine iqsort(x,ix,dir,flag)
    use psb_error_mod
    implicit none 
    integer, intent(inout)           :: x(:) 
    integer, optional, intent(in)    :: dir, flag
    integer, optional, intent(inout) :: ix(:)
    
    integer  :: dir_, flag_, n, err_act
    
    character(len=20)  :: name

    name='psb_qsort'

    if (present(dir)) then 
      dir_ = dir
    else
      dir_= psb_sort_up_
    end if
    select case(dir_) 
    case( psb_sort_up_, psb_sort_down_)
      ! OK keep going
    case default
      call psb_errpush(30,name,i_err=(/3,dir_,0,0,0/))
      goto 9999
    end select
      
    n = size(x)
 
    if (present(ix)) then 
      if (size(ix) < n) then 
        call psb_errpush(35,name,i_err=(/2,size(ix),0,0,0/))
        goto 9999
      end if
      if (present(flag)) then 
        flag_ = flag
      else 
        flag_ = psb_sort_ovw_idx_
      end if
      select case(flag_) 
      case( psb_sort_ovw_idx_, psb_sort_keep_idx_)
        ! OK keep going
      case default
        call psb_errpush(30,name,i_err=(/4,flag_,0,0,0/))
        goto 9999
      end select

      call isrx(n,x,ix,dir_,flag_)
    else
      call isr(n,x,dir_)
    end if

9999 continue 
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
  end subroutine iqsort



end module psb_serial_mod

