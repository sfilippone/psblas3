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
! File: psb_zspcnv.f90
!
! Subroutine: psb_zspcnv
!    converts sparse matrix a into b
! 
! Parameters: 
!    a        - type(<psb_zspmat_type>).          The sparse input matrix.      
!    b        - type(<psb_zspmat_type>).          The sparse output matrix.
!    desc_a   - type(<psb_desc_type>).            The communication descriptor.
!    info     - integer.                          Eventually returns an error code.
!
subroutine psb_zspcnv(a,b,desc_a,info)

  use psb_descriptor_type
  use psb_spmat_type
  use psb_realloc_mod
  use psb_serial_mod
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  interface zcsdp

    subroutine zcsdp(check,trans,m,n,unitd,d,&
         & fida,descra,a,ia1,ia2,infoa,&
         & pl,fidh,descrh,h,ih1,ih2,infoh,pr,lh,lh1,lh2,&
         & work,lwork,ierror)
      integer, intent(in)   :: lh, lwork, lh1, lh2, m, n                 
      integer, intent(out)  :: ierror                 
      character, intent(in) :: check, trans, unitd                               
      complex(kind(1.d0)), intent(in)  :: d(*), a(*)
      complex(kind(1.d0)), intent(out) :: h(*)
      complex(kind(1.d0)), intent(inout) :: work(*)
      integer, intent(in)  :: ia1(*), ia2(*), infoa(*)
      integer, intent(out) :: ih1(*), ih2(*), pl(*),pr(*), infoh(*) 
      character, intent(in) ::  fida*5, descra*11
      character, intent(out) :: fidh*5, descrh*11
    end subroutine zcsdp
  end interface


  interface zcsrp

    subroutine zcsrp(trans,m,n,fida,descra,ia1,ia2,&
         & infoa,p,work,lwork,ierror)
      integer, intent(in)  :: m, n, lwork
      integer, intent(out) :: ierror
      character, intent(in) ::       trans
      complex(kind(1.d0)), intent(inout) :: work(*)                     
      integer, intent(in)    :: p(*)
      integer, intent(inout) :: ia1(*), ia2(*), infoa(*) 
      character, intent(in)  :: fida*5, descra*11
    end subroutine zcsrp
  end interface

  interface zcsprt
    subroutine zcsprt(m,n,fida,descra,a,ia1,ia2,infoa ,iout,ierror)
      integer, intent(in)  ::  iout,m, n                 
      integer, intent(out) ::  ierror                 
      complex(kind(1.d0)), intent(in) :: a(*)
      integer, intent(in)   :: ia1(*), ia2(*), infoa(*)
      character, intent(in) :: fida*5, descra*11
    end subroutine zcsprt
  end interface

  !...parameters....
  type(psb_zspmat_type), intent(in)   :: a
  type(psb_zspmat_type), intent(out)  :: b
  type(psb_desc_type), intent(in)     :: desc_a
  integer, intent(out)                :: info
  !....locals....
  integer                       ::  int_err(5)
  complex(kind(1.d0))           ::  d(1)
  integer,allocatable           ::  i_temp(:)
  complex(kind(1.d0)),allocatable ::  work_dcsdp(:)
  integer                       ::  ia1_size,ia2_size,aspk_size,err_act&
       & ,i,err,np,me,n_col,l_dcsdp
  integer                       ::  lwork_dcsdp,dectype
  integer                       ::  ictxt,n_row
  character                     ::  check*1, trans*1, unitd*1

  logical, parameter :: debug=.false.
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=0
  name = 'psb_zspcnv'
  call psb_erractionsave(err_act)


  ictxt   = psb_cd_get_context(desc_a)
  dectype = psb_cd_get_dectype(desc_a)
  n_row   = psb_cd_get_local_rows(desc_a)
  n_col   = psb_cd_get_local_cols(desc_a)

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.psb_is_ok_dec((dectype))) then
    info = 600
    int_err(1) = dectype
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  endif

  if (debug) write (0, *) name,'   begin matrix assembly...'

  ia1_size  = size(a%ia1)
  ia2_size  = size(a%ia2)
  aspk_size = size(a%aspk)

  if (debug) write (0, *) name,'  sizes',ia1_size,ia2_size,aspk_size

  ! convert only without check
  check='N'
  trans='N'
  unitd='U'

  ! l_dcsdp is the size requested for dcsdp procedure
  l_dcsdp=(ia1_size+100)

  b%m=n_row
  b%k=n_col
  call psb_sp_all(b,ia1_size,ia2_size,aspk_size,info)
  allocate(work_dcsdp(l_dcsdp),stat=info)
  if (info /= 0) then
    info=2025
    int_err(1)=l_dcsdp
    call psb_errpush(info, name, i_err=int_err)
    goto 9999
  endif

  lwork_dcsdp=size(work_dcsdp)
  ! set infoa(1) to nnzero
  b%pl(:)  = 0
  b%pr(:)  = 0

  if (debug) write (0, *) name,'  calling dcsdp',lwork_dcsdp,&
       &size(work_dcsdp)
  ! convert aspk,ia1,ia2 in requested representation mode
  if (debug) then

  endif
  ! result is put in b
  call zcsdp(check,trans,n_row,n_col,unitd,d,a%fida,a%descra,&
       & a%aspk,a%ia1,a%ia2,a%infoa,&
       & b%pl,b%fida,b%descra,b%aspk,b%ia1,b%ia2,b%infoa,b%pr,&
       & size(b%aspk),size(b%ia1),size(b%ia2),&
       & work_dcsdp,size(work_dcsdp),info)

  if(info /= psb_no_err_) then
    info=4010
    ch_err='zcsdp'
    call psb_errpush(info, name, a_err=ch_err)
    goto 9999
  end if

  !
  !  hmmm, have to fix b%pl and b%pr according to a%pl and a%pr!!! 
  !  should work (crossed fingers :-)
  if (a%pr(1) /= 0) then 
    if (b%pr(1) /= 0) then 
      allocate(i_temp(n_col))
      do i=1,  n_col
        i_temp(i) = b%pr(a%pr(i))
      enddo
      call psb_transfer(i_temp,b%pr,info)
    else
      allocate(i_temp(n_col))
      do i=1,  n_col
        i_temp(i) = a%pr(i)
      enddo
      call psb_transfer(i_temp,b%pr,info)
    endif
  endif
  if (a%pl(1) /= 0) then 
    if (b%pr(1) /= 0) then 
      allocate(i_temp(n_row))
      do i=1,  n_row
        i_temp(i) = a%pl(b%pl(i))
      enddo
      call psb_transfer(i_temp,b%pl,info)
    else
      allocate(i_temp(n_row))
      do i=1,  n_row
        i_temp(i) = a%pl(i)
      enddo
      call psb_transfer(i_temp,b%pl,info)
    endif
  endif


  if (debug) write (0, *) me,name,'  from zcsdp ',&
       &b%fida,' pl ', b%pl(:),'pr',b%pr(:)

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_zspcnv
