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
! File:  psb_zgather.f90
!
! Subroutine: psb_zgatherm
!   This subroutine gathers pieces of a distributed dense matrix into a local one.
!
! Parameters:
!   globx     -  real,dimension(:,:).          The local matrix into which gather 
!                                                  the distributed pieces.
!   locx      -  real,dimension(:,:).          The local piece of the distributed 
!                                                  matrix to be gathered.
!   desc_a    -  type(<psb_desc_type>).        The communication descriptor.
!   info      -  integer.                      Eventually returns an error code.
!   iroot     -  integer.                      The process that has to own the 
!                                                  global matrix. If -1 all
!                                              the processes will have a copy.
!   iiglobx   -  integer(optional).            The starting row of the global matrix. 
!   ijglobx   -  integer(optional).            The starting column of the global matrix. 
!   iilocx    -  integer(optional).            The starting row of the local piece 
!                                                  of matrix. 
!   ijlocx    -  integer(optional).            The starting column of the local piece 
!                                                  of matrix.
!   ik        -  integer(optional).            The number of columns to gather. 
!
subroutine  psb_zgatherm(globx, locx, desc_a, info, iroot,&
     & iiglobx, ijglobx, iilocx,ijlocx,ik)
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  complex(kind(1.d0)), intent(in)    :: locx(:,:)
  complex(kind(1.d0)), intent(out)   :: globx(:,:)
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(out)            :: info
  integer, intent(in), optional   :: iroot, iiglobx, ijglobx, iilocx, ijlocx, ik


  ! locals
  integer                  :: int_err(5), ictxt, np, me, &
       & err_act, n, root, iiroot, ilocx, iglobx, jlocx,&
       & jglobx, lda_locx, lda_globx, m, lock, globk, maxk, k, jlx, ilx, i, j, idx

  character(len=20)        :: name, ch_err

  name='psb_zgatherm'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(iroot)) then
    root = iroot
    if((root.lt.-1).or.(root.gt.np)) then
      info=30
      int_err(1:2)=(/5,root/)
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
  else
    root = -1
  end if
  if (root==-1) then
    iiroot=0
  else 
    iiroot = root
  endif

  if (present(iiglobx)) then
    iglobx = iiglobx
  else
    iglobx = 1
  end if

  if (present(ijglobx)) then
    jglobx = ijglobx
  else
    jglobx = 1
  end if

  if (present(iilocx)) then
    ilocx = iilocx
  else
    ilocx = 1
  end if

  if (present(ijlocx)) then
    jlocx = ijlocx
  else
    jlocx = 1
  end if

  lda_globx = size(globx,1)
  lda_locx  = size(locx, 1)

  m = psb_cd_get_global_rows(desc_a)
  n = psb_cd_get_global_cols(desc_a)

  lock=size(locx,2)-jlocx+1
  globk=size(globx,2)-jglobx+1
  maxk=min(lock,globk)

  if(present(ik)) then
    if(ik.gt.maxk) then
      k=maxk
    else
      k=ik
    end if
  else
    k = maxk
  end if

  call psb_bcast(ictxt,k,root=iiroot)

  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,size(globx,1),iglobx,jglobx,desc_a,info)
  if (info == 0) &
       & call psb_chkvect(m,n,size(locx,1),ilocx,jlocx,desc_a,info,ilx,jlx)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chk(glob)vect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((ilx.ne.1).or.(iglobx.ne.1)) then
     info=3040
     call psb_errpush(info,name)
     goto 9999
  end if
  
  globx(:,:)=0.d0

  do j=1,k
     do i=1,psb_cd_get_local_rows(desc_a)
        idx = desc_a%loc_to_glob(i)
        globx(idx,jglobx+j-1) = locx(i,jlx+j-1)
     end do
     ! adjust overlapped elements
     i=1
     do while (desc_a%ovrlap_elem(i).ne.-1)
        idx=desc_a%ovrlap_elem(i+psb_ovrlp_elem_)
        idx=desc_a%loc_to_glob(idx)
        globx(idx,jglobx+j-1) = &
             & globx(idx,jglobx+j-1)/desc_a%ovrlap_elem(i+psb_n_dom_ovr_)
        i=i+2
     end do
  end do

  call psb_sum(ictxt,globx(1:m,jglobx:jglobx+k-1),root=root)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return

end subroutine psb_zgatherm






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
! Subroutine: psb_zgatherv
!   This subroutine gathers pieces of a distributed dense vector into a local one.
!
! Parameters:
!   globx     -  real,dimension(:).            The local vector into which gather 
!                                                  the distributed pieces.
!   locx      -  real,dimension(:).            The local piece of the distributed 
!                                                  vector to be gathered.
!   desc_a    -  type(<psb_desc_type>).        The communication descriptor.
!   info      -  integer.                      Eventually returns an error code.
!   iroot     -  integer.                      The process that has to own the 
!                                                  global vector. If -1 all
!                                              the processes will have a copy.
!   iiglobx   -  integer(optional).            The starting row of the global vector. 
!   iilocx    -  integer(optional).            The starting row of the local piece 
!                                                  of vector. 
!
subroutine  psb_zgatherv(globx, locx, desc_a, info, iroot,&
     & iiglobx, iilocx)
  use psb_descriptor_type
  use psb_check_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none

  complex(kind(1.d0)), intent(in)    :: locx(:)
  complex(kind(1.d0)), intent(out)   :: globx(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(out)            :: info
  integer, intent(in), optional   :: iroot, iiglobx, iilocx


  ! locals
  integer                  :: int_err(5), ictxt, np, me, &
       & err_act, n, root, ilocx, iglobx, jlocx,&
       & jglobx, lda_locx, lda_globx, m, k, jlx, ilx, i, idx

  character(len=20)        :: name, ch_err

  name='psb_zgatherv'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(iroot)) then
    root = iroot
    if((root.lt.-1).or.(root.gt.np)) then
      info=30
      int_err(1:2)=(/5,root/)
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
  else
    root = -1
  end if

  jglobx=1
  if (present(iiglobx)) then
    iglobx = iiglobx
  else
    iglobx = 1
  end if

  jlocx=1
  if (present(iilocx)) then
    ilocx = iilocx
  else
    ilocx = 1
  end if

  lda_globx = size(globx)
  lda_locx  = size(locx)

  m = psb_cd_get_global_rows(desc_a)
  n = psb_cd_get_global_cols(desc_a)

  k = 1


  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,size(globx),iglobx,jglobx,desc_a,info)
  if (info == 0) &
       & call psb_chkvect(m,n,size(locx),ilocx,jlocx,desc_a,info,ilx,jlx)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chk(glob)vect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((ilx.ne.1).or.(iglobx.ne.1)) then
    info=3040
    call psb_errpush(info,name)
    goto 9999
  end if

  globx(:)=0.d0

  do i=1,psb_cd_get_local_rows(desc_a)
    idx = desc_a%loc_to_glob(i)
    globx(idx) = locx(i)
  end do
  ! adjust overlapped elements
  i=1
  do while (desc_a%ovrlap_elem(i).ne.-1)
    idx=desc_a%ovrlap_elem(i+psb_ovrlp_elem_)
    idx=desc_a%loc_to_glob(idx)
    globx(idx) = globx(idx)/desc_a%ovrlap_elem(i+psb_n_dom_ovr_)
    i=i+2
  end do

  call psb_sum(ictxt,globx(1:m),root=root)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_zgatherv
