!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File:  psb_cgather.f90
!
! Subroutine: psb_cgatherm
!   This subroutine gathers pieces of a distributed dense matrix into a local one.
!
! Arguments:
!   globx     -  cplx,dimension(:,:).          The local matrix into which gather 
!                                                  the distributed pieces.
!   locx      -  cplx,dimension(:,:).          The local piece of the distributed 
!                                                  matrix to be gathered.
!   desc_a    -  type(psb_desc_type).        The communication descriptor.
!   info      -  integer.                      Error code.
!   iroot     -  integer.                      The process that has to own the 
!                                              global matrix. If -1 all
!                                              the processes will have a copy.
!
subroutine  psb_cgatherm(globx, locx, desc_a, info, iroot)
  use psb_base_mod, psb_protect_name => psb_cgatherm
  implicit none

  complex(psb_spk_), intent(in)    :: locx(:,:)
  complex(psb_spk_), intent(out)   :: globx(:,:)
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(out)            :: info
  integer, intent(in), optional   :: iroot


  ! locals
  integer                  :: int_err(5), ictxt, np, me, &
       & err_act, n, root, iiroot, ilocx, iglobx, jlocx,&
       & jglobx, lda_locx, lda_globx, m, lock, globk, maxk, k, jlx, ilx, i, j, idx

  character(len=20)        :: name, ch_err

  name='psb_cgatherm'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  ictxt=desc_a%get_context()

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(iroot)) then
    root = iroot
    if((root < -1).or.(root > np)) then
      info=psb_err_input_value_invalid_i_
      int_err(1:2)=(/5,root/)
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
  else
    root = -1
  end if
  if (root == -1) then
    iiroot = psb_root_
  else 
    iiroot = root
  endif

  iglobx = 1
  jglobx = 1
  ilocx = 1
  jlocx = 1

  lda_globx = size(globx,1)
  lda_locx  = size(locx, 1)

  m = desc_a%get_global_rows()
  n = desc_a%get_global_cols()

  lock=size(locx,2)-jlocx+1
  globk=size(globx,2)-jglobx+1
  maxk=min(lock,globk)
  k = maxk

  call psb_bcast(ictxt,k,root=iiroot)

  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,size(globx,1),iglobx,jglobx,desc_a,info)
  if (info == psb_success_) &
       & call psb_chkvect(m,n,size(locx,1),ilocx,jlocx,desc_a,info,ilx,jlx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chk(glob)vect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((ilx /= 1).or.(iglobx /= 1)) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
    goto 9999
  end if

  globx(:,:)=0.d0

  do j=1,k
    do i=1,desc_a%get_local_rows()
      call psb_loc_to_glob(i,idx,desc_a,info)
      globx(idx,jglobx+j-1) = locx(i,jlx+j-1)
    end do
  end do

  do j=1,k
    ! adjust overlapped elements
    do i=1, size(desc_a%ovrlap_elem,1)
      if (me /= desc_a%ovrlap_elem(i,3)) then 
        idx = desc_a%ovrlap_elem(i,1)
        call psb_loc_to_glob(idx,desc_a,info)
        globx(idx,jglobx+j-1) = czero
      end if
    end do
  end do

  call psb_sum(ictxt,globx(1:m,jglobx:jglobx+k-1),root=root)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_cgatherm






!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! Subroutine: psb_cgatherv
!   This subroutine gathers pieces of a distributed dense vector into a local one.
!
! Arguments:
!   globx     -  cplx,dimension(:).            The local vector into which gather 
!                                                  the distributed pieces.
!   locx      -  cplx,dimension(:).            The local piece of the distributed 
!                                                  vector to be gathered.
!   desc_a    -  type(psb_desc_type).        The communication descriptor.
!   info      -  integer.                      Error code.
!   iroot     -  integer.                      The process that has to own the 
!                                              global matrix. If -1 all
!                                              the processes will have a copy.
!                                              default: -1
!
subroutine  psb_cgatherv(globx, locx, desc_a, info, iroot)
  use psb_base_mod, psb_protect_name => psb_cgatherv
  implicit none

  complex(psb_spk_), intent(in)    :: locx(:)
  complex(psb_spk_), intent(out)   :: globx(:)
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(out)            :: info
  integer, intent(in), optional   :: iroot


  ! locals
  integer                  :: int_err(5), ictxt, np, me, &
       & err_act, n, root, ilocx, iglobx, jlocx,&
       & jglobx, lda_locx, lda_globx, m, k, jlx, ilx, i, idx

  character(len=20)        :: name, ch_err

  name='psb_cgatherv'
  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  ictxt=desc_a%get_context()

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (present(iroot)) then
    root = iroot
    if((root < -1).or.(root > np)) then
      info=psb_err_input_value_invalid_i_
      int_err(1:2)=(/5,root/)
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if
  else
    root = -1
  end if

  jglobx=1
  iglobx = 1
  jlocx=1
  ilocx = 1

  lda_globx = size(globx)
  lda_locx  = size(locx)

  m = desc_a%get_global_rows()
  n = desc_a%get_global_cols()

  k = 1


  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,size(globx),iglobx,jglobx,desc_a,info)
  if (info == psb_success_) &
       & call psb_chkvect(m,n,size(locx),ilocx,jlocx,desc_a,info,ilx,jlx)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chk(glob)vect'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if ((ilx /= 1).or.(iglobx /= 1)) then
    info=psb_err_ix_n1_iy_n1_unsupported_
    call psb_errpush(info,name)
    goto 9999
  end if

  globx(:)=0.d0

  do i=1,desc_a%get_local_rows()
    call psb_loc_to_glob(i,idx,desc_a,info)
    globx(idx) = locx(i)
  end do
  
  ! adjust overlapped elements
  do i=1, size(desc_a%ovrlap_elem,1)
    if (me /= desc_a%ovrlap_elem(i,3)) then 
      idx = desc_a%ovrlap_elem(i,1)
      call psb_loc_to_glob(idx,desc_a,info)
      globx(idx) = czero
    end if
  end do
  
  call psb_sum(ictxt,globx(1:m),root=root)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_cgatherv
