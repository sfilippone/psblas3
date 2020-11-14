!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
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
! File:  psb_igather.f90
!
! Subroutine: psb_igather_vect
!   This subroutine gathers pieces of a distributed vector into a local one.
!
! Arguments:
!   globx     -  integer,dimension(:).           The local matrix into which gather 
!                                                  the distributed pieces.
!   locx      -  type(psb_i_vect_type@         The local piece of the distributed 
!                                                  vector to be gathered.
!   desc_a    -  type(psb_desc_type).        The communication descriptor.
!   info      -  integer.                      Error code.
!   iroot     -  integer.                      The process that has to own the 
!                                              global matrix. If -1 all
!                                              the processes will have a copy.
!
subroutine  psb_igather_vect(globx, locx, desc_a, info, iroot)
  use psb_base_mod, psb_protect_name => psb_igather_vect
  implicit none

  type(psb_i_vect_type), intent(inout) :: locx
  integer(psb_ipk_), intent(out), allocatable   :: globx(:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer(psb_ipk_), intent(out)              :: info
  integer(psb_ipk_), intent(in), optional     :: iroot


  ! locals
  type(psb_ctxt_type) :: ictxt
  integer(psb_mpk_) :: np, me, root, iiroot, icomm, myrank, rootrank
  integer(psb_ipk_) :: ierr(5), err_act, jlx, ilx, lda_locx, lda_globx, i
  integer(psb_lpk_) :: m, n, k, ilocx,  jlocx, idx, iglobx, jglobx
  integer(psb_ipk_), allocatable :: llocx(:)
  character(len=20)        :: name, ch_err

  name='psb_cgatherv'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

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
      ierr(1)=5; ierr(2)=root;
      call psb_errpush(info,name,i_err=ierr)
      goto 9999
    end if
  else
    root = -1
  end if

  jglobx=1
  iglobx = 1
  jlocx=1
  ilocx = 1


  m = desc_a%get_global_rows()
  n = desc_a%get_global_cols()
  lda_globx = m
  lda_locx  = locx%get_nrows()
  
  k = 1


  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,lda_globx,iglobx,jglobx,desc_a,info)
  if (info == psb_success_) &
       & call psb_chkvect(m,n,locx%get_nrows(),ilocx,jlocx,desc_a,info,ilx,jlx)
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
  
  call psb_realloc(m,globx,info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if
  
  globx(:) = izero
  llocx    = locx%get_vect()

  do i=1,desc_a%get_local_rows()
    call psb_loc_to_glob(i,idx,desc_a,info)
    globx(idx) = llocx(i)
  end do
  
  ! adjust overlapped elements
  do i=1, size(desc_a%ovrlap_elem,1)
    if (me /= desc_a%ovrlap_elem(i,3)) then 
      idx = desc_a%ovrlap_elem(i,1)
      call psb_loc_to_glob(idx,desc_a,info)
      globx(idx) = izero
    end if
  end do
  
  call psb_sum(ictxt,globx(1:m),root=root)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_igather_vect

! Subroutine: psb_igather_multivect
!   This subroutine gathers pieces of a distributed multivector into a local one.

subroutine  psb_igather_multivect(globx, locx, desc_a, info, iroot)
  use psb_base_mod, psb_protect_name => psb_igather_multivect
  implicit none

  type(psb_i_multivect_type), intent(inout) :: locx
  integer(psb_ipk_), intent(out), allocatable   :: globx(:,:)
  type(psb_desc_type), intent(in)   :: desc_a
  integer(psb_ipk_), intent(out)              :: info
  integer(psb_ipk_), intent(in), optional     :: iroot


  ! locals
  type(psb_ctxt_type) :: ictxt
  integer(psb_mpk_) :: np, me, root, iiroot, icomm, myrank, rootrank
  integer(psb_ipk_) :: ierr(5), err_act, jlx, ilx, lda_locx, lda_globx, i 
  integer(psb_lpk_) :: m, n, k, ilocx, jlocx, idx, iglobx, jglobx
  integer(psb_ipk_), allocatable :: llocx(:,:)
  character(len=20)        :: name, ch_err

  name='psb_cgatherv'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

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
      ierr(1)=5; ierr(2)=root;
      call psb_errpush(info,name,i_err=ierr)
      goto 9999
    end if
  else
    root = -1
  end if

  jglobx=1
  iglobx = 1
  jlocx=1
  ilocx = 1


  m = desc_a%get_global_rows()
  n = desc_a%get_global_cols()
  lda_globx = m
  lda_locx  = locx%get_nrows()
  k = locx%get_ncols()


  !  there should be a global check on k here!!!

  call psb_chkglobvect(m,n,lda_globx,iglobx,jglobx,desc_a,info)
  if (info == psb_success_) &
       & call psb_chkvect(m,n,locx%get_nrows(),ilocx,jlocx,desc_a,info,ilx,jlx)
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
  
  call psb_realloc(m,k,globx,info)
  if (info /= psb_success_) then 
    info=psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if
  
  globx(:,:) = izero
  llocx    = locx%get_vect()

  do i=1,desc_a%get_local_rows()
    call psb_loc_to_glob(i,idx,desc_a,info)
    globx(idx,:) = llocx(i,:)
  end do
  
  ! adjust overlapped elements
  do i=1, size(desc_a%ovrlap_elem,1)
    if (me /= desc_a%ovrlap_elem(i,3)) then 
      idx = desc_a%ovrlap_elem(i,1)
      call psb_loc_to_glob(idx,desc_a,info)
      globx(idx,:) = izero
    end if
  end do
  
  call psb_sum(ictxt,globx(1:m,1:k),root=root)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_igather_multivect
