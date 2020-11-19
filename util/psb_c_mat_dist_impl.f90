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
subroutine psb_cmatdist(a_glob, a, ctxt, desc_a,&
     & info, parts, vg, vsz, inroot,fmt,mold)
  !
  ! an utility subroutine to distribute a matrix among processors
  ! according to a user defined data distribution, using
  ! sparse matrix subroutines.
  !
  !  type(psb_cspmat)                       :: a_glob
  !     on entry: this contains the global sparse matrix as follows:
  !
  !  type(psb_cspmat_type)                            :: a
  !     on exit : this will contain the local sparse matrix.
  !
  !       interface parts
  !         !   .....user passed subroutine.....
  !         subroutine parts(global_indx,n,np,pv,nv)
  !           implicit none
  !           integer(psb_ipk_), intent(in)  :: global_indx, n, np
  !           integer(psb_ipk_), intent(out) :: nv
  !           integer(psb_ipk_), intent(out) :: pv(*)
  !
  !       end subroutine parts
  !       end interface
  !     on entry:  subroutine providing user defined data distribution.
  !        for each global_indx the subroutine should return
  !        the list  pv of all processes owning the row with
  !        that index; the list will contain nv entries.
  !        usually nv=1; if nv >1 then we have an overlap in the data
  !        distribution.
  !
  !  integer(psb_ipk_) :: ctxt
  !     on entry: the PSBLAS parallel environment context.
  !
  !  type (desc_type)                  :: desc_a
  !     on exit : the updated array descriptor
  !
  !  integer(psb_ipk_), optional    :: inroot
  !     on entry: specifies processor holding a_glob. default: 0
  !     on exit : unchanged.
  !
  use psb_base_mod
  use psb_mat_mod
  implicit none

  ! parameters
  type(psb_cspmat_type)      :: a_glob
  type(psb_ctxt_type) :: ctxt
  type(psb_cspmat_type)      :: a
  type(psb_desc_type)        :: desc_a
  integer(psb_ipk_), intent(out)       :: info
  integer(psb_ipk_), optional       :: inroot
  character(len=*), optional :: fmt
  class(psb_c_base_sparse_mat), optional :: mold
  procedure(psb_parts), optional  :: parts
  integer(psb_ipk_), optional     :: vg(:)
  integer(psb_ipk_), optional     :: vsz(:)

  ! local variables
  logical           :: use_parts, use_vg, use_vsz
  integer(psb_ipk_) :: np, iam, np_sharing
  integer(psb_ipk_) :: k_count, root, liwork,  nnzero, nrhs,&
       & i, ll, nz, isize, iproc, nnr, err, err_act
  integer(psb_lpk_) :: i_count, j_count, nrow, ncol, ig, lastigp
  integer(psb_ipk_), allocatable  :: iwork(:), iwrk2(:)
  integer(psb_lpk_), allocatable  :: irow(:),icol(:)
  complex(psb_spk_), allocatable    :: val(:)
  integer(psb_ipk_), parameter    :: nb=30
  real(psb_dpk_)              :: t0, t1, t2, t3, t4, t5
  character(len=20)           :: name, ch_err

  info = psb_success_
  err  = 0
  name = 'psb_c_mat_dist'
  call psb_erractionsave(err_act)

  ! executable statements    
  if (present(inroot)) then
    root = inroot
  else
    root = psb_root_
  end if
  call psb_info(ctxt, iam, np)     

  use_parts = present(parts)
  use_vg    = present(vg)
  use_vsz   = present(vsz)
  if (count((/ use_parts, use_vg, use_vsz /)) /= 1) then 
    info=psb_err_no_optional_arg_
    call psb_errpush(info,name,a_err=" vg, vsz, parts")
    goto 9999 
  endif
  
  if (iam == root) then
    nrow = a_glob%get_nrows()
    ncol = a_glob%get_ncols()
    if (nrow /= ncol) then
      write(psb_err_unit,*) 'a rectangular matrix ? ',nrow,ncol
      info=-1
      call psb_errpush(info,name)
      goto 9999
    endif
    nnzero = a_glob%get_nzeros()
    nrhs   = 1
    if (use_vsz) then 
      if (sum(vsz(1:np)) /= nrow) then
        write(0,*) 'Input data mismatch :',nrow,sum(vsz(1:np))
      end if
    end if
    
  endif
  ! broadcast informations to other processors
  call psb_bcast(ctxt,nrow, root)
  call psb_bcast(ctxt,ncol, root)
  call psb_bcast(ctxt,nnzero, root)
  call psb_bcast(ctxt,nrhs, root)
  liwork = max(np, nrow + ncol)
  allocate(iwork(liwork), iwrk2(np),stat = info)
  if (info /= psb_success_) then
    info=psb_err_alloc_request_
    call psb_errpush(info,name,i_err=(/liwork/),a_err='integer')
    goto 9999
  endif
  if (iam == root) then
    write (*, fmt = *) 'start matdist',root, size(iwork),&
         &nrow, ncol, nnzero,nrhs
  endif
  if (use_parts) then 
    call psb_cdall(ctxt,desc_a,info,mg=nrow,parts=parts)
  else if (use_vg) then 
    call psb_cdall(ctxt,desc_a,info,vg=vg)
  else if (use_vsz) then
    call psb_cdall(ctxt,desc_a,info,nl=vsz(iam+1))
  else
    info = -1
  end if
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_cdall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  call psb_spall(a,desc_a,info,nnz=((nnzero+np-1)/np))
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_spall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  isize = 3*nb*max(((nnzero+nrow)/nrow),nb)
  allocate(val(isize),irow(isize),icol(isize),stat=info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='Allocate'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  i_count   = 1
  if (use_vsz) then
    iproc = 0
    lastigp = vsz(iproc+1)
  end if
  do while (i_count <= nrow)

    if (use_parts) then 
      call parts(i_count,nrow,np,iwork, np_sharing)
      !
      ! np_sharing allows for overlap in the data distribution.
      ! If an index is overlapped, then we have to send its row
      ! to multiple processes. NOTE: we are assuming the output 
      ! from PARTS is coherent, otherwise a deadlock is the most
      ! likely outcome.
      ! 
      j_count = i_count 
      if (np_sharing == 1) then 
        iproc   = iwork(1) 
        do 
          j_count = j_count + 1 
          if (j_count-i_count >= nb) exit
          if (j_count > nrow) exit
          call parts(j_count,nrow,np,iwrk2, np_sharing)
          if (np_sharing /= 1 ) exit
          if (iwrk2(1) /= iproc ) exit
        end do
      end if
    else if (use_vg) then 
      np_sharing = 1
      j_count = i_count 
      iproc   = vg(i_count)
      iwork(1:np_sharing) = iproc
      do 
        j_count = j_count + 1 
        if (j_count-i_count >= nb) exit
        if (j_count > nrow) exit
        if (vg(j_count) /= iproc ) exit
      end do
    else if (use_vsz) then
      np_sharing = 1
      j_count = i_count 
      iwork(1:np_sharing) = iproc
      do 
        j_count = j_count + 1 
        if (j_count-i_count >= nb) exit
        if (j_count > nrow) exit
        if (j_count > lastigp) exit
      end do      
    end if

    ! now we should insert rows i_count..j_count-1
    nnr = j_count - i_count
    
    if (iam == root) then
      
      ll = 0
      do i= i_count, j_count-1
        call a_glob%csget(i,i,nz,&
             & irow,icol,val,info,nzin=ll,append=.true.)
        if (info /= psb_success_) then            
          if (nz >min(size(irow(ll+1:)),size(icol(ll+1:)),size(val(ll+1:)))) then 
            write(psb_err_unit,*) 'Allocation failure? This should not happen!'
          end if
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
        ll = ll + nz
      end do

      do k_count = 1, np_sharing
        iproc = iwork(k_count)
        
        if (iproc == iam) then
          call psb_spins(ll,irow,icol,val,a,desc_a,info)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_spins'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        else
          call psb_snd(ctxt,nnr,iproc)
          call psb_snd(ctxt,ll,iproc)
          call psb_snd(ctxt,irow(1:ll),iproc)
          call psb_snd(ctxt,icol(1:ll),iproc)
          call psb_snd(ctxt,val(1:ll),iproc)
          call psb_rcv(ctxt,ll,iproc)
        endif
      end do
    else if (iam /= root) then
      
      do k_count = 1, np_sharing
        iproc = iwork(k_count)
        if (iproc == iam) then
          call psb_rcv(ctxt,nnr,root)
          call psb_rcv(ctxt,ll,root)
          if (ll > size(irow)) then 
            write(psb_err_unit,*) iam,'need to reallocate ',ll
            deallocate(val,irow,icol)
            allocate(val(ll),irow(ll),icol(ll),stat=info)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='Allocate'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
            
          endif
          call psb_rcv(ctxt,irow(1:ll),root)
          call psb_rcv(ctxt,icol(1:ll),root)
          call psb_rcv(ctxt,val(1:ll),root)
          call psb_snd(ctxt,ll,root)
          call psb_spins(ll,irow,icol,val,a,desc_a,info)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psspins'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        endif
      end do
    endif
    i_count = j_count
    if ((use_vsz).and.(j_count <= nrow)) then 
      if (j_count > lastigp) then
        iproc = iproc + 1
        lastigp = lastigp + vsz(iproc+1)
      end if
    end if
  end do

  call psb_barrier(ctxt)
  t0 = psb_wtime()
  call psb_cdasb(desc_a,info)     
  t1 = psb_wtime()
  if(info /= psb_success_)then
    info=psb_err_from_subroutine_
    ch_err='psb_cdasb'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_barrier(ctxt)
  t2 = psb_wtime()
  call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=fmt,mold=mold)     
  t3 = psb_wtime()
  if(info /= psb_success_)then
    info=psb_err_from_subroutine_
    ch_err='psb_spasb'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if


  if (iam == root) then 
    write(psb_out_unit,*) 'descriptor assembly: ',t1-t0
    write(psb_out_unit,*) 'sparse matrix assembly: ',t3-t2
  end if


  
  deallocate(val,irow,icol,iwork,iwrk2,stat=info)
  if(info /= psb_success_)then
    info=psb_err_from_subroutine_
    ch_err='deallocate'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (iam == root) write (*, fmt = *) 'end matdist'     

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_cmatdist



subroutine psb_lcmatdist(a_glob, a, ctxt, desc_a,&
     & info, parts, vg, vsz, inroot,fmt,mold)
  !
  ! an utility subroutine to distribute a matrix among processors
  ! according to a user defined data distribution, using
  ! sparse matrix subroutines.
  !
  !  type(psb_cspmat)                       :: a_glob
  !     on entry: this contains the global sparse matrix as follows:
  !
  !  type(psb_cspmat_type)                            :: a
  !     on exit : this will contain the local sparse matrix.
  !
  !       interface parts
  !         !   .....user passed subroutine.....
  !         subroutine parts(global_indx,n,np,pv,nv)
  !           implicit none
  !           integer(psb_ipk_), intent(in)  :: global_indx, n, np
  !           integer(psb_ipk_), intent(out) :: nv
  !           integer(psb_ipk_), intent(out) :: pv(*)
  !
  !       end subroutine parts
  !       end interface
  !     on entry:  subroutine providing user defined data distribution.
  !        for each global_indx the subroutine should return
  !        the list  pv of all processes owning the row with
  !        that index; the list will contain nv entries.
  !        usually nv=1; if nv >1 then we have an overlap in the data
  !        distribution.
  !
  !  integer(psb_ipk_) :: ctxt
  !     on entry: the PSBLAS parallel environment context.
  !
  !  type (desc_type)                  :: desc_a
  !     on exit : the updated array descriptor
  !
  !  integer(psb_ipk_), optional    :: inroot
  !     on entry: specifies processor holding a_glob. default: 0
  !     on exit : unchanged.
  !
  use psb_base_mod
  use psb_mat_mod
  implicit none

  ! parameters
  type(psb_lcspmat_type)      :: a_glob
  type(psb_ctxt_type) :: ctxt
  type(psb_cspmat_type)      :: a
  type(psb_desc_type)        :: desc_a
  integer(psb_ipk_), intent(out)       :: info
  integer(psb_ipk_), optional       :: inroot
  character(len=*), optional :: fmt
  class(psb_c_base_sparse_mat), optional :: mold
  procedure(psb_parts), optional  :: parts
  integer(psb_ipk_), optional     :: vg(:) 
  integer(psb_ipk_), optional     :: vsz(:)
  
  ! local variables
  logical           :: use_parts, use_vg, use_vsz
  integer(psb_ipk_) :: np, iam, np_sharing, root, iproc
  integer(psb_ipk_) :: err_act, il, inz
  integer(psb_lpk_) :: k_count, liwork,  nnzero, nrhs,&
       & i, ll, nz, isize, nnr, err
  integer(psb_lpk_) :: i_count, j_count, nrow, ncol, ig, lastigp
  integer(psb_ipk_), allocatable  :: iwork(:), iwrk2(:)
  integer(psb_lpk_), allocatable  :: irow(:),icol(:)
  complex(psb_spk_), allocatable    :: val(:)
  integer(psb_ipk_), parameter    :: nb=30
  real(psb_dpk_)              :: t0, t1, t2, t3, t4, t5
  character(len=20)           :: name, ch_err

  info = psb_success_
  err  = 0
  name = 'psb_c_mat_dist'
  call psb_erractionsave(err_act)

  ! executable statements    
  if (present(inroot)) then
    root = inroot
  else
    root = psb_root_
  end if
  call psb_info(ctxt, iam, np)     
  if (iam == root) then
    nrow = a_glob%get_nrows()
    ncol = a_glob%get_ncols()
    if (nrow /= ncol) then
      write(psb_err_unit,*) 'a rectangular matrix ? ',nrow,ncol
      info=-1
      call psb_errpush(info,name)
      goto 9999
    endif
    nnzero = a_glob%get_nzeros()
    nrhs   = 1
  endif

  use_parts = present(parts)
  use_vg    = present(vg)
  use_vsz   = present(vsz)
  if (count((/ use_parts, use_vg, use_vsz /)) /= 1) then 
    info=psb_err_no_optional_arg_
    call psb_errpush(info,name,a_err=" vg, vsz, parts")
    goto 9999 
  endif
  
  ! broadcast informations to other processors
  call psb_bcast(ctxt,nrow, root)
  call psb_bcast(ctxt,ncol, root)
  call psb_bcast(ctxt,nnzero, root)
  call psb_bcast(ctxt,nrhs, root)
  liwork = max(np, nrow + ncol)
  allocate(iwork(liwork), iwrk2(np),stat = info)
  if (info /= psb_success_) then
    info=psb_err_alloc_request_
    call psb_errpush(info,name,l_err=(/liwork/),a_err='integer')
    goto 9999
  endif
  if (iam == root) then
    write (*, fmt = *) 'start matdist',root, size(iwork),&
         &nrow, ncol, nnzero,nrhs
  endif
  if (use_parts) then 
    call psb_cdall(ctxt,desc_a,info,mg=nrow,parts=parts)
  else if (use_vg) then 
    call psb_cdall(ctxt,desc_a,info,vg=vg)
  else if (use_vsz) then 
    call psb_cdall(ctxt,desc_a,info,nl=vsz(iam+1))
  else
    info = -1
  end if
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_cdall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  inz = ((nnzero+np-1)/np)
  call psb_spall(a,desc_a,info,nnz=inz)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_spall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  isize = 3*nb*max(((nnzero+nrow)/nrow),nb)
  allocate(val(isize),irow(isize),icol(isize),stat=info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='Allocate'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  i_count   = 1
  if (use_vsz) then
    iproc = 0
    lastigp = vsz(iproc+1)
  end if
  do while (i_count <= nrow)

    if (use_parts) then 
      call parts(i_count,nrow,np,iwork, np_sharing)
      !
      ! np_sharing allows for overlap in the data distribution.
      ! If an index is overlapped, then we have to send its row
      ! to multiple processes. NOTE: we are assuming the output 
      ! from PARTS is coherent, otherwise a deadlock is the most
      ! likely outcome.
      ! 
      j_count = i_count 
      if (np_sharing == 1) then 
        iproc   = iwork(1) 
        do 
          j_count = j_count + 1 
          if (j_count-i_count >= nb) exit
          if (j_count > nrow) exit
          call parts(j_count,nrow,np,iwrk2, np_sharing)
          if (np_sharing /= 1 ) exit
          if (iwrk2(1) /= iproc ) exit
        end do
      end if
    else if (use_vg) then 
      np_sharing = 1
      j_count = i_count 
      iproc   = vg(i_count)
      iwork(1:np_sharing) = iproc
      do 
        j_count = j_count + 1 
        if (j_count-i_count >= nb) exit
        if (j_count > nrow) exit
        if (vg(j_count) /= iproc ) exit
      end do
    else if (use_vsz) then
      np_sharing = 1
      j_count = i_count 
      iwork(1:np_sharing) = iproc
      do 
        j_count = j_count + 1 
        if (j_count-i_count >= nb) exit
        if (j_count > nrow) exit
        if (j_count > lastigp) exit
      end do      
    end if

    ! now we should insert rows i_count..j_count-1
    nnr = j_count - i_count
    
    if (iam == root) then
      
      ll = 0
      do i= i_count, j_count-1
        call a_glob%csget(i,i,nz,&
             & irow,icol,val,info,nzin=ll,append=.true.)
        if (info /= psb_success_) then            
          if (nz >min(size(irow(ll+1:)),size(icol(ll+1:)),size(val(ll+1:)))) then 
            write(psb_err_unit,*) 'Allocation failure? This should not happen!'
          end if
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
        ll = ll + nz
      end do

      do k_count = 1, np_sharing
        iproc = iwork(k_count)
        
        if (iproc == iam) then
          il = ll
          call psb_spins(il,irow,icol,val,a,desc_a,info)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psb_spins'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        else
          call psb_snd(ctxt,nnr,iproc)
          call psb_snd(ctxt,ll,iproc)
          call psb_snd(ctxt,irow(1:ll),iproc)
          call psb_snd(ctxt,icol(1:ll),iproc)
          call psb_snd(ctxt,val(1:ll),iproc)
          call psb_rcv(ctxt,ll,iproc)
        endif
      end do
    else if (iam /= root) then
      
      do k_count = 1, np_sharing
        iproc = iwork(k_count)
        if (iproc == iam) then
          call psb_rcv(ctxt,nnr,root)
          call psb_rcv(ctxt,ll,root)
          if (ll > size(irow)) then 
            write(psb_err_unit,*) iam,'need to reallocate ',ll
            deallocate(val,irow,icol)
            allocate(val(ll),irow(ll),icol(ll),stat=info)
            if(info /= psb_success_) then
              info=psb_err_from_subroutine_
              ch_err='Allocate'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
            
          endif
          call psb_rcv(ctxt,irow(1:ll),root)
          call psb_rcv(ctxt,icol(1:ll),root)
          call psb_rcv(ctxt,val(1:ll),root)
          call psb_snd(ctxt,ll,root)
          il = ll 
          call psb_spins(il,irow,icol,val,a,desc_a,info)
          if(info /= psb_success_) then
            info=psb_err_from_subroutine_
            ch_err='psspins'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        endif
      end do
    endif
    i_count = j_count
    if ((use_vsz).and.(j_count <= nrow)) then 
      if (j_count > lastigp) then
        iproc = iproc + 1
        lastigp = lastigp + vsz(iproc+1)
      end if
    end if
  end do

  call psb_barrier(ctxt)
  t0 = psb_wtime()
  call psb_cdasb(desc_a,info)     
  t1 = psb_wtime()
  if(info /= psb_success_)then
    info=psb_err_from_subroutine_
    ch_err='psb_cdasb'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_barrier(ctxt)
  t2 = psb_wtime()
  call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=fmt,mold=mold)     
  t3 = psb_wtime()
  if(info /= psb_success_)then
    info=psb_err_from_subroutine_
    ch_err='psb_spasb'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if


  if (iam == root) then 
    write(psb_out_unit,*) 'descriptor assembly: ',t1-t0
    write(psb_out_unit,*) 'sparse matrix assembly: ',t3-t2
  end if


  
  deallocate(val,irow,icol,iwork,iwrk2,stat=info)
  if(info /= psb_success_)then
    info=psb_err_from_subroutine_
    ch_err='deallocate'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (iam == root) write (*, fmt = *) 'end matdist'     

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_lcmatdist
