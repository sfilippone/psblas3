!   
!                Parallel Sparse BLAS  version 3.4
!      (C) Copyright 2006, 2010, 2015
!                         Salvatore Filippone    University of Rome Tor Vergata
!                         Alfredo Buttari        CNRS-IRIT, Toulouse
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
subroutine psb_zmatdist(a_glob, a, ictxt, desc_a,&
     & info, parts, v, inroot,fmt,mold)
  !
  ! an utility subroutine to distribute a matrix among processors
  ! according to a user defined data distribution, using
  ! sparse matrix subroutines.
  !
  !  type(psb_zspmat)                       :: a_glob
  !     on entry: this contains the global sparse matrix as follows:
  !
  !  type(psb_zspmat_type)                            :: a
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
  !  integer(psb_ipk_) :: ictxt
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
  type(psb_zspmat_type)      :: a_glob
  integer(psb_ipk_) :: ictxt
  type(psb_zspmat_type)      :: a
  type(psb_desc_type)        :: desc_a
  integer(psb_ipk_), intent(out)       :: info
  integer(psb_ipk_), optional       :: inroot
  character(len=*), optional :: fmt
  class(psb_z_base_sparse_mat), optional :: mold
  procedure(psb_parts), optional  :: parts
  integer(psb_ipk_), optional     :: v(:)

  ! local variables
  logical           :: use_parts, use_v
  integer(psb_ipk_) :: np, iam, np_sharing
  integer(psb_ipk_) :: i_count, j_count,&
       & k_count, root, liwork, nrow, ncol, nnzero, nrhs,&
       & i, ll, nz, isize, iproc, nnr, err, err_act, int_err(5)
  integer(psb_ipk_), allocatable  :: iwork(:)
  integer(psb_ipk_), allocatable  :: irow(:),icol(:)
  complex(psb_dpk_), allocatable    :: val(:)
  integer(psb_ipk_), parameter    :: nb=30
  real(psb_dpk_)              :: t0, t1, t2, t3, t4, t5
  character(len=20)           :: name, ch_err

  info = psb_success_
  err  = 0
  name = 'psb_z_mat_dist'
  call psb_erractionsave(err_act)

  ! executable statements    
  if (present(inroot)) then
    root = inroot
  else
    root = psb_root_
  end if
  call psb_info(ictxt, iam, np)     
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
  use_v     = present(v)
  if (count((/ use_parts, use_v /)) /= 1) then 
    info=psb_err_no_optional_arg_
    call psb_errpush(info,name,a_err=" v, parts")
    goto 9999 
  endif
  
  ! broadcast informations to other processors
  call psb_bcast(ictxt,nrow, root)
  call psb_bcast(ictxt,ncol, root)
  call psb_bcast(ictxt,nnzero, root)
  call psb_bcast(ictxt,nrhs, root)
  liwork = max(np, nrow + ncol)
  allocate(iwork(liwork), stat = info)
  if (info /= psb_success_) then
    info=psb_err_alloc_request_
    int_err(1)=liwork
    call psb_errpush(info,name,i_err=int_err,a_err='integer')
    goto 9999
  endif
  if (iam == root) then
    write (*, fmt = *) 'start matdist',root, size(iwork),&
         &nrow, ncol, nnzero,nrhs
  endif
  if (use_parts) then 
    call psb_cdall(ictxt,desc_a,info,mg=nrow,parts=parts)
  else if (use_v) then 
    call psb_cdall(ictxt,desc_a,info,vg=v)
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
          call parts(j_count,nrow,np,iwork, np_sharing)
          if (np_sharing /= 1 ) exit
          if (iwork(1) /= iproc ) exit
        end do
      end if
    else 
      np_sharing = 1
      j_count = i_count 
      iproc   = v(i_count)
      iwork(1:np_sharing) = iproc
      do 
        j_count = j_count + 1 
        if (j_count-i_count >= nb) exit
        if (j_count > nrow) exit
        if (v(j_count) /= iproc ) exit
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
          call psb_snd(ictxt,nnr,iproc)
          call psb_snd(ictxt,ll,iproc)
          call psb_snd(ictxt,irow(1:ll),iproc)
          call psb_snd(ictxt,icol(1:ll),iproc)
          call psb_snd(ictxt,val(1:ll),iproc)
          call psb_rcv(ictxt,ll,iproc)
        endif
      end do
    else if (iam /= root) then
      
      do k_count = 1, np_sharing
        iproc = iwork(k_count)
        if (iproc == iam) then
          call psb_rcv(ictxt,nnr,root)
          call psb_rcv(ictxt,ll,root)
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
          call psb_rcv(ictxt,irow(1:ll),root)
          call psb_rcv(ictxt,icol(1:ll),root)
          call psb_rcv(ictxt,val(1:ll),root)
          call psb_snd(ictxt,ll,root)
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
  end do

  call psb_barrier(ictxt)
  t0 = psb_wtime()
  call psb_cdasb(desc_a,info)     
  t1 = psb_wtime()
  if(info /= psb_success_)then
    info=psb_err_from_subroutine_
    ch_err='psb_cdasb'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  call psb_barrier(ictxt)
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


  
  deallocate(val,irow,icol,stat=info)
  if(info /= psb_success_)then
    info=psb_err_from_subroutine_
    ch_err='deallocate'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  deallocate(iwork)   
  if (iam == root) write (*, fmt = *) 'end matdist'     

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

end subroutine psb_zmatdist
