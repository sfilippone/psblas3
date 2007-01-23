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
module psb_mat_dist_mod

  interface psb_matdist
     module procedure dmatdistf, dmatdistv, zmatdistf, zmatdistv
  end interface

contains

  subroutine dmatdistf (a_glob, a, parts, ictxt, desc_a,&
       & b_glob, b, info, inroot,fmt)
    !
    ! an utility subroutine to distribute a matrix among processors
    ! according to a user defined data distribution, using 
    ! sparse matrix subroutines.
    !
    !  type(d_spmat)                            :: a_glob
    !     on entry: this contains the global sparse matrix as follows:
    !        a%fida =='csr'
    !        a%aspk for coefficient values
    !        a%ia1  for column indices
    !        a%ia2  for row pointers
    !        a%m    for number of global matrix rows
    !        a%k    for number of global matrix columns
    !     on exit : undefined, with unassociated pointers.
    !
    !  type(d_spmat)                            :: a
    !     on entry: fresh variable.
    !     on exit : this will contain the local sparse matrix.
    !
    !       interface parts
    !         !   .....user passed subroutine.....
    !         subroutine parts(global_indx,n,np,pv,nv)
    !           implicit none
    !           integer, intent(in)  :: global_indx, n, np
    !           integer, intent(out) :: nv
    !           integer, intent(out) :: pv(*)
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
    !  integer                                  :: ictxt
    !     on entry: blacs context.
    !     on exit : unchanged.
    !
    !  type (desc_type)                  :: desc_a
    !     on entry: fresh variable.
    !     on exit : the updated array descriptor
    !
    !  real(kind(1.d0)), optional      :: b_glob(:)
    !     on entry: this contains right hand side.
    !     on exit :
    !
    !  real(kind(1.d0)), allocatable, optional      :: b(:)
    !     on entry: fresh variable.
    !     on exit : this will contain the local right hand side.
    !
    !  integer, optional    :: inroot
    !     on entry: specifies processor holding a_glob. default: 0
    !     on exit : unchanged.
    !
    use psb_base_mod
    implicit none

    ! parameters
    type(psb_dspmat_type)      :: a_glob
    real(kind(1.d0))           :: b_glob(:)
    integer                    :: ictxt
    type(psb_dspmat_type)      :: a
    real(kind(1.d0)), allocatable  :: b(:)
    type (psb_desc_type)       :: desc_a
    integer, intent(out)       :: info
    integer, optional          :: inroot
    character(len=5), optional :: fmt
    interface 

      !   .....user passed subroutine.....
      subroutine parts(global_indx,n,np,pv,nv)
        implicit none
        integer, intent(in)  :: global_indx, n, np
        integer, intent(out) :: nv
        integer, intent(out) :: pv(*) 
      end subroutine parts
    end interface

    ! local variables
    integer                     :: np, iam 
    integer                     :: ircode, length_row, i_count, j_count,&
         & k_count, blockdim, root, liwork, nrow, ncol, nnzero, nrhs,&
         & i,j,k, ll, nz, isize, iproc, nnr, err, err_act, int_err(5)
    integer, allocatable            :: iwork(:)
    character                   :: afmt*5, atyp*5
    integer, allocatable          :: irow(:),icol(:)
    real(kind(1.d0)), allocatable :: val(:)
    integer, parameter          :: nb=30
    real(kind(1.d0))            :: t0, t1, t2, t3, t4, t5
    character(len=20)           :: name, ch_err

    info = 0
    err  = 0
    name = 'mat_distf'
    call psb_erractionsave(err_act)

    ! executable statements    
    if (present(inroot)) then
      root = inroot
    else
      root = 0
    end if
    call psb_info(ictxt, iam, np)     

    if (iam == root) then
      ! extract information from a_glob
      if (a_glob%fida.ne. 'CSR') then
        info=135
        ch_err='CSR'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      endif
      nrow = a_glob%m
      ncol = a_glob%k
      if (nrow /= ncol) then
        write(0,*) 'a rectangular matrix ? ',nrow,ncol
        info=-1
        call psb_errpush(info,name)
        goto 9999
      endif
      nnzero = size(a_glob%aspk)
      nrhs   = 1
      ! broadcast informations to other processors
    endif
    call psb_bcast(ictxt, nrow,root)
    call psb_bcast(ictxt, ncol,root)
    call psb_bcast(ictxt, nnzero,root)
    call psb_bcast(ictxt, nrhs,root)

    liwork = max(np, nrow + ncol)
    allocate(iwork(liwork), stat = info)
    if (info /= 0) then
      info=2025
      int_err(1)=liwork
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif
    if (iam == root) then
      write (*, fmt = *) 'start matdist',root, size(iwork),&
           &nrow, ncol, nnzero,nrhs
    endif
    call psb_cdall(ictxt,desc_a,info,mg=nrow,parts=parts)
    if(info/=0) then
      info=4010
      ch_err='psb_cdall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_spall(a,desc_a,info,nnz=nnzero/np)
    if(info/=0) then
      info=4010
      ch_err='psb_psspall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    call psb_geall(b,desc_a,info)   
    if(info/=0) then
      info=4010
      ch_err='psb_psdsall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    isize = max(3*nb,ncol)


    allocate(val(nb*ncol),irow(nb*ncol),icol(nb*ncol),stat=info)
    if(info/=0) then
      info=4010
      ch_err='Allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    i_count   = 1

    do while (i_count.le.nrow)

      call parts(i_count,nrow,np,iwork, length_row)

      if (length_row.eq.1) then 
        j_count = i_count 
        iproc   = iwork(1) 
        do 
          j_count = j_count + 1 
          if (j_count-i_count >= nb) exit
          if (j_count > nrow) exit
          call parts(j_count,nrow,np,iwork, length_row)
          if (length_row /= 1 ) exit
          if (iwork(1) /= iproc ) exit
        end do

        ! now we should insert rows i_count..j_count-1
        nnr = j_count - i_count

        if (iam == root) then

          ll = 0
          do i= i_count, j_count-1
            call psb_sp_getrow(i,a_glob,nz,&
                 & irow(ll+1:),icol(ll+1:),val(ll+1:), info)
            if (info /= 0) then            
              if (nz >min(size(irow(ll+1:)),size(icol(ll+1:)),size(val(ll+1:)))) then 
                write(0,*) 'Allocation failure? This should not happen!'
              end if
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
            ll = ll + nz
          end do

          if (iproc == iam) then
            call psb_spins(ll,irow,icol,val,a,desc_a,info)
            if(info/=0) then
              info=4010
              ch_err='psb_spins'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
            call psb_geins(nnr,(/(i,i=i_count,j_count-1)/),&
                 & b_glob(i_count:j_count-1),b,desc_a,info)
            if(info/=0) then
              info=4010
              ch_err='psb_ins'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
          else
            call psb_snd(ictxt,nnr,iproc)
            call psb_snd(ictxt,ll,iproc)
            call psb_snd(ictxt,irow(1:ll),iproc)
            call psb_snd(ictxt,icol(1:ll),iproc)
            call psb_snd(ictxt,val(1:ll),iproc)
            call psb_snd(ictxt,b_glob(i_count:j_count-1),iproc)
            call psb_rcv(ictxt,ll,iproc)
          endif
        else if (iam /= root) then

          if (iproc == iam) then
            call psb_rcv(ictxt,nnr,root)
            call psb_rcv(ictxt,ll,root)
            if (ll > size(irow)) then 
              write(0,*) iam,'need to reallocate ',ll
              deallocate(val,irow,icol)
              allocate(val(ll),irow(ll),icol(ll),stat=info)
              if(info/=0) then
                info=4010
                ch_err='Allocate'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if

            endif
            call psb_rcv(ictxt,irow(1:ll),root)
            call psb_rcv(ictxt,icol(1:ll),root)
            call psb_rcv(ictxt,val(1:ll),root)
            call psb_rcv(ictxt,b_glob(i_count:i_count+nnr-1),root)
            call psb_snd(ictxt,ll,root)
            call psb_spins(ll,irow,icol,val,a,desc_a,info)
            if(info/=0) then
              info=4010
              ch_err='psspins'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
            call psb_geins(nnr,(/(i,i=i_count,i_count+nnr-1)/),&
                 & b_glob(i_count:i_count+nnr-1),b,desc_a,info)
            if(info/=0) then
              info=4010
              ch_err='psdsins'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
          endif
        endif

        i_count = j_count

      else
        write(0,*) iam,'unexpected turn'
        ! here processors are counted 1..np
        do j_count = 1, length_row
          k_count = iwork(j_count)
          if (iam == root) then

            ll = 0
            do i= i_count, i_count
              call psb_sp_getrow(i,a_glob,nz,&
                   & irow(ll+1:),icol(ll+1:),val(ll+1:), info)
              if (info /= 0) then            
                if (nz >min(size(irow(ll+1:)),size(icol(ll+1:)),size(val(ll+1:)))) then 
                  write(0,*) 'Allocation failure? This should not happen!'
                end if
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
              ll = ll + nz
            end do

            if (k_count == iam) then

              call psb_spins(ll,irow,icol,val,a,desc_a,info)
              if(info/=0) then
                info=4010
                ch_err='psspins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
              call psb_geins(1,(/i_count/),b_glob(i_count:i_count),&
                   & b,desc_a,info)
              if(info/=0) then
                info=4010
                ch_err='psdsins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
            else
              call psb_snd(ictxt,ll,k_count)
              call psb_snd(ictxt,irow(1:ll),k_count)
              call psb_snd(ictxt,icol(1:ll),k_count)
              call psb_snd(ictxt,val(1:ll),k_count)
              call psb_snd(ictxt,b_glob(i_count),k_count)
              call psb_rcv(ictxt,ll,k_count)
            endif
          else if (iam /= root) then
            if (k_count == iam) then

              call psb_rcv(ictxt,ll,root)
              call psb_rcv(ictxt,irow(1:ll),root)
              call psb_rcv(ictxt,icol(1:ll),root)
              call psb_rcv(ictxt,val(1:ll),root)
              call psb_rcv(ictxt,b_glob(i_count),root)
              call psb_snd(ictxt,ll,root)
              call psb_spins(ll,irow,icol,val,a,desc_a,info)
              if(info/=0) then
                info=4010
                ch_err='psspins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
              call psb_geins(1,(/i_count/),b_glob(i_count:i_count),&
                   & b,desc_a,info)
              if(info/=0) then
                info=4010
                ch_err='psdsins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
            endif
          endif
        end do
        i_count = i_count + 1
      endif
    end do

    if (present(fmt)) then 
      afmt=fmt
    else
      afmt = 'CSR'
    endif
    call psb_barrier(ictxt)
    t0 = psb_wtime()
    call psb_cdasb(desc_a,info)     
    t1 = psb_wtime()
    if(info/=0)then
      info=4010
      ch_err='psb_cdasb'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_barrier(ictxt)
    t2 = psb_wtime()
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)     
    t3 = psb_wtime()
    if(info/=0)then
      info=4010
      ch_err='psb_spasb'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if


    if (iam == root) then 
      write(*,*) 'descriptor assembly: ',t1-t0
      write(*,*) 'sparse matrix assembly: ',t3-t2
    end if



    call psb_geasb(b,desc_a,info)     
    if(info/=0)then
      info=4010
      ch_err='psdsasb'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    deallocate(val,irow,icol,stat=info)
    if(info/=0)then
      info=4010
      ch_err='deallocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    deallocate(iwork)   
    if (iam == root) write (*, fmt = *) 'end matdist'     

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return

  end subroutine dmatdistf


  subroutine dmatdistv (a_glob, a, v, ictxt, desc_a,&
       & b_glob, b, info, inroot,fmt)
    !
    ! an utility subroutine to distribute a matrix among processors
    ! according to a user defined data distribution, using 
    ! sparse matrix subroutines.
    !
    !  type(d_spmat)                            :: a_glob
    !     on entry: this contains the global sparse matrix as follows:
    !        a%fida =='csr'
    !        a%aspk for coefficient values
    !        a%ia1  for column indices
    !        a%ia2  for row pointers
    !        a%m    for number of global matrix rows
    !        a%k    for number of global matrix columns
    !     on exit : undefined, with unassociated pointers.
    !
    !  type(d_spmat)                            :: a
    !     on entry: fresh variable.
    !     on exit : this will contain the local sparse matrix.
    !
    !       interface parts
    !         !   .....user passed subroutine.....
    !         subroutine parts(global_indx,n,np,pv,nv)
    !           implicit none
    !           integer, intent(in)  :: global_indx, n, np
    !           integer, intent(out) :: nv
    !           integer, intent(out) :: pv(*)
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
    !  integer                                  :: ictxt
    !     on entry: blacs context.
    !     on exit : unchanged.
    !
    !  type (desc_type)                  :: desc_a
    !     on entry: fresh variable.
    !     on exit : the updated array descriptor
    !
    !  real(kind(1.d0)),  optional      :: b_glob(:)
    !     on entry: this contains right hand side.
    !     on exit :
    !
    !  real(kind(1.d0)), allocatable, optional      :: b(:)
    !     on entry: fresh variable.
    !     on exit : this will contain the local right hand side.
    !
    !  integer, optional    :: inroot
    !     on entry: specifies processor holding a_glob. default: 0
    !     on exit : unchanged.
    !
    use psb_base_mod
    implicit none   ! parameters
    type(psb_dspmat_type)      :: a_glob
    real(kind(1.d0))           :: b_glob(:)
    integer                    :: ictxt, v(:)
    type(psb_dspmat_type)      :: a
    real(kind(1.d0)), allocatable  :: b(:)
    type (psb_desc_type)       :: desc_a
    integer, intent(out)       :: info
    integer, optional          :: inroot
    character(len=5), optional :: fmt

    integer                     :: np, iam
    integer                     :: ircode, length_row, i_count, j_count,&
         & k_count, blockdim, root, liwork, nrow, ncol, nnzero, nrhs,&
         & i,j,k, ll, nz, isize, iproc, nnr, err, err_act, int_err(5)
    integer, allocatable            :: iwork(:)
    character                   :: afmt*5, atyp*5
    integer, allocatable          :: irow(:),icol(:)
    real(kind(1.d0)), allocatable :: val(:)
    integer, parameter          :: nb=30
    real(kind(1.d0))            :: t0, t1, t2, t3, t4, t5
    character(len=20)  :: name, ch_err

    info = 0
    err  = 0
    name = 'mat_distv'
    call psb_erractionsave(err_act)

    ! executable statements    
    if (present(inroot)) then
      root = inroot
    else
      root = 0
    end if

    call psb_info(ictxt, iam, np)     
    if (iam == root) then
      ! extract information from a_glob
      if (a_glob%fida.ne. 'CSR') then
        info=135
        ch_err='CSR'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      endif

      nrow = a_glob%m
      ncol = a_glob%k
      if (nrow /= ncol) then
        write(0,*) 'a rectangular matrix ? ',nrow,ncol
        info=-1
        call psb_errpush(info,name)
        goto 9999
      endif

      nnzero = size(a_glob%aspk)
      nrhs   = 1
    end if
    ! broadcast informations to other processors
    call psb_bcast(ictxt,nrow, root)
    call psb_bcast(ictxt,ncol, root)
    call psb_bcast(ictxt,nnzero, root)
    call psb_bcast(ictxt,nrhs, root)
    liwork = max(np, nrow + ncol)
    allocate(iwork(liwork), stat = info)
    if (info /= 0) then
      write(0,*) 'matdist allocation failed'
      info=2025
      int_err(1)=liwork
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif

    call psb_cdall(ictxt,desc_a,info,vg=v)
    if(info/=0) then
      info=4010
      ch_err='psb_cdall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_spall(a,desc_a,info,nnz=((nnzero+np-1)/np))
    if(info/=0) then
      info=4010
      ch_err='psb_psspall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    call psb_geall(b,desc_a,info)   
    if(info/=0) then
      info=4010
      ch_err='psb_psdsall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    isize = max(3*nb,ncol)


    allocate(val(nb*ncol),irow(nb*ncol),icol(nb*ncol),stat=info)
    if(info/=0) then
      info=4010
      ch_err='Allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    i_count   = 1

    do while (i_count <= nrow)

      j_count = i_count 
      iproc   = v(i_count)

      do 
        j_count = j_count + 1 
        if (j_count-i_count >= nb) exit
        if (j_count > nrow) exit
        if (v(j_count) /= iproc ) exit
      end do

      ! now we should insert rows i_count..j_count-1
      nnr = j_count - i_count

      if (iam == root) then
        
        ll = 0
        do i= i_count, j_count-1
          call psb_sp_getrow(i,a_glob,nz,&
               & irow(ll+1:),icol(ll+1:),val(ll+1:), info)
          if (info /= 0) then            
            if (nz >min(size(irow(ll+1:)),size(icol(ll+1:)),size(val(ll+1:)))) then 
              write(0,*) 'Allocation failure? This should not happen!'
            end if
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
          ll = ll + nz
        end do

        if (iproc == iam) then
          call psb_spins(ll,irow,icol,val,a,desc_a,info)
          if(info/=0) then
            info=4010
            ch_err='psb_spins'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if

          call psb_geins(nnr,(/(i,i=i_count,j_count-1)/),b_glob(i_count:j_count-1),&
               & b,desc_a,info)
          if(info/=0) then
            info=4010
            ch_err='dsins'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        else
          call psb_snd(ictxt,nnr,iproc)
          call psb_snd(ictxt,ll,iproc)
          call psb_snd(ictxt,irow(1:ll),iproc)
          call psb_snd(ictxt,icol(1:ll),iproc)
          call psb_snd(ictxt,val(1:ll),iproc)
          call psb_snd(ictxt,b_glob(i_count:j_count-1),iproc)
          call psb_rcv(ictxt,ll,iproc)
        endif
      else if (iam /= root) then

        if (iproc == iam) then
          call psb_rcv(ictxt,nnr,root)
          call psb_rcv(ictxt,ll,root)
          if (ll > size(val)) then 
            write(0,*) iam,'need to reallocate ',ll
            deallocate(val,irow,icol)
            allocate(val(ll),irow(ll),icol(ll),stat=info)
            if(info/=0) then
              info=4010
              ch_err='Allocate'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
          endif
          call psb_rcv(ictxt,irow(1:ll),root)
          call psb_rcv(ictxt,icol(1:ll),root)
          call psb_rcv(ictxt,val(1:ll),root)
          call psb_rcv(ictxt,b_glob(i_count:i_count+nnr-1),root)
          call psb_snd(ictxt,ll,root)

          call psb_spins(ll,irow,icol,val,a,desc_a,info)
          if(info/=0) then
            info=4010
            ch_err='spins'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
          call psb_geins(nnr,(/(i,i=i_count,i_count+nnr-1)/),&
               & b_glob(i_count:i_count+nnr-1),b,desc_a,info)
          if(info/=0) then
            info=4010
            ch_err='psdsins'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        endif
      endif
      i_count = j_count

    end do

    ! default storage format for sparse matrix; we do not
    ! expect duplicated entries.

    if (present(fmt)) then  
      afmt=fmt
    else
      afmt = 'CSR'
    endif
    call psb_barrier(ictxt)
    t0 = psb_wtime()
    call psb_cdasb(desc_a,info)     
    t1 = psb_wtime()
    if(info/=0)then
      info=4010
      ch_err='psb_cdasb'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_barrier(ictxt)
    t2 = psb_wtime()
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)     
    t3 = psb_wtime()
    if(info/=0)then
      info=4010
      ch_err='psb_spasb'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_geasb(b,desc_a,info)

    if (iam == root) then 
      write(*,'("Descriptor assembly   : ",es10.4)')t1-t0
      write(*,'("Sparse matrix assembly: ",es10.4)')t3-t2
    end if

    if(info/=0)then
      info=4010
      ch_err='psdsasb'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    deallocate(iwork)   

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return

  end subroutine dmatdistv


  subroutine zmatdistf (a_glob, a, parts, ictxt, desc_a,&
       & b_glob, b, info, inroot,fmt)
    !
    ! an utility subroutine to distribute a matrix among processors
    ! according to a user defined data distribution, using
    ! sparse matrix subroutines.
    !
    !  type(d_spmat)                            :: a_glob
    !     on entry: this contains the global sparse matrix as follows:
    !        a%fida =='csr'
    !        a%aspk for coefficient values
    !        a%ia1  for column indices
    !        a%ia2  for row pointers
    !        a%m    for number of global matrix rows
    !        a%k    for number of global matrix columns
    !     on exit : undefined, with unassociated pointers.
    !
    !  type(d_spmat)                            :: a
    !     on entry: fresh variable.
    !     on exit : this will contain the local sparse matrix.
    !
    !       interface parts
    !         !   .....user passed subroutine.....
    !         subroutine parts(global_indx,n,np,pv,nv)
    !           implicit none
    !           integer, intent(in)  :: global_indx, n, np
    !           integer, intent(out) :: nv
    !           integer, intent(out) :: pv(*)
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
    !  integer                                  :: ictxt
    !     on entry: blacs context.
    !     on exit : unchanged.
    !
    !  type (desc_type)                  :: desc_a
    !     on entry: fresh variable.
    !     on exit : the updated array descriptor
    !
    !  real(kind(1.d0)),  optional      :: b_glob(:)
    !     on entry: this contains right hand side.
    !     on exit :
    !
    !  real(kind(1.d0)), allocatable, optional      :: b(:)
    !     on entry: fresh variable.
    !     on exit : this will contain the local right hand side.
    !
    !  integer, optional    :: inroot
    !     on entry: specifies processor holding a_glob. default: 0
    !     on exit : unchanged.
    !
    use psb_base_mod
    implicit none

    ! parameters
    type(psb_zspmat_type)      :: a_glob
    complex(kind(1.d0))        :: b_glob(:)
    integer                    :: ictxt
    type(psb_zspmat_type)      :: a
    complex(kind(1.d0)), allocatable  :: b(:)
    type (psb_desc_type)       :: desc_a
    integer, intent(out)       :: info
    integer, optional          :: inroot
    character(len=5), optional :: fmt
    interface 

      !   .....user passed subroutine.....
      subroutine parts(global_indx,n,np,pv,nv)
        implicit none
        integer, intent(in)  :: global_indx, n, np
        integer, intent(out) :: nv
        integer, intent(out) :: pv(*) 
      end subroutine parts
    end interface

    ! local variables
    integer                     :: np, iam
    integer                     :: ircode, length_row, i_count, j_count,&
         & k_count, blockdim, root, liwork, nrow, ncol, nnzero, nrhs,&
         & i,j,k, ll, nz, isize, iproc, nnr, err, err_act, int_err(5)
    integer, allocatable            :: iwork(:)
    character                   :: afmt*5, atyp*5
    integer, allocatable          :: irow(:),icol(:)
    complex(kind(1.d0)), allocatable :: val(:)
    integer, parameter          :: nb=30
    real(kind(1.d0))            :: t0, t1, t2, t3, t4, t5
    character(len=20)           :: name, ch_err

    info = 0
    err  = 0
    name = 'mat_distf'
    call psb_erractionsave(err_act)

    ! executable statements    
    if (present(inroot)) then
      root = inroot
    else
      root = 0
    end if
    call psb_info(ictxt, iam, np)     
    if (iam == root) then
      ! extract information from a_glob
      if (a_glob%fida.ne. 'CSR') then
        info=135
        ch_err='CSR'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      endif
      nrow = a_glob%m
      ncol = a_glob%k
      if (nrow /= ncol) then
        write(0,*) 'a rectangular matrix ? ',nrow,ncol
        info=-1
        call psb_errpush(info,name)
        goto 9999
      endif
      nnzero = size(a_glob%aspk)
      nrhs   = 1
    endif
    ! broadcast informations to other processors
    call psb_bcast(ictxt,nrow, root)
    call psb_bcast(ictxt,ncol, root)
    call psb_bcast(ictxt,nnzero, root)
    call psb_bcast(ictxt,nrhs, root)
    liwork = max(np, nrow + ncol)
    allocate(iwork(liwork), stat = info)
    if (info /= 0) then
      info=2025
      int_err(1)=liwork
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif
    if (iam == root) then
      write (*, fmt = *) 'start matdist',root, size(iwork),&
           &nrow, ncol, nnzero,nrhs
    endif
    call psb_cdall(ictxt,desc_a,info,mg=nrow,parts=parts)
    if(info/=0) then
      info=4010
      ch_err='psb_cdall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    call psb_spall(a,desc_a,info,nnz=nnzero/np)
    if(info/=0) then
      info=4010
      ch_err='psb_psspall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    call psb_geall(b,desc_a,info)   
    if(info/=0) then
      info=4010
      ch_err='psb_psdsall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    isize = max(3*nb,ncol)


    allocate(val(nb*ncol),irow(nb*ncol),icol(nb*ncol),stat=info)
    if(info/=0) then
      info=4010
      ch_err='Allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    i_count   = 1

    do while (i_count.le.nrow)

      call parts(i_count,nrow,np,iwork, length_row)

      if (length_row.eq.1) then 
        j_count = i_count 
        iproc   = iwork(1) 
        do 
          j_count = j_count + 1 
          if (j_count-i_count >= nb) exit
          if (j_count > nrow) exit
          call parts(j_count,nrow,np,iwork, length_row)
          if (length_row /= 1 ) exit
          if (iwork(1) /= iproc ) exit
        end do

        ! now we should insert rows i_count..j_count-1
        nnr = j_count - i_count

        if (iam == root) then

          ll = 0
          do i= i_count, j_count-1
            call psb_sp_getrow(i,a_glob,nz,&
                 & irow(ll+1:),icol(ll+1:),val(ll+1:), info)
            if (info /= 0) then            
              if (nz >min(size(irow(ll+1:)),size(icol(ll+1:)),size(val(ll+1:)))) then 
                write(0,*) 'Allocation failure? This should not happen!'
              end if
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
            ll = ll + nz
          end do

          if (iproc == iam) then
            call psb_spins(ll,irow,icol,val,a,desc_a,info)
            if(info/=0) then
              info=4010
              ch_err='psb_spins'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
            call psb_geins(nnr,(/(i,i=i_count,j_count-1)/),b_glob(i_count:j_count-1),&
                 & b,desc_a,info)
            if(info/=0) then
              info=4010
              ch_err='psb_ins'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
          else
            call psb_snd(ictxt,nnr,iproc)
            call psb_snd(ictxt,ll,iproc)
            call psb_snd(ictxt,irow(1:ll),iproc)
            call psb_snd(ictxt,icol(1:ll),iproc)
            call psb_snd(ictxt,val(1:ll),iproc)
            call psb_snd(ictxt,b_glob(i_count:j_count-1),iproc)
            call psb_rcv(ictxt,ll,iproc)
          endif
        else if (iam /= root) then

          if (iproc == iam) then
            call psb_rcv(ictxt,nnr,root)
            call psb_rcv(ictxt,ll,root)
            if (ll > size(irow)) then 
              write(0,*) iam,'need to reallocate ',ll
              deallocate(val,irow,icol)
              allocate(val(ll),irow(ll),icol(ll),stat=info)
              if(info/=0) then
                info=4010
                ch_err='Allocate'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if

            endif
            call psb_rcv(ictxt,irow(1:ll),root)
            call psb_rcv(ictxt,icol(1:ll),root)
            call psb_rcv(ictxt,val(1:ll),root)
            call psb_rcv(ictxt,b_glob(i_count:i_count+nnr-1),root)
            call psb_snd(ictxt,ll,root)
            call psb_spins(ll,irow,icol,val,a,desc_a,info)
            if(info/=0) then
              info=4010
              ch_err='psspins'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
            call psb_geins(nnr,(/(i,i=i_count,i_count+nnr-1)/),&
                 & b_glob(i_count:i_count+nnr-1),b,desc_a,info)
            if(info/=0) then
              info=4010
              ch_err='psdsins'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
          endif
        endif

        i_count = j_count

      else
        write(0,*) iam,'unexpected turn'
        ! here processors are counted 1..np
        do j_count = 1, length_row
          k_count = iwork(j_count)
          if (iam == root) then

            ll = 0
            do i= i_count, i_count
              call psb_sp_getrow(i,a_glob,nz,&
                   & irow(ll+1:),icol(ll+1:),val(ll+1:), info)
              if (info /= 0) then            
                if (nz >min(size(irow(ll+1:)),size(icol(ll+1:)),size(val(ll+1:)))) then 
                  write(0,*) 'Allocation failure? This should not happen!'
                end if
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
              ll = ll + nz
            end do

            if (k_count == iam) then

              call psb_spins(ll,irow,icol,val,a,desc_a,info)
              if(info/=0) then
                info=4010
                ch_err='psspins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
              call psb_geins(1,(/i_count/),b_glob(i_count:i_count),&
                   & b,desc_a,info)
              if(info/=0) then
                info=4010
                ch_err='psdsins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
            else
              call psb_snd(ictxt,ll,k_count)
              call psb_snd(ictxt,irow(1:ll),k_count)
              call psb_snd(ictxt,icol(1:ll),k_count)
              call psb_snd(ictxt,val(1:ll),k_count)
              call psb_snd(ictxt,b_glob(i_count),k_count)
              call psb_rcv(ictxt,ll,k_count)
            endif
          else if (iam /= root) then
            if (k_count == iam) then
              call psb_rcv(ictxt,ll,root)
              call psb_rcv(ictxt,irow(1:ll),root)
              call psb_rcv(ictxt,icol(1:ll),root)
              call psb_rcv(ictxt,val(1:ll),root)
              call psb_rcv(ictxt,b_glob(i_count),root)
              call psb_snd(ictxt,ll,root)
              call psb_spins(ll,irow,icol,val,a,desc_a,info)
              if(info/=0) then
                info=4010
                ch_err='psspins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
              call psb_geins(1,(/i_count/),b_glob(i_count:i_count),&
                   & b,desc_a,info)
              if(info/=0) then
                info=4010
                ch_err='psdsins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
            endif
          endif
        end do
        i_count = i_count + 1
      endif
    end do

    if (present(fmt)) then 
      afmt=fmt
    else
      afmt = 'CSR'
    endif

    call psb_barrier(ictxt)
    t0 = psb_wtime()
    call psb_cdasb(desc_a,info)     
    t1 = psb_wtime()
    if(info/=0)then
      info=4010
      ch_err='psb_cdasb'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_barrier(ictxt)
    t2 = psb_wtime()
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)     
    t3 = psb_wtime()
    if(info/=0)then
      info=4010
      ch_err='psb_spasb'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if


    if (iam == root) then 
      write(*,*) 'descriptor assembly: ',t1-t0
      write(*,*) 'sparse matrix assembly: ',t3-t2
    end if

    call psb_geasb(b,desc_a,info)     
    if(info/=0)then
      info=4010
      ch_err='psdsasb'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    deallocate(val,irow,icol,stat=info)
    if(info/=0)then
      info=4010
      ch_err='deallocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    deallocate(iwork)   
    if (iam == root) write (*, fmt = *) 'end matdist'     

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return

  end subroutine zmatdistf


  subroutine zmatdistv (a_glob, a, v, ictxt, desc_a,&
       & b_glob, b, info, inroot,fmt)
    !
    ! an utility subroutine to distribute a matrix among processors
    ! according to a user defined data distribution, using 
    ! sparse matrix subroutines.
    !
    !  type(d_spmat)                            :: a_glob
    !     on entry: this contains the global sparse matrix as follows:
    !        a%fida =='csr'
    !        a%aspk for coefficient values
    !        a%ia1  for column indices
    !        a%ia2  for row pointers
    !        a%m    for number of global matrix rows
    !        a%k    for number of global matrix columns
    !     on exit : undefined, with unassociated pointers.
    !
    !  type(d_spmat)                            :: a
    !     on entry: fresh variable.
    !     on exit : this will contain the local sparse matrix.
    !
    !       interface parts
    !         !   .....user passed subroutine.....
    !         subroutine parts(global_indx,n,np,pv,nv)
    !           implicit none
    !           integer, intent(in)  :: global_indx, n, np
    !           integer, intent(out) :: nv
    !           integer, intent(out) :: pv(*)
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
    !  integer                                  :: ictxt
    !     on entry: blacs context.
    !     on exit : unchanged.
    !
    !  type (desc_type)                  :: desc_a
    !     on entry: fresh variable.
    !     on exit : the updated array descriptor
    !
    !  real(kind(1.d0)),  optional      :: b_glob(:)
    !     on entry: this contains right hand side.
    !     on exit :
    !
    !  real(kind(1.d0)), allocatable, optional      :: b(:)
    !     on entry: fresh variable.
    !     on exit : this will contain the local right hand side.
    !
    !  integer, optional    :: inroot
    !     on entry: specifies processor holding a_glob. default: 0
    !     on exit : unchanged.
    !
    use psb_base_mod
    implicit none   ! parameters
    type(psb_zspmat_type)      :: a_glob
    complex(kind(1.d0))        :: b_glob(:)
    integer                    :: ictxt, v(:)
    type(psb_zspmat_type)      :: a
    complex(kind(1.d0)), allocatable  :: b(:)
    type(psb_desc_type)       :: desc_a
    integer, intent(out)       :: info
    integer, optional          :: inroot
    character(len=5), optional :: fmt

    integer                     :: np, iam
    integer                     :: ircode, length_row, i_count, j_count,&
         & k_count, blockdim, root, liwork, nrow, ncol, nnzero, nrhs,&
         & i,j,k, ll, nz, isize, iproc, nnr, err, err_act, int_err(5)
    integer, allocatable            :: iwork(:)
    character                   :: afmt*5, atyp*5
    integer, allocatable          :: irow(:),icol(:)
    complex(kind(1.d0)), allocatable :: val(:)
    integer, parameter          :: nb=30
    logical, parameter          :: newt=.true.
    real(kind(1.d0))            :: t0, t1, t2, t3, t4, t5
    character(len=20)  :: name, ch_err

    info = 0
    err  = 0
    name = 'mat_distv'
    call psb_erractionsave(err_act)

    ! executable statements    
    if (present(inroot)) then
      root = inroot
    else
      root = 0
    end if

    call psb_info(ictxt, iam, np)     
    if (iam == root) then
      ! extract information from a_glob
      if (toupper(a_glob%fida) /=  'CSR') then
        info=135
        ch_err='CSR'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      endif

      nrow = a_glob%m
      ncol = a_glob%k
      if (nrow /= ncol) then
        write(0,*) 'a rectangular matrix ? ',nrow,ncol
        info=-1
        call psb_errpush(info,name)
        goto 9999
      endif

      nnzero = size(a_glob%aspk)
      nrhs   = 1
    end if
    ! broadcast informations to other processors
    call psb_bcast(ictxt,nrow, root)
    call psb_bcast(ictxt,ncol, root)
    call psb_bcast(ictxt,nnzero, root)
    call psb_bcast(ictxt,nrhs, root)
    liwork = max(np, nrow + ncol)
    allocate(iwork(liwork), stat = info)
    if (info /= 0) then
      write(0,*) 'matdist allocation failed'
      info=2025
      int_err(1)=liwork
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif

    call psb_cdall(ictxt,desc_a,info,vg=v)
    if(info/=0) then
      info=4010
      ch_err='psb_cdall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_spall(a,desc_a,info,nnz=((nnzero+np-1)/np))
    if(info/=0) then
      info=4010
      ch_err='psb_psspall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    call psb_geall(b,desc_a,info)   
    if(info/=0) then
      info=4010
      ch_err='psb_psdsall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    isize = max(3*nb,ncol)


    allocate(val(nb*ncol),irow(nb*ncol),icol(nb*ncol),stat=info)
    if(info/=0) then
      info=4010
      ch_err='Allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    i_count   = 1

    do while (i_count <= nrow)

      j_count = i_count 
      iproc   = v(i_count)

      do 
        j_count = j_count + 1 
        if (j_count-i_count >= nb) exit
        if (j_count > nrow) exit
        if (v(j_count) /= iproc ) exit
      end do

      ! now we should insert rows i_count..j_count-1
      nnr = j_count - i_count

      if (iam == root) then
        ll = 0
        do i= i_count, j_count-1
          call psb_sp_getrow(i,a_glob,nz,&
               & irow(ll+1:),icol(ll+1:),val(ll+1:), info)
          if (info /= 0) then            
            if (nz >min(size(irow(ll+1:)),size(icol(ll+1:)),size(val(ll+1:)))) then 
              write(0,*) 'Allocation failure? This should not happen!'
            end if
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
          ll = ll + nz
        end do

        if (iproc == iam) then
          call psb_spins(ll,irow,icol,val,a,desc_a,info)
          if(info/=0) then
            info=4010
            ch_err='psb_spins'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if

          call psb_geins(nnr,(/(i,i=i_count,j_count-1)/),b_glob(i_count:j_count-1),&
               & b,desc_a,info)
          if(info/=0) then
            info=4010
            ch_err='dsins'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        else
          call psb_snd(ictxt,nnr,iproc)
          call psb_snd(ictxt,ll,iproc)
          call psb_snd(ictxt,irow(1:ll),iproc)
          call psb_snd(ictxt,icol(1:ll),iproc)
          call psb_snd(ictxt,val(1:ll),iproc)
          call psb_snd(ictxt,b_glob(i_count:j_count-1),iproc)
          call psb_rcv(ictxt,ll,iproc)
        endif
      else if (iam /= root) then

        if (iproc == iam) then
          call psb_rcv(ictxt,nnr,root)
          call psb_rcv(ictxt,ll,root)
          if (ll > size(val)) then 
            write(0,*) iam,'need to reallocate ',ll
            deallocate(val,irow,icol)
            allocate(val(ll),irow(ll),icol(ll),stat=info)
            if(info/=0) then
              info=4010
              ch_err='Allocate'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
          endif
          call psb_rcv(ictxt,irow(1:ll),root)
          call psb_rcv(ictxt,icol(1:ll),root)
          call psb_rcv(ictxt,val(1:ll),root)
          call psb_rcv(ictxt,b_glob(i_count:i_count+nnr-1),root)
          call psb_snd(ictxt,ll,root)

          call psb_spins(ll,irow,icol,val,a,desc_a,info)
          if(info/=0) then
            info=4010
            ch_err='spins'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
          call psb_geins(nnr,(/(i,i=i_count,i_count+nnr-1)/),&
               & b_glob(i_count:i_count+nnr-1),b,desc_a,info)
          if(info/=0) then
            info=4010
            ch_err='psdsins'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        endif
      endif
      i_count = j_count

    end do

    ! default storage format for sparse matrix; we do not
    ! expect duplicated entries.

    if (present(fmt)) then  
      afmt=fmt
    else
      afmt = 'CSR'
    endif
    call psb_barrier(ictxt)
    t0 = psb_wtime()
    call psb_cdasb(desc_a,info)     
    t1 = psb_wtime()
    if(info/=0)then
      info=4010
      ch_err='psb_cdasb'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_barrier(ictxt)
    t2 = psb_wtime()
    call psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)     
    t3 = psb_wtime()
    if(info/=0)then
      info=4010
      ch_err='psb_spasb'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_geasb(b,desc_a,info)

    if (iam == root) then 
      write(*,'("Descriptor assembly   : ",es10.4)')t1-t0
      write(*,'("Sparse matrix assembly: ",es10.4)')t3-t2
    end if

    if(info/=0)then
      info=4010
      ch_err='psdsasb'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    deallocate(iwork)   

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error(ictxt)
      return
    end if
    return

  end subroutine zmatdistv

end module psb_mat_dist_mod
