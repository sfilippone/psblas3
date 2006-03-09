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
module mat_dist

  interface matdist
     module procedure dmatdistf, dmatdistv
  end interface

contains

  subroutine dmatdistf (a_glob, a, parts, icontxt, desc_a,&
       & b_glob, b, info, inroot,fmt)
    !
    ! an utility subroutine to distribute a matrix among processors
    ! according to a user defined data distribution, using pessl
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
    !  integer                                  :: icontxt
    !     on entry: blacs context.
    !     on exit : unchanged.
    !
    !  type (desc_type)                  :: desc_a
    !     on entry: fresh variable.
    !     on exit : the updated array descriptor
    !
    !  real(kind(1.d0)), pointer, optional      :: b_glob(:)
    !     on entry: this contains right hand side.
    !     on exit :
    !
    !  real(kind(1.d0)), pointer, optional      :: b(:)
    !     on entry: fresh variable.
    !     on exit : this will contain the local right hand side.
    !
    !  integer, optional    :: inroot
    !     on entry: specifies processor holding a_glob. default: 0
    !     on exit : unchanged.
    !
    use psb_sparse_mod
    implicit none

    ! parameters
    type(psb_dspmat_type)      :: a_glob
    real(kind(1.d0)), pointer  :: b_glob(:)
    integer                    :: icontxt
    type(psb_dspmat_type)      :: a
    real(kind(1.d0)), pointer  :: b(:)
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
    integer                     :: nprow, npcol, myprow, mypcol
    integer                     :: ircode, length_row, i_count, j_count,&
         & k_count, blockdim, root, liwork, nrow, ncol, nnzero, nrhs,&
         & i,j,k, ll, isize, iproc, nnr, err, err_act, int_err(5)
    integer, pointer            :: iwork(:)
    character                   :: afmt*5, atyp*5
    type(psb_dspmat_type)       :: blck   
    integer, parameter          :: nb=30
    real(kind(1.d0))            :: t0, t1, t2, t3, t4, t5, mpi_wtime
    external                    :: mpi_wtime
    logical, parameter          :: newt=.true.
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
    call blacs_gridinfo(icontxt, nprow, npcol, myprow, mypcol)     
    if (myprow == root) then
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
      call gebs2d(icontxt, 'a', nrow)
      call gebs2d(icontxt, 'a', ncol)
      call gebs2d(icontxt, 'a', nnzero)
      call gebs2d(icontxt, 'a', nrhs)
    else !(myprow /= root)
      ! receive informations
      call gebr2d(icontxt, 'a', nrow)
      call gebr2d(icontxt, 'a', ncol)
      call gebr2d(icontxt, 'a', nnzero)
      call gebr2d(icontxt, 'a', nrhs)
    end if   ! allocate integer work area
    liwork = max(nprow, nrow + ncol)
    allocate(iwork(liwork), stat = info)
    if (info /= 0) then
      info=2025
      int_err(1)=liwork
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif
    if (myprow == root) then
      write (*, fmt = *) 'start matdist',root, size(iwork),&
           &nrow, ncol, nnzero,nrhs
    endif
    if (newt) then 
      call psb_cdall(nrow,nrow,parts,icontxt,desc_a,info)
      if(info/=0) then
        info=4010
        ch_err='psb_cdall'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    else
      call psb_cdall(nrow,nrow,parts,icontxt,desc_a,info)
      if(info/=0) then
        info=4010
        ch_err='psb_psdscall'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    endif
    call psb_spall(a,desc_a,info,nnz=nnzero/nprow)
    if(info/=0) then
      info=4010
      ch_err='psb_psspall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    call psb_geall(nrow,b,desc_a,info)   
    if(info/=0) then
      info=4010
      ch_err='psb_psdsall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    isize = max(3*nb,ncol)


    blck%m = nb
    blck%k = ncol
    call psb_sp_all(blck,nb*ncol,info)
    if(info/=0) then
      info=4010
      ch_err='spall'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    blck%fida = 'CSR'
    i_count   = 1

    do while (i_count.le.nrow)

      call parts(i_count,nrow,nprow,iwork, length_row)

      if (length_row.eq.1) then 
        j_count = i_count 
        iproc   = iwork(1) 
        do 
          j_count = j_count + 1 
          if (j_count-i_count >= nb) exit
          if (j_count > nrow) exit
          call parts(j_count,nrow,nprow,iwork, length_row)
          if (length_row /= 1 ) exit
          if (iwork(1) /= iproc ) exit
        end do

        ! now we should insert rows i_count..j_count-1
        nnr = j_count - i_count

        if (myprow == root) then

          do j = i_count, j_count
            blck%ia2(j-i_count+1) = a_glob%ia2(j) - &
                 & a_glob%ia2(i_count) + 1
          enddo

          k = a_glob%ia2(i_count)
          do j = k, a_glob%ia2(j_count)-1
            blck%aspk(j-k+1) = a_glob%aspk(j)
            blck%ia1(j-k+1) = a_glob%ia1(j)
          enddo

          ll     = blck%ia2(nnr+1) - 1
          blck%m = nnr
          blck%k = nrow
          if (iproc == myprow) then
            call psb_spins(ll,blck%ia1,blck%ia2,blck%aspk,a,desc_a,info)
            if(info/=0) then
              info=4010
              ch_err='psb_spins'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
            call psb_geins(nnr,b,i_count,b_glob(i_count:j_count-1),&
                 &desc_a,info)
            if(info/=0) then
              info=4010
              ch_err='psb_ins'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
          else
            call igesd2d(icontxt,1,1,nnr,1,iproc,0)
            call igesd2d(icontxt,1,1,ll,1,iproc,0)
            call igesd2d(icontxt,nnr+1,1,blck%ia2,nnr+1,iproc,0)
            call igesd2d(icontxt,ll,1,blck%ia1,ll,iproc,0)
            call dgesd2d(icontxt,ll,1,blck%aspk,ll,iproc,0)
            call dgesd2d(icontxt,nnr,1,b_glob(i_count:j_count-1),nnr,iproc,0)
            call igerv2d(icontxt,1,1,ll,1,iproc,0)
          endif
        else if (myprow /= root) then

          if (iproc == myprow) then
            call igerv2d(icontxt,1,1,nnr,1,root,0)
            call igerv2d(icontxt,1,1,ll,1,root,0)
            if (ll > size(blck%ia1)) then 
              write(0,*) myprow,'need to reallocate ',ll
              call psb_sp_reall(blck,ll,info)
              if(info/=0) then
                info=4010
                ch_err='psb_sp_reall'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if

            endif
            call igerv2d(icontxt,ll,1,blck%ia1,ll,root,0)
            call igerv2d(icontxt,nnr+1,1,blck%ia2,nnr+1,root,0)
            call dgerv2d(icontxt,ll,1,blck%aspk,ll,root,0)
            call dgerv2d(icontxt,nnr,1,b_glob(i_count:i_count+nnr-1),nnr,root,0)
            call igesd2d(icontxt,1,1,ll,1,root,0)
            blck%m = nnr
            blck%k = nrow
            blck%infoa(psb_nnz_) = ll
            call psb_spins(ll,blck%ia1,blck%ia2,blck%aspk,a,desc_a,info)
            if(info/=0) then
              info=4010
              ch_err='psspins'
              call psb_errpush(info,name,a_err=ch_err)
              goto 9999
            end if
            call psb_geins(nnr,b,i_count,b_glob(i_count:i_count+nnr-1),&
                 &desc_a,info)
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
        write(0,*) myprow,'unexpected turn'
        ! here processors are counted 1..nprow
        do j_count = 1, length_row
          k_count = iwork(j_count)
          if (myprow == root) then
            blck%ia2(1) = 1
            blck%ia2(2) = 1
            do j = a_glob%ia2(i_count), a_glob%ia2(i_count+1)-1
              blck%aspk(blck%ia2(2)) = a_glob%aspk(j)
              blck%ia1(blck%ia2(2)) = a_glob%ia1(j)
              blck%ia2(2) =blck%ia2(2) + 1
            enddo
            ll = blck%ia2(2) - 1
            if (k_count == myprow) then
              blck%infoa(1) = ll
              blck%infoa(2) = ll
              blck%infoa(3) = 2
              blck%infoa(4) = 1
              blck%infoa(5) = 1
              blck%infoa(6) = 1
              blck%m = 1
              blck%k = nrow

              call psb_spins(ll,blck%ia1,blck%ia2,blck%aspk,a,desc_a,info)
              if(info/=0) then
                info=4010
                ch_err='psspins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
              call psb_geins(1,b,i_count,b_glob(i_count:i_count),&
                   &desc_a,info)
              if(info/=0) then
                info=4010
                ch_err='psdsins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
            else
              call igesd2d(icontxt,1,1,ll,1,k_count,0)
              call igesd2d(icontxt,ll,1,blck%ia1,ll,k_count,0)
              call dgesd2d(icontxt,ll,1,blck%aspk,ll,k_count,0)
              call dgesd2d(icontxt,1,1,b_glob(i_count),1,k_count,0)
              call igerv2d(icontxt,1,1,ll,1,k_count,0)
            endif
          else if (myprow /= root) then
            if (k_count == myprow) then
              call igerv2d(icontxt,1,1,ll,1,root,0)
              blck%ia2(1) = 1
              blck%ia2(2) = ll+1
              call igerv2d(icontxt,ll,1,blck%ia1,ll,root,0)
              call dgerv2d(icontxt,ll,1,blck%aspk,ll,root,0)
              call dgerv2d(icontxt,1,1,b_glob(i_count),1,root,0)
              call igesd2d(icontxt,1,1,ll,1,root,0)
              blck%m = 1
              blck%k = nrow          
              call psb_spins(ll,blck%ia1,blck%ia2,blck%aspk,a,desc_a,info)
              if(info/=0) then
                info=4010
                ch_err='psspins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
              end if
              call psb_geins(1,b,i_count,b_glob(i_count:i_count),&
                   &desc_a,info)
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
    if (newt) then 

      call blacs_barrier(icontxt,'all')
      t0 = mpi_wtime()
      call psb_cdasb(desc_a,info)     
      t1 = mpi_wtime()
      if(info/=0)then
        info=4010
        ch_err='psb_cdasb'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

      call blacs_barrier(icontxt,'all')
      t2 = mpi_wtime()
      call psb_spasb(a,desc_a,info,dup=1,afmt=afmt)     
      t3 = mpi_wtime()
      if(info/=0)then
        info=4010
        ch_err='psb_spasb'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if


      if (myprow == root) then 
        write(*,*) 'descriptor assembly: ',t1-t0
        write(*,*) 'sparse matrix assembly: ',t3-t2
      end if


    else
      call psb_spasb(a,desc_a,info,afmt=afmt,dup=1)     
      if(info/=0)then
        info=4010
        ch_err='psspasb'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    endif


    call psb_geasb(b,desc_a,info)     
    if(info/=0)then
      info=4010
      ch_err='psdsasb'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    call psb_spfree(blck,info)
    if(info/=0)then
      info=4010
      ch_err='spfree'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    deallocate(iwork)   
    if (myprow == root) write (*, fmt = *) 'end matdist'     

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.act_abort) then
      call psb_error(icontxt)
      return
    end if
    return

  end subroutine dmatdistf


  subroutine dmatdistv (a_glob, a, v, icontxt, desc_a,&
       & b_glob, b, info, inroot,fmt)
    !
    ! an utility subroutine to distribute a matrix among processors
    ! according to a user defined data distribution, using pessl
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
    !  integer                                  :: icontxt
    !     on entry: blacs context.
    !     on exit : unchanged.
    !
    !  type (desc_type)                  :: desc_a
    !     on entry: fresh variable.
    !     on exit : the updated array descriptor
    !
    !  real(kind(1.d0)), pointer, optional      :: b_glob(:)
    !     on entry: this contains right hand side.
    !     on exit :
    !
    !  real(kind(1.d0)), pointer, optional      :: b(:)
    !     on entry: fresh variable.
    !     on exit : this will contain the local right hand side.
    !
    !  integer, optional    :: inroot
    !     on entry: specifies processor holding a_glob. default: 0
    !     on exit : unchanged.
    !
    use psb_sparse_mod
    implicit none   ! parameters
    type(psb_dspmat_type)      :: a_glob
    real(kind(1.d0)), pointer  :: b_glob(:)
    integer                    :: icontxt, v(:)
    type(psb_dspmat_type)      :: a
    real(kind(1.d0)), pointer  :: b(:)
    type (psb_desc_type)       :: desc_a
    integer, intent(out)       :: info
    integer, optional          :: inroot
    character(len=5), optional :: fmt

    integer                     :: nprow, npcol, myprow, mypcol
    integer                     :: ircode, length_row, i_count, j_count,&
         & k_count, blockdim, root, liwork, nrow, ncol, nnzero, nrhs,&
         & i,j,k, ll, isize, iproc, nnr, err, err_act, int_err(5)
    integer, pointer            :: iwork(:)
    character                   :: afmt*5, atyp*5
    type(psb_dspmat_type)       :: blck   
    integer, parameter          :: nb=30
    logical, parameter          :: newt=.true.
    real(kind(1.d0))            :: t0, t1, t2, t3, t4, t5, mpi_wtime
    external                    :: mpi_wtime
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

    call blacs_gridinfo(icontxt, nprow, npcol, myprow, mypcol)     
    if (myprow == root) then
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
       call igebs2d(icontxt, 'a', ' ', 1, 1, nrow, 1)
       call igebs2d(icontxt, 'a', ' ', 1, 1, ncol, 1)
       call igebs2d(icontxt, 'a', ' ', 1, 1, nnzero, 1)
       call igebs2d(icontxt, 'a', ' ', 1, 1, nrhs, 1)
    else !(myprow /= root)
       ! receive informations
       call igebr2d(icontxt, 'a', ' ', 1, 1, nrow, 1, root, 0)
       call igebr2d(icontxt, 'a', ' ', 1, 1, ncol, 1, root, 0)
       call igebr2d(icontxt, 'a', ' ', 1, 1, nnzero, 1, root, 0)
       call igebr2d(icontxt, 'a', ' ', 1, 1, nrhs, 1, root, 0)
    end if   ! allocate integer work area
    liwork = max(nprow, nrow + ncol)
    allocate(iwork(liwork), stat = info)
    if (info /= 0) then
       write(0,*) 'matdist allocation failed'
       info=2025
       int_err(1)=liwork
       call psb_errpush(info,name,i_err=int_err)
       goto 9999
    endif

    call psb_cdall(nrow,v,icontxt,desc_a,info)
    if(info/=0) then
       info=4010
       ch_err='psb_cdall'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if

    call psb_spall(a,desc_a,info,nnz=((nnzero+nprow-1)/nprow))
    if(info/=0) then
       info=4010
       ch_err='psb_psspall'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if
    call psb_geall(nrow,b,desc_a,info)   
    if(info/=0) then
       info=4010
       ch_err='psb_psdsall'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if
    isize = max(3*nb,ncol)


    blck%m = nb
    blck%k = ncol
    call psb_sp_all(blck,nb*ncol,info)
    if(info/=0) then
       info=4010
       ch_err='spall'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if
    blck%fida = 'COO'

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

       if (myprow == root) then
          ll = a_glob%ia2(j_count)-a_glob%ia2(i_count)
          if (ll > size(blck%aspk)) then 
             call psb_sp_reall(blck,ll,info)
             if(info/=0) then
                info=4010
                ch_err='spreall'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
             end if

          endif
          k = a_glob%ia2(i_count)
          do i= i_count, j_count-1
             do j = a_glob%ia2(i),a_glob%ia2(i+1)-1
                blck%ia1(j-k+1)  = i
                blck%ia2(j-k+1)  = a_glob%ia1(j) 
                blck%aspk(j-k+1) = a_glob%aspk(j)
             end do
          enddo

          blck%m = nnr
          blck%k = nrow
          blck%infoa(psb_nnz_) = ll
          if (iproc == myprow) then
             call psb_spins(ll,blck%ia1,blck%ia2,blck%aspk,a,desc_a,info)
             if(info/=0) then
                info=4010
                ch_err='psb_spins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
             end if

             call psb_geins(nnr,b,i_count,b_glob(i_count:j_count-1),&
                  &desc_a,info)
             if(info/=0) then
                info=4010
                ch_err='dsins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
             end if
          else
             call igesd2d(icontxt,1,1,nnr,1,iproc,0)
             call igesd2d(icontxt,1,1,ll,1,iproc,0)
             call igesd2d(icontxt,ll,1,blck%ia1,ll,iproc,0)
             call igesd2d(icontxt,ll,1,blck%ia2,ll,iproc,0)
             call dgesd2d(icontxt,ll,1,blck%aspk,ll,iproc,0)
             call dgesd2d(icontxt,nnr,1,b_glob(i_count:j_count-1),nnr,iproc,0)
             call igerv2d(icontxt,1,1,ll,1,iproc,0)
          endif
       else if (myprow /= root) then

          if (iproc == myprow) then
             call igerv2d(icontxt,1,1,nnr,1,root,0)
             call igerv2d(icontxt,1,1,ll,1,root,0)
             if (ll > size(blck%aspk)) then 
                write(0,*) myprow,'need to reallocate ',ll
                call psb_sp_reall(blck,ll,info)
                if(info/=0) then
                   info=4010
                   ch_err='spreall'
                   call psb_errpush(info,name,a_err=ch_err)
                   goto 9999
                end if
             endif
             call igerv2d(icontxt,ll,1,blck%ia1,ll,root,0)
             call igerv2d(icontxt,ll,1,blck%ia2,ll,root,0)
             call dgerv2d(icontxt,ll,1,blck%aspk,ll,root,0)
             call dgerv2d(icontxt,nnr,1,b_glob(i_count:i_count+nnr-1),nnr,root,0)
             call igesd2d(icontxt,1,1,ll,1,root,0)
             blck%m = nnr
             blck%k = nrow
             blck%infoa(psb_nnz_) = ll
             call psb_spins(ll,blck%ia1,blck%ia2,blck%aspk,a,desc_a,info)
             if(info/=0) then
                info=4010
                ch_err='spins'
                call psb_errpush(info,name,a_err=ch_err)
                goto 9999
             end if
             call psb_ins(nnr,b,i_count,b_glob(i_count:i_count+nnr-1),&
                  &desc_a,info)
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
    call blacs_barrier(icontxt,'all')
    t0 = mpi_wtime()
    call psb_cdasb(desc_a,info)     
    t1 = mpi_wtime()
    if(info/=0)then
       info=4010
       ch_err='psb_cdasb'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if

    call blacs_barrier(icontxt,'all')
    t2 = mpi_wtime()
    call psb_spasb(a,desc_a,info,dup=1,afmt=afmt)     
    t3 = mpi_wtime()
    if(info/=0)then
       info=4010
       ch_err='psb_spasb'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if

    call psb_geasb(b,desc_a,info)

    if (myprow == root) then 
       write(*,'("Descriptor assembly   : ",es10.4)')t1-t0
       write(*,'("Sparse matrix assembly: ",es10.4)')t3-t2
    end if

    if(info/=0)then
       info=4010
       ch_err='psdsasb'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if
    call psb_spfree(blck,info)
    if(info/=0)then
       info=4010
       ch_err='spfree'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if

    deallocate(iwork)   

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.act_abort) then
       call psb_error(icontxt)
       return
    end if
    return

  end subroutine dmatdistv

end module mat_dist
