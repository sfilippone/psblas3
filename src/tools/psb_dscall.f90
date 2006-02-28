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
! File: psb_dscall.f90
!
! Subroutine: psb_dscall
!    Allocate descriptor
!    and checks correctness of PARTS subroutine
! 
! Parameters: 
!    m       - integer.                       The number of rows.
!    n       - integer.                       The number of columns.
!    parts   - external subroutine.           The routine that contains the partitioning scheme.
!    icontxt - integer.                       The communication context.
!    desc_a  - type(<psb_desc_type>).         The communication descriptor.
!    info    - integer.                       Eventually returns an error code
subroutine psb_dscall(m, n, parts, icontxt, desc_a, info)
  use psb_error_mod
  use psb_descriptor_type
  use psb_realloc_mod
  use psb_serial_mod
  use psb_const_mod
  implicit None
  include 'parts.fh'
  !....Parameters...
  Integer, intent(in)                 :: M,N,ICONTXT
  Type(psb_desc_type), intent(out)    :: desc_a
  integer, intent(out)                :: info

  !locals
  Integer             :: counter,i,j,nprow,npcol,myrow,mycol,&
       & loc_row,err,loc_col,nprocs,&
       & l_ov_ix,l_ov_el,idx, err_act, itmpov, k
  Integer             :: INT_ERR(5),TEMP(1),EXCH(2)
  Real(Kind(1.d0))    :: REAL_ERR(5)
  Integer, Pointer    :: PRC_V(:), TEMP_OVRLAP(:), OV_IDX(:),OV_EL(:)
  logical, parameter  :: debug=.false.
  character(len=20)   :: name, char_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  err=0
  name = 'psb_dscall'
  call psb_erractionsave(err_act)

  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
  if (debug) write(*,*) 'psb_dscall: ',nprow,npcol,myrow,mycol
  !     ....verify blacs grid correctness..
  if (npcol /= 1) then
     info = 2030
     err=info
     int_err(1) = npcol
     call psb_errpush(err,name,int_err)
     goto 9999
  endif


  !... check m and n parameters....
  if (m.lt.1) then
     info = 10
     err=info
     int_err(1) = 1
     int_err(2) = m
     call psb_errpush(err,name,int_err)
     goto 9999
  else if (n.lt.1) then
     info = 10
     err=info
     int_err(1) = 2
     int_err(2) = n
     call psb_errpush(err,name,int_err)
     goto 9999
  endif


  if (debug) write(*,*) 'psb_dscall:  doing global checks'  
  !global check on m and n parameters
  if (myrow.eq.psb_root_) then
     exch(1)=m
     exch(2)=n
     call igebs2d(icontxt,psb_all_,psb_topdef_, itwo,ione, exch, itwo)
  else
     call igebr2d(icontxt,psb_all_,psb_topdef_, itwo,ione, exch, itwo, psb_root_,&
          & 0)
     if (exch(1) /= m) then
        err=550
        int_err(1)=1
        call psb_errpush(err,name,int_err)
        goto 9999
     else if (exch(2) /= n) then
        err=550
        int_err(1)=2
        call psb_errpush(err,name,int_err)
        goto 9999
     endif
  endif

  call psb_nullify_desc(desc_a)

  !count local rows number
  ! allocate work vector
  allocate(prc_v(nprow),desc_a%glob_to_loc(m),&
       &desc_a%matrix_data(psb_mdata_size_),temp_ovrlap(m),stat=info)
  if (info /= no_err) then     
     info=2025
     err=info
     int_err(1)=m
     call psb_errpush(err,name,int_err)
     goto 9999
  endif


  if (debug) write(*,*) 'PSB_DSCALL:  starting main loop' ,info
  counter = 0
  itmpov  = 0
  temp_ovrlap(:) = -1
  do i=1,m
     if (info.eq.0) then
        call parts(i,m,nprow,prc_v,nprocs)
        if (nprocs.gt.nprow) then
           info=570
           int_err(1)=3
           int_err(2)=nprow
           int_err(3)=nprocs
           int_err(4)=i
           err=info
           call psb_errpush(err,name,int_err)
           goto 9999
        else if (nprocs.le.0) then
           info=575
           int_err(1)=3
           int_err(2)=nprocs
           int_err(3)=i
           err=info
           call psb_errpush(err,name,int_err)
           goto 9999
        else
           do j=1,nprocs
              if ((prc_v(j).gt.nprow-1).or.(prc_v(j).lt.0)) then
                 info=580
                 int_err(1)=3
                 int_err(2)=prc_v(j)
                 int_err(3)=i
                 err=info
                 call psb_errpush(err,name,int_err)
                 goto 9999
              end if
           end do
        endif
        desc_a%glob_to_loc(i) = -(nprow+prc_v(1)+1)
        j=1
!!$      do while ((j.le.nprocs).and.(prc_v(j).ne.myrow))
        do 
           if (j > nprocs) exit
           if (prc_v(j) == myrow) exit
           j=j+1
        enddo
        if (j.le.nprocs) then 
           if (prc_v(j).eq.myrow) then
              ! this point belongs to myrow
              counter=counter+1
              desc_a%glob_to_loc(i) = counter
              if (nprocs.gt.1)  then
                 if ((itmpov+2+nprocs).gt.m)  then
                    info=2025
                    int_err(1)=m
                    err=info
                    call psb_errpush(err,name,int_err)
                    goto 9999
                 else
                    itmpov = itmpov + 1
                    temp_ovrlap(itmpov) = i
                    itmpov = itmpov + 1
                    temp_ovrlap(itmpov) = nprocs
                    temp_ovrlap(itmpov+1:itmpov+nprocs) = prc_v(1:nprocs)
                    itmpov = itmpov + nprocs
                 endif
              endif
           end if
        end if
     endif
  enddo

  loc_row=counter
  ! check on parts function
  if (debug) write(*,*) 'PSB_DSCALL:  End main loop:' ,loc_row,itmpov,info


  if (debug) write(*,*) 'PSB_DSCALL:  error check:' ,err

  l_ov_ix=0
  l_ov_el=0
  i = 1
  do while (temp_ovrlap(i) /= -1) 
     idx = temp_ovrlap(i)
     i=i+1
     nprocs = temp_ovrlap(i)
     i = i + 1
     l_ov_ix = l_ov_ix+3*(nprocs-1)
     l_ov_el = l_ov_el + 2
     i = i + nprocs     
  enddo

  l_ov_ix = l_ov_ix+3  
  l_ov_el = l_ov_el+3

  if (debug) write(*,*) 'PSB_DSCALL: Ov len',l_ov_ix,l_ov_el
  allocate(ov_idx(l_ov_ix),ov_el(l_ov_el), stat=info)
  if (info /= no_err) then
     info=4010
     char_err='psb_realloc'
     err=info
     call psb_errpush(err,name,a_err=char_err)
     goto 9999
  end if

  l_ov_ix=0
  l_ov_el=0
  i = 1
  do while (temp_ovrlap(i) /= -1) 
     idx = temp_ovrlap(i)
     i   = i+1
     nprocs = temp_ovrlap(i)
     ov_el(l_ov_el+1)  = idx
     ov_el(l_ov_el+2)  = nprocs
     l_ov_el           = l_ov_el+2
     do j=1, nprocs
        if (temp_ovrlap(i+j) /= myrow) then
           ov_idx(l_ov_ix+1) = temp_ovrlap(i+j)
           ov_idx(l_ov_ix+2) = 1
           ov_idx(l_ov_ix+3) = idx
           l_ov_ix = l_ov_ix+3
        endif
     enddo
     i = i + nprocs +1
  enddo
  l_ov_el         = l_ov_el + 1
  ov_el(l_ov_el)  = -1
  l_ov_ix         = l_ov_ix + 1
  ov_idx(l_ov_ix) = -1

  desc_a%ovrlap_index => ov_idx
  desc_a%ovrlap_elem  => ov_el
  deallocate(prc_v,temp_ovrlap,stat=info)
  if (info /= no_err) then 
     info=4000
     err=info
     call psb_errpush(err,name)
     Goto 9999
  endif
  ! estimate local cols number 
  loc_col=int((psb_colrow_+1.d0)*loc_row)+1  
  allocate(desc_a%loc_to_glob(loc_col),&
       &desc_a%lprm(1),stat=info)  
  if (info /= 0) then 
    call psb_errpush(4010,name,a_err='Allocate')
    goto 9999      
  end if

  ! set LOC_TO_GLOB array to all "-1" values
  desc_a%lprm(1) = 0
  desc_a%loc_to_glob(:) = -1
  do i=1,m
     k = desc_a%glob_to_loc(i) 
     if (k.gt.0) then 
        desc_a%loc_to_glob(k) = i
     endif
  enddo
  nullify(desc_a%bnd_elem,desc_a%halo_index)

!!$  if (debug) write(*,*) 'PSB_DSCALL:  Last bits in desc_a', loc_row,k
  ! set fields in desc_a%MATRIX_DATA....
  desc_a%matrix_data(psb_n_row_)  = loc_row
  desc_a%matrix_data(psb_n_col_)  = loc_row

  call psb_realloc(1,desc_a%halo_index, info)
  if (info /= no_err) then
     info=2025
     char_err='psb_realloc'
     call psb_errpush(err,name,a_err=char_err)
     Goto 9999
  end if

  desc_a%halo_index(:) = -1


  desc_a%matrix_data(psb_m_)        = m
  desc_a%matrix_data(psb_n_)        = n
  desc_a%matrix_data(psb_dec_type_) = psb_desc_bld_
  desc_a%matrix_data(psb_ctxt_)     = icontxt
  call blacs_get(icontxt,10,desc_a%matrix_data(psb_mpi_c_))

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error(icontxt)
     return
  end if
  return

end subroutine psb_dscall
