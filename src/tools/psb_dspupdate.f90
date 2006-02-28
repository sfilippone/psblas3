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
! File: psb_dspupdate.f90
!
! Subroutine: psb_dspupdate
!    Updates a sparse matrix.
! 
! Parameters: 
!    a        - type(<psb_dspmat_type>).              
!    ia       - integer, dimension(:).
!    ja       - integer, dimension(:).
!    blck     - type(<psb_dspmat_type>).              
!    desc_a   - type(<psb_desc_type>).            
!    info     - integer.                          
!    ix       - integer(optional).
!    jx       - integer(optional).
!    updflag  - integer(optional).
!
subroutine psb_dspupdate(a, ia, ja, blck, desc_a,info,ix,jx,updflag)

  use psb_descriptor_type
  use psb_spmat_type
  use psb_serial_mod
  use psb_error_mod
  implicit none

  !....parameters...
  type(psb_desc_type), intent(in)      ::  desc_a
  type(psb_dspmat_type), intent(inout) ::  a
  integer, intent(in)                  ::  ia,ja
  type(psb_dspmat_type), intent(in)    ::  blck
  integer, intent(out)                 ::  info
  integer, optional, intent(in)        ::  ix,jx
  integer, optional, intent(in)        ::  updflag

  !locals.....


  interface
     subroutine dcsupd(m,n,fida,descra,a,ia1,ia2,infoa,ia,ja,&
          & fidh,descrh,h,ih1,ih2,infoh,ih,jh,&
          & flag,glob_to_loc,iwork,liwork,ierror)
       implicit none
       !      .. scalar arguments ..
       integer, intent(in) :: m, n, liwork,ia,ja,ih,jh, flag
       integer, intent(out) :: ierror
       !      .. array arguments ..
       double precision, intent(in)    :: h(*)
       double precision, intent(inout) :: a(*)
       integer, intent(in) :: ih1(*), ih2(*), infoh(10), glob_to_loc(*)
       integer, intent(inout) :: ia1(*), ia2(*), infoa(10), iwork(*)
       character, intent(in)  :: fida*5, fidh*5,descra*11, descrh*11

     end subroutine dcsupd
  end interface

  integer :: icontxt,i,loc_row,prec_loc_row ,glob_row,row,&
       & k ,start_row,end_row,first_loc_row,n_row,j,int_err(5),&
       &locix,locjx,allocated_prcv, dectype, flag,err_act,err
  integer,pointer        :: prcv(:),gtl(:), ltg(:)
  integer                :: nprow,npcol, myrow ,mycol, lr, lc, nrow,ncol
  integer                ::  m,n, iupdflag
  integer,pointer        ::  iworkaux(:)    
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='psb_dspupdate'
  call psb_erractionsave(err_act)

  if (present(ix)) then
     locix=ix
  else
     locix=1
  endif

  if (present(updflag)) then
     iupdflag = updflag
  else
     iupdflag = psb_upd_glb_
  endif

  if (present(jx)) then
     locjx=jx
  else
     locjx=1
  endif
  icontxt=desc_a%matrix_data(psb_ctxt_)

  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
  if (nprow.eq.-1) then
     info = 2010
     call psb_errpush(info,name)
     goto 9999
  else if (npcol.ne.1) then
     info = 2030
     int_err(1) = npcol
     call psb_errpush(info,name,int_err)
     goto 9999
  endif

  gtl => desc_a%glob_to_loc
  ltg => desc_a%loc_to_glob
  nrow = desc_a%matrix_data(psb_n_row_)
  ncol = desc_a%matrix_data(psb_n_col_)
  dectype = desc_a%matrix_data(psb_dec_type_)
  ! check if a is already allocated (called psdalloc)
  if (.not.psb_is_upd_dec(dectype)) then
     info = 290
     int_err(1) = dectype
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  endif
  allocate(prcv(nprow),iworkaux(3*ncol+4),stat=info)
  if (info.ne.0) then
     info = 2023
     int_err(1) = max(1,nprow,3*ncol+4)
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  endif

  flag = 2

  m = blck%m
  n = blck%k

  if (iupdflag == psb_upd_glb_) then 

     row = ia
     i   = 1
     blckr: do while (i.le.m)
        !loop over all blck's rows

        ! row actual block row 
        row      = locix+i-1
        glob_row = ia+i-1

        lr = gtl(glob_row)

        if ((1 <= lr) .and. (lr <= nrow)) then
           ! at least one row belongs to me 

           start_row=row
           do 
              ! loop until actual row belong to me
              ! and all actual row to insert are ordered

              prec_loc_row=loc_row

              ! --if  loc_row is != -1 is already assigned
              !   local index to globrow whith value loc_row
              ! --if loc:row == -1 it isn't assigned local row to
              !   glob_row
              loc_row=gtl(glob_row)
              if (start_row.eq.i) first_loc_row=loc_row
              ! next blck's row
              i=i+1
              if (i.le.m) then
                 row=locix+i-1
                 glob_row=ia+i-1
                 k = gtl(glob_row)
                 if ((.not.((1 <= lr) .and. (lr <= nrow)))&
                      & .or.((prec_loc_row+1.ne.loc_row).and.&
                      & (start_row+1.ne.i)).or.(i.gt.m)) exit
              else
                 exit
              endif
           enddo

           end_row=i-1
           ! insert blck submatrix
           call dcsupd(end_row-start_row+1,n,a%fida,a%descra,a%aspk,&
       	        & a%ia1,a%ia2,a%infoa,first_loc_row, ja, blck%fida ,&
		& blck%descra,blck%aspk,blck&
 		& %ia1,blck%ia2,blck%infoa,start_row, locjx, flag,&
		& desc_a%glob_to_loc,&
		& iworkaux, size(iworkaux),info)
           if (info.ne.0) exit blckr
        endif
        ! next blck's row
        i=i+1
     enddo blckr

     if (info.ne.0) then
        info = 4010
        ch_err='dcsupd'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     endif


  else if (iupdflag == psb_upd_loc_) then

     ! insert blck submatrix
     call dcsupd(m,n,a%fida,a%descra,a%aspk,&
          & a%ia1,a%ia2,a%infoa,ia, ja, blck%fida ,&
          & blck%descra,blck%aspk,blck&
          & %ia1,blck%ia2,blck%infoa,locix,locjx, flag,&
          & desc_a%glob_to_loc,&
          & iworkaux, size(iworkaux),info)

     if (info.ne.0) then
        info = 4010
        ch_err='dcsupd'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     endif
  else
     ! fix next error code
     info = 999
     call psb_errpush(info,name)
     goto 9999
  endif
  deallocate(prcv,iworkaux,stat=info)
  if (info.ne.0) then
     info=2040
     call psb_errpush(info,name)
     goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_ret) then
     return
  else
     call psb_error(icontxt)
  end if
  return

end subroutine psb_dspupdate








