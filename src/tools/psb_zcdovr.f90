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
! File: psb_zcdovr.f90
!
! Subroutine: psb_zcdovr
!    This routine takes a matrix A with its descriptor, and builds the 
!    auxiliary descriptor corresponding to the number of overlap levels
!    specified on input. It really is just a size estimation/allocation
!    front end for <psb_zcdovrbld>.
! 
! Parameters: 
!    a        - type(<psb_zspmat_type>).       The input sparse matrix.
!    desc_a   - type(<psb_desc_type>).         The input communication descriptor.
!    norv     - integer.                       The number of overlap levels.
!    desc_ov  - type(<psb_desc_type>).         The auxiliary output communication descriptor.
!    info     - integer.                       Eventually returns an error code.
!
Subroutine psb_zcdovr(a,desc_a,novr,desc_ov,info)

  use psb_serial_mod
  use psb_descriptor_type
  Use psb_prec_type
  Use psb_prec_mod
  use psb_error_mod
  Implicit None

  !     .. Array Arguments ..
  integer, intent(in)                :: novr
  Type(psb_zspmat_type), Intent(in)  ::  a
  Type(psb_desc_type), Intent(in)    :: desc_a
  Type(psb_desc_type), Intent(inout) :: desc_ov
  integer, intent(out)               :: info

  real(kind(1.d0)) :: t1,t2,t3,mpi_wtime
  external  mpi_wtime
  integer   idscb,idsce,iovrb,iovre, ierr, irank, icomm, err_act
!!$  integer mpe_log_get_event_number,mpe_Describe_state,mpe_log_event

  interface psb_cdcpy
     subroutine psb_cdcpy(desc_in,desc_out,info)
       use psb_descriptor_type
       type(psb_desc_type), intent(in)  :: desc_in
       type(psb_desc_type), intent(out) :: desc_out
       integer, intent(out)             :: info
     end subroutine psb_cdcpy
  end interface

  interface psb_cdovrbld
     subroutine psb_zcdovrbld(n_ovr,desc_p,desc_a,a,l_tmp_halo,&
          & l_tmp_ovr_idx,lworks,lworkr,info)
       use psb_prec_type
       use psb_spmat_type
       type(psb_zspmat_type),intent(in)     :: a
       type(psb_desc_type),intent(in)       :: desc_a
       type(psb_desc_type),intent(inout)    :: desc_p
       integer,intent(in)                   :: n_ovr
       integer, intent(in)                  :: l_tmp_halo,l_tmp_ovr_idx
       integer, intent(inout)               :: lworks, lworkr
       integer, intent(out)                 :: info
     end subroutine psb_zcdovrbld
  end interface



  !     .. Local Scalars ..
  Integer ::  i, j, k, nprow,npcol, me, mycol,m,nnzero,&
       &  ictxt, lovr, lelem,lworks,lworkr, n_col, int_err(5),&
       &  n_row,index_dim,elem_dim, l_tmp_ovr_idx,l_tmp_halo, nztot,nhalo
  Logical,Parameter :: debug=.false.
  character(len=20)                 :: name, ch_err

  name='psb_cdovr'
  info  = 0
  call psb_erractionsave(err_act)

  ictxt=desc_a%matrix_data(psb_ctxt_)

  Call blacs_gridinfo(ictxt,nprow,npcol,me,mycol)

  If(debug) Write(0,*)'in psb_cdovr',novr

  m=desc_a%matrix_data(psb_n_row_)
  nnzero=Size(a%aspk)
  n_col=desc_a%matrix_data(psb_n_col_)
  nhalo = n_col-m
  If(debug) Write(0,*)'IN CDOVR1',novr ,m,nnzero,n_col
  if (novr<0) then
     info=10
     int_err(1)=1
     int_err(2)=novr
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  endif

  if (novr==0) then 
    !
    ! Just copy the input.  
    !
    if (debug) write(0,*) 'Calling desccpy'
    call psb_cdcpy(desc_a,desc_ov,info)
    if (info /= 0) then
       info=4010
       ch_err='psb_cdcpy'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if
    if (debug) write(0,*) 'From desccpy'
    return
  endif

  call blacs_get(ictxt,10,icomm )
!!$    call MPI_Comm_rank(icomm,irank,ierr)
!!$    idscb  = mpe_log_get_event_number()
!!$    idsce  = mpe_log_get_event_number()
!!$    iovrb  = mpe_log_get_event_number()
!!$    iovre  = mpe_log_get_event_number()
!!$    if (irank==0) then 
!!$      info = mpe_describe_state(idscb,idsce,"CDASB ","NavyBlue")
!!$      info = mpe_describe_state(iovrb,iovre,"CDOVRR ","DeepPink")
!!$    endif
  If(debug)then 
    Write(0,*)'BEGIN cdovr',me,nhalo
    call blacs_barrier(ictxt,'All')
  endif
  t1 = mpi_wtime()



!!$      ierr = MPE_Log_event( idscb, 0, "st CDASB" )
  !
  ! Ok, since we are only estimating, do it as follows: 
  ! LOVR= (NNZ/NROW)*N_HALO*N_OVR  This assumes that the local average 
  ! nonzeros per row is the same as the global. 
  !
  call psb_spinfo(psb_nztotreq_,a,nztot,info)
  if (info /= 0) then
     info=4010
     ch_err='psb_spinfo'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  if (nztot>0) then 
     lovr   = ((nztot+m-1)/m)*nhalo*novr
     lworks = ((nztot+m-1)/m)*nhalo
     lworkr = ((nztot+m-1)/m)*nhalo
  else
     info=-1
     call psb_errpush(info,name)
     goto 9999
  endif
  If(debug)Write(0,*)'ovr_est done',me,novr,lovr
  index_dim = size(desc_a%halo_index)
  elem_dim  = size(desc_a%halo_index)

  allocate(desc_ov%ovrlap_index(novr*(Max(2*index_dim,1)+1)),&
       &   desc_ov%ovrlap_elem(novr*(Max(elem_dim,1)+3)),&
       &   desc_ov%matrix_data(psb_mdata_size_),&
       &   desc_ov%halo_index(novr*(Size(desc_a%halo_index)+3)),STAT=INFO)
  if (info /= 0) then
     info=4000
     call psb_errpush(info,name)
     goto 9999
  end if

  l_tmp_ovr_idx=novr*(3*Max(2*index_dim,1)+1)
  l_tmp_halo=novr*(3*Size(desc_a%halo_index))

  desc_ov%ovrlap_index(:)   = -1
  desc_ov%ovrlap_elem(:)    = -1
  desc_ov%halo_index(:)     = -1
  desc_ov%matrix_data(1:10) = desc_a%matrix_data(1:10)
  desc_ov%matrix_data(psb_dec_type_) = psb_desc_bld_ 

  Allocate(desc_ov%loc_to_glob(Size(desc_a%loc_to_glob)),&
       & desc_ov%glob_to_loc(Size(desc_a%glob_to_loc)),stat=info)
  if (info /= 0) then
     info=4000
     call psb_errpush(info,name)
     goto 9999
  end if

  desc_ov%loc_to_glob(:) = desc_a%loc_to_glob(:)
  desc_ov%glob_to_loc(:) = desc_a%glob_to_loc(:)
  If(debug)then 
    Write(0,*)'Start cdovrbld',me,lworks,lworkr
    call blacs_barrier(ictxt,'All')
  endif
  
  !
  ! The real work goes on in here....
  !
  Call psb_cdovrbld(novr,desc_ov,desc_a,a,&
       & l_tmp_halo,l_tmp_ovr_idx,lworks,lworkr,info) 
  if (info /= 0) then
     info=4010
     ch_err='psb_cdovrbld'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  desc_ov%matrix_data(psb_dec_type_) = psb_desc_asb_
  If(debug)then 
    Write(0,*)'Done cdovrbld',me,lworks,lworkr
    call blacs_barrier(ictxt,'All')
  endif
!!$      ierr = MPE_Log_event( idsce, 0, "st CDASB" )

  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == act_abort) then
     call psb_error(ictxt)
     return
  end if
  Return

End Subroutine psb_zcdovr

