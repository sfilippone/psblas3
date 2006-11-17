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
  use psb_penv_mod
  Implicit None

  !     .. Array Arguments ..
  integer, intent(in)                :: novr
  Type(psb_zspmat_type), Intent(in)  :: a
  Type(psb_desc_type), Intent(in)    :: desc_a
  Type(psb_desc_type), Intent(inout) :: desc_ov
  integer, intent(out)               :: info

  real(kind(1.d0)) :: t1,t2,t3,mpi_wtime
  external  mpi_wtime
  integer   icomm, err_act

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
  Integer ::  i, j, k, np, me,m,nnzero,&
       &  ictxt, lovr, lworks,lworkr, n_col, int_err(5),&
       &  index_dim,elem_dim, l_tmp_ovr_idx,l_tmp_halo, nztot,nhalo
  Logical, parameter :: debug=.false.
  character(len=20)  :: name, ch_err

  name='psb_cdovr'
  info  = 0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)

  Call psb_info(ictxt, me, np)

  If(debug) Write(0,*)'in psb_cdovr',novr

  m=psb_cd_get_local_rows(desc_a)
  nnzero=Size(a%aspk)
  n_col=psb_cd_get_local_cols(desc_a)
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

  call psb_get_mpicomm(ictxt,icomm )

  If(debug)then 
    Write(0,*)'BEGIN cdovr',me,nhalo
    call psb_barrier(ictxt)
  endif
  t1 = mpi_wtime()


  !
  ! Ok, since we are only estimating, do it as follows: 
  ! LOVR= (NNZ/NROW)*N_HALO*N_OVR  This assumes that the local average 
  ! nonzeros per row is the same as the global. 
  !
  nztot = psb_sp_get_nnzeros(a)
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

  call psb_realloc(psb_mdata_size_,desc_ov%matrix_data,info)
  if (info==0) call psb_realloc(novr*(Max(elem_dim,1)+3),desc_ov%ovrlap_elem,info)
  if (info /= 0) then
    info=4000
    call psb_errpush(info,name)
    goto 9999
  end if

  l_tmp_ovr_idx=novr*(3*Max(2*index_dim,1)+1)
  l_tmp_halo=novr*(3*Size(desc_a%halo_index))

  desc_ov%matrix_data(:)    = desc_a%matrix_data(:)
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
  If(debug) then
    Write(0,*)'Start cdovrbld',me,lworks,lworkr
    call psb_barrier(ictxt)
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
  If(debug) then
    Write(0,*)'Done cdovrbld',me,lworks,lworkr
    call psb_barrier(ictxt)
  endif

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

