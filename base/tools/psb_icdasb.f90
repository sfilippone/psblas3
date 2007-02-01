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
! File: psb_cdasb.f90
!
! Subroutine: psb_cdasb
!   Assembly the psblas communications descriptor.
! 
! Parameters: 
!    desc_a  - type(<psb_desc_type>).         The communication descriptor.
!    info    - integer.                       Eventually returns an error code.
subroutine psb_icdasb(desc_a,info,ext_hv)
  use mpi
  use psb_descriptor_type
  use psb_serial_mod
  use psb_const_mod
  use psi_mod
  use psb_error_mod
  use psb_penv_mod
  implicit none
  !...Parameters....
  type(psb_desc_type), intent(inout) :: desc_a
  integer, intent(out)               :: info
  logical, intent(in), optional      :: ext_hv

  !....Locals....
  integer          ::  int_err(5), itemp(2)
  integer,allocatable ::  ovrlap_index(:),halo_index(:), ext_index(:)

  integer          ::  i,j,err,np,me,lovrlap,lhalo,nhalo,novrlap,max_size,&
       & max_halo,n_col,ldesc_halo, ldesc_ovrlap, dectype, err_act, &
       & key, ih, nh, idx, nk,icomm,hsize
  integer                       :: ictxt,n_row
  logical                       :: ext_hv_
  logical, parameter            :: debug=.false., debugwrt=.false.
  character(len=20)             :: name,ch_err

  info = 0
  int_err(1) = 0
  name = 'psb_cdasb'

  call psb_erractionsave(err_act)

  ictxt   = psb_cd_get_context(desc_a)
  dectype = psb_cd_get_dectype(desc_a)
  n_row   = psb_cd_get_local_rows(desc_a)
  n_col   = psb_cd_get_local_cols(desc_a)
  call psb_get_mpicomm(ictxt,icomm )

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.psb_is_ok_desc(desc_a)) then 
    info = 600
    int_err(1) = dectype
    call psb_errpush(info,name)
    goto 9999
  endif
  
  if (present(ext_hv)) then 
    ext_hv_ = ext_hv
  else
    ext_hv_ = .false.
  end if
  if (debug) write (0, *) '   Begin matrix assembly...'

  if (psb_is_bld_desc(desc_a)) then 
    if (debug) write(0,*) 'psb_cdasb: Checking rows insertion'
    ! check if all local row are inserted
    do i=1,psb_cd_get_local_cols(desc_a)
      if (desc_a%loc_to_glob(i) < 0) then
        info=3100
        exit
      endif
    enddo

    if (info /= psb_no_err_) then    
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    endif

    call psb_realloc(psb_cd_get_local_cols(desc_a),desc_a%loc_to_glob,info)

    if (psb_is_large_desc(desc_a)) then 
      call psi_ldsc_pre_halo(desc_a,ext_hv_,info)
    end if

    call psb_transfer(desc_a%ovrlap_index,ovrlap_index,info)
    call psb_transfer(desc_a%halo_index,halo_index,info)
    call psb_transfer(desc_a%ext_index,ext_index,info)

    call psi_cnv_dsc(halo_index,ovrlap_index,ext_index,desc_a,info) 
    if (info /= 0) then
      call psb_errpush(4010,name,a_err='psi_cnv_dsc')
      goto 9999
    end if


    deallocate(ovrlap_index, halo_index, ext_index, stat=info)
    if (info /= 0) then
      info =4000
      call psb_errpush(info,name)
      goto 9999
    end if
    ! Finally, cleanup the AVL tree, as it is really only needed 
    ! when building. 
    if (allocated(desc_a%ptree)) then 
      call FreePairSearchTree(desc_a%ptree)   
      deallocate(desc_a%ptree,stat=info)
      if (info /= 0) then 
        info=2059
        call psb_errpush(info,name)
        goto 9999
      end if
    end if
    ! Ok, register into MATRIX_DATA &  free temporary work areas
    desc_a%matrix_data(psb_dec_type_) = psb_desc_asb_
  else
    info = 600
    call psb_errpush(info,name)
    goto 9999
    if (debug) write(0,*) 'dectype 2 :',psb_cd_get_dectype(desc_a),&
         &psb_desc_bld_,psb_desc_asb_,psb_desc_upd_
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_icdasb
