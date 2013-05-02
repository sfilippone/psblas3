!!$  
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010, 2012
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
! File: psb_icdasb.f90
!
! Subroutine: psb_icdasb
!   Assemble the psblas communications descriptor: inner part.
!   The user callable routine is defined in the psb_tools_mod module.
! 
! Arguments: 
!    desc  - type(psb_desc_type).    The communication descriptor.
!    info    - integer.                return code.
!    ext_hv  - logical                 Essentially this distinguishes a call 
!                                      coming from the build of an extended
!                                      halo descriptor with respect to a normal call. 
!
subroutine psb_icdasb(desc,info,ext_hv)
  use psb_base_mod, psb_protect_name => psb_icdasb
  use psi_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif
  !...Parameters....
  type(psb_desc_type), intent(inout) :: desc
  integer(psb_ipk_), intent(out)               :: info
  logical, intent(in), optional      :: ext_hv

  !....Locals....
  integer(psb_ipk_) ::  int_err(5)
  integer(psb_ipk_),allocatable ::  ovrlap_index(:),halo_index(:), ext_index(:)

  integer(psb_ipk_)  ::  i, n_col, dectype, err_act, n_row,j
  integer(psb_mpik_) ::  np,me, icomm, ictxt,proc_to_comm,iret,bfsz
  logical             :: ext_hv_
  integer(psb_ipk_) :: debug_level, debug_unit
  integer     	:: totxch, idxr, idxs, data_, pnti, snd_pt, rcv_pt,nerv,nesd,idx_pt
  integer(psb_mpik_), allocatable :: blens(:), new_idx(:)
  integer(psb_ipk_), pointer             :: idx(:)
  character(len=20)   :: name

  info = psb_success_
  int_err(1) = 0
  name = 'psb_cdasb'

  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt   = desc%get_context()
  dectype = desc%get_dectype()
  n_row   = desc%get_local_rows()
  n_col   = desc%get_local_cols()
  call psb_get_mpicomm(ictxt,icomm )

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.psb_is_ok_desc(desc)) then 
    info = psb_err_spmat_invalid_state_
    int_err(1) = dectype
    call psb_errpush(info,name)
    goto 9999
  endif

  info = psb_get_errstatus()
  if (info /= psb_success_) then 
    ! Something went wrong in cdins/spins
    ! signal and exit
    info = psb_err_wrong_ins_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (present(ext_hv)) then 
    ext_hv_ = ext_hv
  else
    ext_hv_ = .false.
  end if
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit, *) me,' ',trim(name),': start'

  if (allocated(desc%indxmap)) then 
    call psi_ldsc_pre_halo(desc,ext_hv_,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='ldsc_pre_halo')
      goto 9999
    end if

    ! Take out the lists for ovrlap, halo and ext...
    call psb_move_alloc(desc%ovrlap_index,ovrlap_index,info)
    call psb_move_alloc(desc%halo_index,halo_index,info)
    call psb_move_alloc(desc%ext_index,ext_index,info)

    if (debug_level >= psb_debug_ext_) &
         & write(debug_unit,*) me,' ',trim(name),': Final conversion'
    ! Then convert and put them back where they belong.    
    call psi_cnv_dsc(halo_index,ovrlap_index,ext_index,desc,info) 

    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='psi_cnv_dsc')
      goto 9999
    end if

    deallocate(ovrlap_index, halo_index, ext_index, stat=info)
    if (info /= psb_success_) then
      info =psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    end if

    call desc%indxmap%asb(info)
    if (info == psb_success_) then 
      if (allocated(desc%indxmap%tempvg)) &
           & deallocate(desc%indxmap%tempvg,stat=info)
    end if
    if (info /= psb_success_) then 
      write(0,*) 'Error from internal indxmap asb ',info
      info = psb_success_
    end if

  else
    info = psb_err_spmat_invalid_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  !datatypes allocation
  data_ = psb_comm_halo_
  call desc%get_list(data_,idx,totxch,idxr,idxs,info)
  call psb_realloc(max(1,totxch),psb_nkidx_,desc%sendtypes,info)
  if (info == 0) call psb_realloc(max(1,totxch),psb_nkidx_,desc%recvtypes,info)
  if (info /= 0) then 
    write(0,*) 'Failed alloc  send/recvtypes',totxch,psb_nkidx_,info
    info =psb_err_alloc_dealloc_
    call psb_errpush(info,name)
    goto 9999
  end if

  ! Init here, they will be filled in upon request
  desc%sendtypes(:,:) = mpi_datatype_null
  desc%recvtypes(:,:) = mpi_datatype_null
  pnti   = 1
  bfsz   = 0
  do i=1, totxch
    proc_to_comm = idx(pnti+psb_proc_id_)
    nerv = idx(pnti+psb_n_elem_recv_)
    nesd = idx(pnti+nerv+psb_n_elem_send_)
    bfsz = max(bfsz,nesd,nerv)
    pnti   = pnti + nerv + nesd + 3      
  end do
  bfsz = max(1,bfsz)
  call psb_realloc(bfsz,blens,info)
  if (info == 0) call psb_realloc(bfsz,new_idx,info)
  if(info /= psb_success_) then
    write(0,*) 'Failed alloc  blens/new_idx',bfsz,info
    call psb_errpush(psb_err_alloc_dealloc_,name)
    goto 9999
  end if

  !We've got to set the derivate datatypes
  !Send/Gather
  pnti   = 1
  snd_pt = 1
  do i=1, totxch
    nerv = idx(pnti+psb_n_elem_recv_)
    nesd = idx(pnti+nerv+psb_n_elem_send_)
    idx_pt = 1+pnti+psb_n_elem_recv_
    do j=1, nerv
      blens(j)   = 1
      new_idx(j) = idx(idx_pt+j-1)-1
    end do
    call psb_mpi_type(nerv,blens,new_idx,&
         & psb_mpi_ipk_integer,desc%recvtypes(i,psb_ipkidx_),iret)
    call psb_mpi_type(nerv,blens,new_idx,&
         & psb_mpi_def_integer,desc%recvtypes(i,psb_mpikidx_),iret)
    call psb_mpi_type(nerv,blens,new_idx,&
         & psb_mpi_lng_integer,desc%recvtypes(i,psb_lngkidx_),iret)
    call psb_mpi_type(nerv,blens,new_idx,&
         & psb_mpi_r_spk_,desc%recvtypes(i,psb_rspkidx_),iret)
    call psb_mpi_type(nerv,blens,new_idx,&
         & psb_mpi_r_dpk_,desc%recvtypes(i,psb_rdpkidx_),iret)
    call psb_mpi_type(nerv,blens,new_idx,&
         & psb_mpi_c_spk_,desc%recvtypes(i,psb_cspkidx_),iret)
    call psb_mpi_type(nerv,blens,new_idx,&
         & psb_mpi_c_dpk_,desc%recvtypes(i,psb_cdpkidx_),iret)


    idx_pt = 1+pnti+nerv+psb_n_elem_send_
    do j=1,nesd
      blens(j)   = 1
      new_idx(j) = idx(idx_pt+j-1)-1
    end do
    call psb_mpi_type(nesd,blens,new_idx,&
         & psb_mpi_ipk_integer,desc%sendtypes(i,psb_ipkidx_),iret)
    call psb_mpi_type(nesd,blens,new_idx,&
         & psb_mpi_def_integer,desc%sendtypes(i,psb_mpikidx_),iret)
    call psb_mpi_type(nesd,blens,new_idx,&
         & psb_mpi_lng_integer,desc%sendtypes(i,psb_lngkidx_),iret)
    call psb_mpi_type(nesd,blens,new_idx,&
         & psb_mpi_r_spk_,desc%sendtypes(i,psb_rspkidx_),iret)
    call psb_mpi_type(nesd,blens,new_idx,&
         & psb_mpi_r_dpk_,desc%sendtypes(i,psb_rdpkidx_),iret)
    call psb_mpi_type(nesd,blens,new_idx,&
         & psb_mpi_c_spk_,desc%sendtypes(i,psb_cspkidx_),iret)
    call psb_mpi_type(nesd,blens,new_idx,&
         & psb_mpi_c_dpk_,desc%sendtypes(i,psb_cdpkidx_),iret)

    pnti   = pnti + nerv + nesd + 3
  end do

  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': Done'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return

contains
  subroutine psb_mpi_type(nitem,disp,idx,type,newtype,iret)
    integer(psb_mpik_) :: nitem, disp(:),idx(:),type,newtype,iret
    call mpi_type_indexed(nitem,disp,idx,type,newtype,iret)
    if (iret /= 0) &
         &  write(0,*) 'From mpi_type_indexed: ',iret,type
    call mpi_type_commit(newtype,iret) 
    if (iret /= 0) &
         &  write(0,*) 'From mpi_type_commit: ',iret,newtype
  end subroutine psb_mpi_type

end subroutine psb_icdasb
