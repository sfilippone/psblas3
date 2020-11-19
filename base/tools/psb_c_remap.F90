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
!
! Subroutine: psb_c_remap
! 
! Arguments: 
!    desc_in  - type(psb_desc_type).       The communication descriptor to be cloned.
!    desc_out - type(psb_desc_type).       The output communication descriptor.
!    info     - integer.                       Return code.
subroutine psb_c_remap(np_remap, desc_in, a_in, ipd, isrc, nrsrc, naggr, &
     & desc_out, a_out, info)

  use psb_base_mod, psb_protect_name => psb_c_remap

  implicit none
  !....parameters...
  integer(psb_ipk_), intent(in)        :: np_remap
  type(psb_desc_type), intent(inout)   :: desc_in
  type(psb_cspmat_type), intent(inout) :: a_in
  type(psb_cspmat_type), intent(out)   :: a_out
  type(psb_desc_type), intent(out)     :: desc_out
  integer(psb_ipk_), intent(out)       :: ipd 
  integer(psb_ipk_), allocatable, intent(out) :: isrc(:), nrsrc(:), naggr(:)
  integer(psb_ipk_), intent(out)       :: info


  ! locals
  type(psb_ctxt_type) :: ctxt, newctxt
  integer(psb_ipk_) :: np, me, err_act
  integer(psb_ipk_) :: rnp, rme
  integer(psb_ipk_) :: ipdest, id1, id2, imd, i, nsrc
  integer(psb_ipk_), allocatable :: newnl(:), nzsrc(:), ids(:) 
  type(psb_lc_coo_sparse_mat) :: acoo_snd, acoo_rcv
  integer(psb_ipk_) :: debug_level, debug_unit  
  character(len=20)   :: name

  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  if (psb_get_errstatus() /= 0) return 
  info = psb_success_
  call psb_erractionsave(err_act)
  name = 'psb_cdcpy'

  ctxt = desc_in%get_context()

  ! check on blacs grid 
  call psb_info(ctxt, me, np)
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': Entered'
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

!!$  write(0,*) ' Remapping from ',np,' onto ', np_remap

  if (desc_in%get_fmt() == 'BLOCK') then
    !
    ! Should we spread the processes in the new context,
    ! or should we keep them close? 
    ! 
    if (.true.) then 
      allocate(ids(0:np_remap-1))
      if (np_remap <= np/2) then
        ids(0) = 0
        do ipdest=1,np_remap -1
          ids(ipdest) = ids(ipdest-1) + np/np_remap
        end do
!!$        write(0,*) ' IDS ',ids(:) 
      else
        do ipdest = 0, np_remap-1
          ids(ipdest) = ipdest
        end do
      end if
      call psb_init(newctxt,np=np_remap,basectxt=ctxt,ids=ids)
    else
      call psb_init(newctxt,np=np_remap,basectxt=ctxt)
    end if

    call psb_info(newctxt,rme,rnp)
!!$    write(0,*) 'Old context: ',me,np,' New context: ',rme,rnp
    call psb_bcast(ctxt,rnp)
    allocate(newnl(rnp),naggr(np),stat=info)
    if (info /= 0) then
      info = psb_err_alloc_dealloc_
      call psb_errpush(info,name)
      goto 9999
    endif
    if (rnp >= np) then 
      write(0,*) ' No remapping on larger proc count now'
      info = psb_err_internal_error_
      call psb_errpush(info,name)
      goto 9999
    end if
    naggr = 0

    !
    ! Compute destination for my data.
    ! Simplistic reallocation: divide the NP processes
    ! across the new ones (as balanced as possible),
    ! then send all data from old to new process
    !
    id2 = np/rnp
    id1 = id2+1
    imd = mod(np,rnp)
    if (me < (imd*id1)) then
      ipdest = (me/id1)
    else
      ipdest = ( ((me-imd*id1)/id2) +  imd)
    end if
    if (allocated(ids)) then
      ipd = ids(ipdest)
    else
      ipd = ipdest
    end if
!!$    write(0,*) ' Sending my data from ',me,' to ', &
!!$         & ipd, 'out of ',rnp,rnp-1

    !
    ! Compute local rows for all new
    ! processes; will have a BLOCK distribution
    ! 
    newnl = 0
    newnl(ipdest+1) = desc_in%get_local_rows()    
    call psb_sum(ctxt,newnl)

    if (rme>=0) then
      ! 
      if (rme < imd) then
        isrc = [ (i, i=rme*id1,min(rme*id1+id1-1,np-1)) ]
      else
        isrc = [ (i, i=  imd*id1+((rme-imd))*id2,&
             & min(imd*id1+(rme-imd)*id2+id2-1,np-1) ) ]            
      end if
!!$      write(0,*) me,rme,imd,' ISRC: ',isrc(:)
      nsrc = size(isrc)
!!$      write(0,*) me,rme,'In ',desc_in%get_local_rows(),desc_in%get_global_rows(),&
!!$           & ' out ',desc_out%get_local_rows(),desc_out%get_global_rows()
    else 
!!$      write(0,*) me,rme,'In ',desc_in%get_local_rows(),desc_in%get_global_rows(),&
!!$           & ' out ',0,0
    end if
  else
    write(0,*) 'Right now only BLOCK on input '
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  endif


  !
  ! Collect matrices on their destinations
  !
  block
    integer(psb_ipk_) :: nzsnd, nzrcv, ip
    integer(psb_ipk_) :: nrl, ncl, nzl, nzp
    call a_in%cp_to(acoo_snd)
    nzsnd = acoo_snd%get_nzeros()
    call psb_snd(ctxt,nzsnd,ipd)
    call psb_snd(ctxt,desc_in%get_local_rows(),ipd)
    ! Convert to global numbering
    call psb_loc_to_glob(acoo_snd%ia(1:nzsnd),desc_in,info)
    call psb_loc_to_glob(acoo_snd%ja(1:nzsnd),desc_in,info)

    call psb_snd(ctxt,acoo_snd%ia(1:nzsnd),ipd)
    call psb_snd(ctxt,acoo_snd%ja(1:nzsnd),ipd)
    call psb_snd(ctxt,acoo_snd%val(1:nzsnd),ipd)

    if (rme>=0) then
      ! prepare to receive
      nzsrc = isrc
      nrsrc = isrc
      nzl = 0
      do ip=1, nsrc
        call psb_rcv(ctxt,nzsrc(ip),isrc(ip))
        call psb_rcv(ctxt,nrsrc(ip),isrc(ip))
        nzl = nzl + nzsrc(ip)
      end do
!!$      write(0,*) rme,' Check on NR:',newnl(rme+1),sum(nrsrc)
      call acoo_rcv%allocate(newnl(rme+1),newnl(rme+1),nzl)
      nrl = acoo_rcv%get_nrows()
      ncl = acoo_rcv%get_ncols()
      nzp = 0
      do ip=1, nsrc
        call psb_rcv(ctxt,acoo_rcv%ia(nzp+1:nzp+nzsrc(ip)),isrc(ip))
        call psb_rcv(ctxt,acoo_rcv%ja(nzp+1:nzp+nzsrc(ip)),isrc(ip))
        call psb_rcv(ctxt,acoo_rcv%val(nzp+1:nzp+nzsrc(ip)),isrc(ip))
        nzp = nzp + nzsrc(ip)
      end do
      call acoo_rcv%set_nzeros(nzp)
!!$      write(0,*) rme,' Collected: ',&
!!$           & acoo_rcv%get_nrows(),acoo_rcv%get_ncols(),acoo_rcv%get_nzeros()

      !
      !  New descriptor
      !    
      call psb_cdall(newctxt,desc_out,info,nl=newnl(rme+1))
      ! Insert
      call psb_spall(a_out,desc_out,info)
      call psb_spins(nzp,acoo_rcv%ia(1:nzp),acoo_rcv%ja(1:nzp),&
           & acoo_rcv%val(1:nzp),a_out,desc_out,info)      
      ! Assemble
      call psb_cdasb(desc_out,info)
      call psb_spasb(a_out,desc_out,info)

!!$      write(0,*) rme,' Regenerated: ',&
!!$           & desc_out%get_local_rows(), desc_out%get_local_cols(),&
!!$           & a_out%get_nrows(),a_out%get_ncols(),a_out%get_nzeros()
      naggr(me+1) = desc_out%get_local_rows()
    else
      naggr(me+1) = 0
    end if
    call psb_sum(ctxt,naggr)
    
  end block

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return  

end subroutine psb_c_remap
