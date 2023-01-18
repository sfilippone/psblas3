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
! File: psb_sspasb.f90
!
! Subroutine: psb_sspasb
!    Assemble sparse matrix
!
! Arguments: 
!    a        - type(psb_sspmat_type).     The sparse matrix to be allocated.      
!    desc_a   - type(psb_desc_type).       The communication descriptor.
!    info     - integer.                     return code.
!    afmt     - character(optional)          The desired output storage format.
!    upd      - character(optional).         How will the matrix be updated? 
!                                            psb_upd_srch_    Simple strategy  
!                                            psb_upd_perm_    Permutation(more memory)
! 
!
subroutine psb_sspasb(a,desc_a, info, afmt, upd, mold)
  use psb_base_mod, psb_protect_name => psb_sspasb
  use psb_sort_mod
  use psi_mod
  implicit none


  !...Parameters....
  type(psb_sspmat_type), intent (inout)    :: a
  type(psb_desc_type), intent(inout)       :: desc_a
  integer(psb_ipk_), intent(out)           :: info
  integer(psb_ipk_), optional, intent(in)  :: upd
  character(len=*), optional, intent(in)   :: afmt
  class(psb_s_base_sparse_mat), intent(in), optional :: mold
  !....Locals....
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me, err_act
  integer(psb_ipk_) :: n_row,n_col, dupl_
  integer(psb_ipk_) :: debug_level, debug_unit
  character(len=20)     :: name, ch_err
  class(psb_i_base_vect_type), allocatable :: ivm

  info = psb_success_
  name = 'psb_spasb'
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt    = desc_a%get_context()
  n_row    = desc_a%get_local_rows()
  n_col    = desc_a%get_local_cols()

  ! check on BLACS grid 
  call psb_info(ctxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.desc_a%is_asb()) then
    info = psb_err_invalid_cd_state_
    call psb_errpush(info,name)
    goto 9999
  end if

  if (debug_level >= psb_debug_ext_)&
       & write(debug_unit, *) me,' ',trim(name),&
       & '   Begin matrix assembly...'

  !check on errors encountered in psdspins

  if (a%is_bld()) then
    dupl_ = a%get_dupl()
    !
    ! First case: we come from a fresh build. 
    ! 
    if (a%is_remote_build()) then 
      !write(0,*) me,name,' Size of rmta:',a%rmta%get_nzeros()
      block
        type(psb_ls_coo_sparse_mat) :: a_add
        integer(psb_ipk_), allocatable :: ila(:), jla(:)
        integer(psb_ipk_) :: nz, nzt,k
        call psb_remote_mat(a%rmta,desc_a,a_add,info)
        nz  = a_add%get_nzeros()
        nzt = nz
        call psb_sum(ctxt,nzt)
        if (nzt>0) then
          allocate(ivm, mold=desc_a%v_halo_index%v)
          call psb_cd_reinit(desc_a, info)
        end if
        if (nz > 0) then
          !
          ! Should we check for new indices here?
          !
          call psb_realloc(nz,ila,info)
          call psb_realloc(nz,jla,info)
          call desc_a%indxmap%g2l(a_add%ia(1:nz),ila(1:nz),info,owned=.true.)    
          if (info == 0) call desc_a%indxmap%g2l_ins(a_add%ja(1:nz),jla(1:nz),info)
          !write(0,*) me,name,' Check before insert',a%get_nzeros()
          n_row = desc_a%get_local_rows()
          n_col = desc_a%get_local_cols()
          call a%set_ncols(desc_a%get_local_cols())
          call a%csput(nz,ila,jla,a_add%val,ione,n_row,ione,n_col,info)
          !write(0,*) me,name,' Check after insert',a%get_nzeros(),nz
        end if
        if (nzt > 0) call psb_cdasb(desc_a,info,mold=ivm)

      end block
    end if
    call a%set_ncols(desc_a%get_local_cols())    
    call a%cscnv(info,type=afmt,mold=mold,dupl=dupl_)
  else if (a%is_upd()) then 
    if (a%is_remote_build()) then 
      !write(0,*) me,name,' Size of rmta:',a%rmta%get_nzeros()
      block
        type(psb_ls_coo_sparse_mat) :: a_add
        integer(psb_ipk_), allocatable :: ila(:), jla(:)
        integer(psb_ipk_) :: nz, nzt,k
        call psb_remote_mat(a%rmta,desc_a,a_add,info)
        nz = a_add%get_nzeros()
!!$        write(0,*) me,name,' Nz to be added',nz
        if (nz > 0) then
          !
          ! Should we check for new indices here?
          !
          call psb_realloc(nz,ila,info)
          call psb_realloc(nz,jla,info)
          call desc_a%indxmap%g2l(a_add%ia(1:nz),ila(1:nz),info,owned=.true.)    
          if (info == 0) call desc_a%indxmap%g2l_ins(a_add%ja(1:nz),jla(1:nz),info)
          !write(0,*) me,name,' Check before insert',a%get_nzeros()
          n_row = desc_a%get_local_rows()
          n_col = desc_a%get_local_cols()
          call a%set_ncols(desc_a%get_local_cols())
          call a%csput(nz,ila,jla,a_add%val,ione,n_row,ione,n_col,info)
          !write(0,*) me,name,' Check after insert',a%get_nzeros(),nz
        end if
      end block
    end if
    call a%asb(mold=mold)
  else
    info = psb_err_invalid_mat_state_
    call psb_errpush(info,name)
    goto 9999
    
  end if

  if (.true.) then
    block
      character(len=1024) :: fname
      type(psb_s_coo_sparse_mat) :: acoo
      type(psb_s_csr_sparse_mat), allocatable :: aclip, andclip
      allocate(aclip,andclip)
      call a%a%csclip(acoo,info,jmax=n_row,rscale=.false.,cscale=.false.)
      call aclip%mv_from_coo(acoo,info)
      call a%a%csclip(acoo,info,jmin=n_row+1,jmax=n_col,rscale=.false.,cscale=.false.)
      call andclip%mv_from_coo(acoo,info)
      call move_alloc(aclip,a%ad)
      call move_alloc(andclip,a%and)
      if (.false.) then 
        write(fname,'(a,i2.2,a)') 'adclip_',me,'.mtx'
        open(25,file=fname)
        call a%ad%print(25)
        close(25)
        write(fname,'(a,i2.2,a)') 'andclip_',me,'.mtx'
        open(25,file=fname)
        call a%and%print(25)
        close(25)
        !call andclip%set_cols(n_col)
        write(*,*) me,' ',trim(name),' ad  ',&
             &a%ad%get_nrows(),a%ad%get_ncols(),n_row,n_col
        write(*,*) me,' ',trim(name),' and ',&
             &a%and%get_nrows(),a%and%get_ncols(),n_row,n_col
      end if
    end block
  end if
  if (debug_level >= psb_debug_ext_) then 
    ch_err=a%get_fmt()
    write(debug_unit, *) me,' ',trim(name),':  From SPCNV',&
         & info,' ',ch_err
  end if
  
  if (psb_errstatus_fatal()) then    
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='cscnv')
    goto 9999
  endif

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ctxt,err_act)

  return

end subroutine psb_sspasb
