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
! File: psb_s_par_csr_spspmm.f90
!
! Subroutine: psb_s_par_csr_spspmm
! Version:    real
!
!  This routine computes a  parallel product of two sparse matrices
!
!         C = A * B 
!
!  where all the matrices are stored in CSR. On input and output the matrices
!  are stored with column indices in local numbering, but intermediate quantities
!  are in global numbering because gathering the halo of B to multiply it
!  by A implies a potential enlargement of the support.
!  Also, B may have a column index space different from its row index space,
!  which is obviously the same as the column  space of A. 
! 
!
! Arguments:
!    acsr       -  type(psb_s_csr_sparse_mat), input.     
!                  The sparse matrix structure A 
!    desc_a     -  type(psb_desc_type), input.
!                  The communication descriptor of the column space of A
!    bcsr       -  type(psb_s_csr_sparse_mat), input/output.     
!                  The sparse matrix structure B, gets row-extended on output 
!    ccsr       -  type(psb_s_csr_sparse_mat), output      
!                  The sparse matrix structure C 
!    desc_c     -  type(psb_desc_type), input/output.
!                  The communication descriptor of the column space of B
!
!    info       -  integer, output.
!                  Error code.
!
Subroutine psb_s_par_csr_spspmm(acsr,desc_a,bcsr,ccsr,desc_c,info,data)
  use psb_base_mod, psb_protect_name => psb_s_par_csr_spspmm
  Implicit None

  type(psb_s_csr_sparse_mat),intent(in)    :: acsr
  type(psb_s_csr_sparse_mat),intent(inout) :: bcsr
  type(psb_s_csr_sparse_mat),intent(out)   :: ccsr      
  type(psb_desc_type),intent(in)           :: desc_a
  type(psb_desc_type),intent(inout)        :: desc_c
  integer(psb_ipk_), intent(out)           :: info
  integer(psb_ipk_), intent(in), optional  :: data
  !     ...local scalars....
  integer(psb_ipk_) :: ictxt, np,me
  integer(psb_ipk_) :: ncol, nnz
  type(psb_ls_csr_sparse_mat) :: ltcsr
  type(psb_s_csr_sparse_mat) :: tcsr
  logical           :: update_desc_c
  integer(psb_ipk_) :: debug_level, debug_unit, err_act
  character(len=20) :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name='psb_s_p_csr_spspmm'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': Start'

  update_desc_c = desc_c%is_bld()

  !
  ! This is a bit tricky.
  ! DESC_A is the descriptor of (the columns of) A, and therefore
  ! of the rows of B; the columns of B, in the intended usage, span
  ! a different space for which we have DESC_C. 
  ! We are gathering the halo rows of B to multiply by A;
  ! now, the columns of B would ideally be kept in
  ! global numbering, so that we can call this repeatedly to accumulate
  ! the product of multiple operators, and convert to local numbering
  ! at the last possible moment. However, this would imply calling
  ! the serial SPSPMM with a matrix B with the GLOBAL number of columns
  ! and this could be very expensive in memory. The solution is to keep B
  ! in local numbering, so that only columns really appearing count, but to
  ! expand the descriptor when gathering the halo, because by performing
  ! the products we are extending the support of the operator; hence
  ! this routine is intended to be called with a temporary descriptor
  ! DESC_C which is in the BUILD state, to allow for such expansion
  ! across multiple products. 
  ! The caller will at some later point finalize the descriptor DESC_C. 
  ! 

  ncol = desc_a%get_local_cols()
  call psb_sphalo(bcsr,desc_a,ltcsr,info,&
       & colcnv=.true.,rowscale=.true.,outcol_glob=.true.,col_desc=desc_c,data=data)
  nnz    = ltcsr%get_nzeros()
  if (update_desc_c) then 
    call desc_c%indxmap%g2lip_ins(ltcsr%ja(1:nnz),info)
  else
    call desc_c%indxmap%g2lip(ltcsr%ja(1:nnz),info)
  end if
  call ltcsr%mv_to_ifmt(tcsr,info)
  if (info == psb_success_) call psb_rwextd(ncol,bcsr,info,b=tcsr)      
  if (info == psb_success_) call tcsr%free()
  if(info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Extend am3')
    goto 9999
  end if
  call bcsr%set_ncols(desc_c%get_local_cols())


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting spspmm 3'
  if (debug_level >= psb_debug_outer_) write(debug_unit,*) me,' ',trim(name),&
       & 'starting spspmm ',acsr%get_nrows(),acsr%get_ncols(),bcsr%get_nrows(),bcsr%get_ncols()
  call psb_spspmm(acsr,bcsr,ccsr,info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

End Subroutine psb_s_par_csr_spspmm

Subroutine psb_ls_par_csr_spspmm(acsr,desc_a,bcsr,ccsr,desc_c,info,data)
  use psb_base_mod, psb_protect_name => psb_ls_par_csr_spspmm
  Implicit None

  type(psb_ls_csr_sparse_mat),intent(in)    :: acsr
  type(psb_ls_csr_sparse_mat),intent(inout) :: bcsr
  type(psb_ls_csr_sparse_mat),intent(out)   :: ccsr      
  type(psb_desc_type),intent(in)           :: desc_a
  type(psb_desc_type),intent(inout)        :: desc_c
  integer(psb_ipk_), intent(out)           :: info
  integer(psb_ipk_), intent(in), optional  :: data
  !     ...local scalars....
  integer(psb_ipk_) :: ictxt, np,me
  integer(psb_lpk_) :: nacol, nccol, nnz
  type(psb_ls_csr_sparse_mat) :: tcsr1
  logical           :: update_desc_c
  integer(psb_ipk_) :: debug_level, debug_unit, err_act
  character(len=20) :: name, ch_err

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  name='psb_ls_p_csr_spspmm'
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ictxt = desc_a%get_context()

  call psb_info(ictxt, me, np)

  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),': Start'

  update_desc_c = desc_c%is_bld()

  !
  ! This is a bit tricky.
  ! DESC_A is the descriptor of (the columns of) A, and therefore
  ! of the rows of B; the columns of B, in the intended usage, span
  ! a different space for which we have DESC_C. 
  ! We are gathering the halo rows of B to multiply by A;
  ! now, the columns of B would ideally be kept in
  ! global numbering, so that we can call this repeatedly to accumulate
  ! the product of multiple operators, and convert to local numbering
  ! at the last possible moment. However, this would imply calling
  ! the serial SPSPMM with a matrix B with the GLOBAL number of columns
  ! and this could be very expensive in memory. The solution is to keep B
  ! in local numbering, so that only columns really appearing count, but to
  ! expand the descriptor when gathering the halo, because by performing
  ! the products we are extending the support of the operator; hence
  ! this routine is intended to be called with a temporary descriptor
  ! DESC_C which is in the BUILD state, to allow for such expansion
  ! across multiple products. 
  ! The caller will at some later point finalize the descriptor DESC_C. 
  ! 

  nacol = desc_a%get_local_cols()
  call psb_sphalo(bcsr,desc_a,tcsr1,info,&
       & colcnv=.true.,rowscale=.true.,outcol_glob=.true.,col_desc=desc_c,data=data)
  nnz    = tcsr1%get_nzeros()
  if (update_desc_c) then 
    call desc_c%indxmap%g2lip_ins(tcsr1%ja(1:nnz),info)
  else
    call desc_c%indxmap%g2lip(tcsr1%ja(1:nnz),info)
  end if
  if (info == psb_success_) call psb_rwextd(nacol,bcsr,info,b=tcsr1)      
  if (info == psb_success_) call tcsr1%free()
  if(info /= psb_success_) then
    call psb_errpush(psb_err_internal_error_,name,a_err='Extend am3')
    goto 9999
  end if
  nccol = desc_c%get_local_cols()
  call bcsr%set_ncols(nccol)


  if (debug_level >= psb_debug_outer_) &
       & write(debug_unit,*) me,' ',trim(name),&
       & 'starting spspmm 3'
  if (debug_level >= psb_debug_outer_) write(debug_unit,*) me,' ',trim(name),&
       & 'starting spspmm ',acsr%get_nrows(),acsr%get_ncols(),bcsr%get_nrows(),bcsr%get_ncols()
  call psb_spspmm(acsr,bcsr,ccsr,info)

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(ictxt,err_act)

  return

End Subroutine psb_ls_par_csr_spspmm
