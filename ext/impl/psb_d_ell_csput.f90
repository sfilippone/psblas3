!                Parallel Sparse BLAS   GPU plugin 
!      (C) Copyright 2013
!  
!                         Salvatore Filippone
!                         Alessandro Fanfarillo
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
  

subroutine psb_d_ell_csput_a(nz,ia,ja,val,a,imin,imax,jmin,jmax,info) 
  
  use psb_base_mod
  use psb_d_ell_mat_mod, psb_protect_name => psb_d_ell_csput_a
  implicit none 

  class(psb_d_ell_sparse_mat), intent(inout) :: a
  real(psb_dpk_), intent(in)               :: val(:)
  integer(psb_ipk_), intent(in)             :: nz, ia(:), ja(:), imin,imax,jmin,jmax
  integer(psb_ipk_), intent(out)            :: info


  integer(psb_ipk_)  :: err_act
  character(len=20)  :: name='d_ell_csput_a'
  logical, parameter :: debug=.false.
  integer(psb_ipk_)  :: nza, i,j,k, nzl, isza, int_err(5), debug_level, debug_unit


  call psb_erractionsave(err_act)
  info = psb_success_
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  if (nz <= 0) then 
    info = psb_err_iarg_neg_
    int_err(1)=1
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if
  if (size(ia) < nz) then 
    info = psb_err_input_asize_invalid_i_
    int_err(1)=2
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  if (size(ja) < nz) then 
    info = psb_err_input_asize_invalid_i_
    int_err(1)=3
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if
  if (size(val) < nz) then 
    info = psb_err_input_asize_invalid_i_
    int_err(1)=4
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  if (nz == 0) return

  nza  = a%get_nzeros()

  if (a%is_bld()) then 
    ! Build phase should only ever be in COO
    info = psb_err_invalid_mat_state_

  else  if (a%is_upd()) then 
    if (a%is_dev()) call a%sync()
    call  psb_d_ell_srch_upd(nz,ia,ja,val,a,&
         & imin,imax,jmin,jmax,info)

    if (info < 0) then 
      info = psb_err_internal_error_
    else if (info > 0) then 
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*) trim(name),&
           & ': Discarded entries not  belonging to us.'                    
      info = psb_success_
    end if
    call a%set_host()
  else 
    ! State is wrong.
    info = psb_err_invalid_mat_state_
  end if
  if (info /= psb_success_) then
    call psb_errpush(info,name)
    goto 9999
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return


contains

  subroutine psb_d_ell_srch_upd(nz,ia,ja,val,a,&
       & imin,imax,jmin,jmax,info)

    implicit none 

    class(psb_d_ell_sparse_mat), intent(inout) :: a
    integer(psb_ipk_), intent(in)  :: nz, imin,imax,jmin,jmax
    integer(psb_ipk_), intent(in)  :: ia(:),ja(:)
    real(psb_dpk_), intent(in)    :: val(:)
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_)    :: i,ir,ic, ilr, ilc, ip, &
         & i1,i2,nr,nc,nnz,dupl
    integer(psb_ipk_)    :: debug_level, debug_unit
    character(len=20)    :: name='d_ell_srch_upd'

    info = psb_success_
    debug_unit  = psb_get_debug_unit()
    debug_level = psb_get_debug_level()

    dupl = a%get_dupl()

    if (.not.a%is_sorted()) then 
      info = -4
      return
    end if

    ilr = -1 
    ilc = -1 
    nnz = a%get_nzeros()
    nr  = a%get_nrows()
    nc  = a%get_ncols()

    select case(dupl)
    case(psb_dupl_ovwrt_,psb_dupl_err_)
      ! Overwrite.
      ! Cannot test for error, should have been caught earlier.

      ilr = -1 
      ilc = -1 
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 

        if ((ir > 0).and.(ir <= nr)) then 

          nc = a%irn(ir)
          ip = psb_bsrch(ic,nc,a%ja(ir,1:nc))    
          if (ip>0) then 
            a%val(ir,ip) = val(i)
          else
            info = max(info,3)
          end if
        else
          info = max(info,2)
        end if

      end do

    case(psb_dupl_add_)
      ! Add
      ilr = -1 
      ilc = -1 
      do i=1, nz
        ir = ia(i)
        ic = ja(i) 
        if ((ir > 0).and.(ir <= nr)) then 
          nc = a%irn(ir)
          ip = psb_bsrch(ic,nc,a%ja(ir,1:nc))    
          if (ip>0) then 
            a%val(ir,ip) = a%val(ir,ip) + val(i)
          else
            info = max(info,3)
          end if
        else
          info = max(info,2)
        end if
      end do

    case default
      info = -3
      if (debug_level >= psb_debug_serial_) &
           & write(debug_unit,*) trim(name),&
           & ': Duplicate handling: ',dupl
    end select

  end subroutine psb_d_ell_srch_upd
end subroutine psb_d_ell_csput_a
