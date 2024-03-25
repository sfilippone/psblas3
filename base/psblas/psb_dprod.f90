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
! File: psb_dprod.f90
!
! Function: psb_dprod_multivect
!    psb_dprod computes the product of two distributed multivectors,
!
!    prod := ( X ) * ( Y ) or
!    prod := ( X )**C * ( Y )
!
!
! Arguments:
!    x      -  type(psb_d_multivect_type) The input vector containing the entries of sub( X ).
!    y      -  type(psb_d_multivect_type) The input vector containing the entries of sub( Y ).
!    desc_a -  type(psb_desc_type).       The communication descriptor.
!    info   -  integer.                   Return code
!    trans  -  logical(optional)          Whether multivector X is transposed, default: .false.
!    global -  logical(optional)          Whether to perform the global reduce, default: .true.
!
!  Note: from a functional point of view, X and Y are input, but here
!        they are declared INOUT because of the sync() methods.
!
!
function psb_dprod_multivect(x,y,desc_a,info,trans,global) result(res)
   use psb_desc_mod
   use psb_d_base_mat_mod
   use psb_check_mod
   use psb_error_mod
   use psb_penv_mod
   use psb_d_vect_mod
   use psb_d_psblas_mod, psb_protect_name => psb_dprod_multivect
   implicit none
   real(psb_dpk_), allocatable               :: res(:,:)
   type(psb_d_multivect_type), intent(inout) :: x, y
   type(psb_desc_type), intent(in)           :: desc_a
   integer(psb_ipk_), intent(out)            :: info
   logical, intent(in), optional             :: trans
   logical, intent(in), optional             :: global

   ! locals
   type(psb_ctxt_type) :: ctxt
   integer(psb_ipk_) :: np, me, idx, ndm,&
   & err_act, iix, jjx, iiy, jjy, i, j, nr
   integer(psb_lpk_) :: ix, ijx, iy, ijy, m
   logical :: global_, trans_
   character(len=20)      :: name, ch_err

   name='psb_dprod_multivect'
   info=psb_success_
   call psb_erractionsave(err_act)
   if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_ ;    goto 9999
   end if

   ctxt=desc_a%get_context()
   call psb_info(ctxt, me, np)
   if (np == -ione) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
   endif
   if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      call psb_errpush(info,name)
      goto 9999
   endif
   if (.not.allocated(y%v)) then
      info = psb_err_invalid_vect_state_
      call psb_errpush(info,name)
      goto 9999
   endif

   if (present(trans)) then
      trans_ = trans
   else
      trans_ = .false.
   end if

   if (present(global)) then
      global_ = global
   else
      global_ = .true.
   end if

   ix = ione
   ijx = ione

   iy = ione
   ijy = ione

   m = desc_a%get_global_rows()

   ! check vector correctness
   call psb_chkvect(m,x%get_ncols(),x%get_nrows(),ix,ijx,desc_a,info,iix,jjx)
   if (info == psb_success_) &
   & call psb_chkvect(m,y%get_ncols(),y%get_nrows(),iy,ijy,desc_a,info,iiy,jjy)
   if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_chkvect'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
   end if

   if ((iix /= ione).or.(iiy /= ione)) then
      info=psb_err_ix_n1_iy_n1_unsupported_
      call psb_errpush(info,name)
      goto 9999
   end if

   nr = desc_a%get_local_rows()
   if (nr > 0) then
      res = x%prod(nr,y,trans_)
      ! adjust dot_local because overlapped elements are computed more than once
      if (size(desc_a%ovrlap_elem,1)>0) then
         if (x%v%is_dev()) call x%sync()
         if (y%v%is_dev()) call y%sync()
         do j=1,x%get_ncols()
            do i=1,size(desc_a%ovrlap_elem,1)
               idx = desc_a%ovrlap_elem(i,1)
               ndm = desc_a%ovrlap_elem(i,2)
               res(j,:) = res(j,:) - (real(ndm-1)/real(ndm))*(x%v%v(idx,:)*y%v%v(idx,:))
            end do
         end do
      end if
   else
      res = dzero
   end if

   ! compute global sum
   if (global_) call psb_sum(ctxt, res)

   call psb_erractionrestore(err_act)
   return

9999 call psb_error_handler(ctxt,err_act)

   return

end function psb_dprod_multivect
!
! Function: psb_dprod_multivect_a
!    psb_dprod computes the product of two distributed multivectors,
!
!    prod := ( X ) * ( Y ) or
!    prod := ( X )**C * ( Y )
!
!
! Arguments:
!    x      -  type(psb_d_multivect_type) The input vector containing the entries of sub( X ).
!    y      -  real(:,:)                  The input vector containing the entries of sub( Y ).
!    desc_a -  type(psb_desc_type).       The communication descriptor.
!    info   -  integer.                   Return code
!    trans  -  logical(optional)          Whether multivector X is transposed, default: .false.
!    global -  logical(optional)          Whether to perform the global reduce, default: .true.
!
!  Note: from a functional point of view, X and Y are input, but here
!        they are declared INOUT because of the sync() methods.
!
!
function psb_dprod_multivect_a(x,y,desc_a,info,trans,global) result(res)
   use psb_desc_mod
   use psb_d_base_mat_mod
   use psb_check_mod
   use psb_error_mod
   use psb_penv_mod
   use psb_d_vect_mod
   use psb_d_psblas_mod, psb_protect_name => psb_dprod_multivect_a
   implicit none
   real(psb_dpk_), allocatable               :: res(:,:)
   type(psb_d_multivect_type), intent(inout) :: x
   real(psb_dpk_), intent(in)                :: y(:,:)
   type(psb_desc_type), intent(in)           :: desc_a
   integer(psb_ipk_), intent(out)            :: info
   logical, intent(in), optional             :: trans
   logical, intent(in), optional             :: global

   ! locals
   type(psb_ctxt_type) :: ctxt
   integer(psb_ipk_) :: np, me, idx, ndm,&
   & err_act, iix, jjx, iiy, jjy, i, j, nr
   integer(psb_lpk_) :: ix, ijx, iy, ijy, m
   logical :: global_, trans_
   character(len=20)      :: name, ch_err

   name='psb_dprod_multivect'
   info=psb_success_
   call psb_erractionsave(err_act)
   if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_ ;    goto 9999
   end if

   ctxt=desc_a%get_context()
   call psb_info(ctxt, me, np)
   if (np == -ione) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
   endif
   if (.not.allocated(x%v)) then
      info = psb_err_invalid_vect_state_
      call psb_errpush(info,name)
      goto 9999
   endif

   if (present(trans)) then
      trans_ = trans
   else
      trans_ = .false.
   end if

   if (present(global)) then
      global_ = global
   else
      global_ = .true.
   end if

   ix = ione
   ijx = ione

   m = desc_a%get_global_rows()

   ! check vector correctness
   call psb_chkvect(m,x%get_ncols(),x%get_nrows(),ix,ijx,desc_a,info,iix,jjx)
   if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_chkvect'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
   end if

   if ((iix /= ione)) then
      info=psb_err_ix_n1_iy_n1_unsupported_
      call psb_errpush(info,name)
      goto 9999
   end if

   nr = desc_a%get_local_rows()
   if (nr > 0) then
      res = x%prod(nr,y,trans_)
      ! adjust dot_local because overlapped elements are computed more than once
      if (size(desc_a%ovrlap_elem,1)>0) then
         if (x%v%is_dev()) call x%sync()
         do j=1,x%get_ncols()
            do i=1,size(desc_a%ovrlap_elem,1)
               idx = desc_a%ovrlap_elem(i,1)
               ndm = desc_a%ovrlap_elem(i,2)
               res(j,:) = res(j,:) - (real(ndm-1)/real(ndm))*(x%v%v(idx,:)*y(idx,:))
            end do
         end do
      end if
   else
      res = dzero
   end if

   ! compute global sum
   if (global_) call psb_sum(ctxt, res)

   call psb_erractionrestore(err_act)
   return

9999 call psb_error_handler(ctxt,err_act)

   return

end function psb_dprod_multivect_a
!
! Function: psb_dprod_m
!    psb_dprod computes the product of two distributed multivectors,
!
!    prod := ( X ) * ( Y ) or
!    prod := ( X )**C * ( Y )
!
!
! Arguments:
!    x      -  real(:,:)                  The input vector containing the entries of sub( X ).
!    y      -  real(:,:)                  The input vector containing the entries of sub( Y ).
!    desc_a -  type(psb_desc_type).       The communication descriptor.
!    info   -  integer.                   Return code
!    trans  -  logical(optional)          Whether multivector X is transposed, default: .false.
!    global -  logical(optional)          Whether to perform the global reduce, default: .true.
!
!  Note: from a functional point of view, X and Y are input, but here
!        they are declared INOUT because of the sync() methods.
!
!
function psb_dprod_m(x,y,desc_a,info,trans,global) result(res)
   use psb_desc_mod
   use psb_d_base_mat_mod
   use psb_check_mod
   use psb_error_mod
   use psb_penv_mod
   use psb_d_vect_mod
   use psb_d_psblas_mod, psb_protect_name => psb_dprod_m
   implicit none
   real(psb_dpk_), allocatable               :: res(:,:)
   real(psb_dpk_), intent(in)                :: x(:,:), y(:,:)
   type(psb_desc_type), intent(in)           :: desc_a
   integer(psb_ipk_), intent(out)            :: info
   logical, intent(in), optional             :: trans
   logical, intent(in), optional             :: global

   ! locals
   type(psb_ctxt_type) :: ctxt
   integer(psb_ipk_) :: np, me, idx, ndm,&
   & err_act, iix, jjx, iiy, jjy, i, j, nr, x_n, y_n, lda, ldb
   integer(psb_lpk_) :: ix, ijx, iy, ijy, m
   logical :: global_, trans_
   character(len=20)      :: name, ch_err

   name='psb_dprod_multivect'
   info=psb_success_
   call psb_erractionsave(err_act)
   if (psb_errstatus_fatal()) then
      info = psb_err_internal_error_ ;    goto 9999
   end if

   ctxt=desc_a%get_context()
   call psb_info(ctxt, me, np)
   if (np == -ione) then
      info = psb_err_context_error_
      call psb_errpush(info,name)
      goto 9999
   endif

   if (present(trans)) then
      trans_ = trans
   else
      trans_ = .false.
   end if

   if (present(global)) then
      global_ = global
   else
      global_ = .true.
   end if

   ix = ione
   ijx = ione

   m = desc_a%get_global_rows()

   ! check vector correctness
   call psb_chkvect(m,size(x,2),size(x,1),ix,ijx,desc_a,info,iix,jjx)
   if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_chkvect'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
   end if

   if ((iix /= ione)) then
      info=psb_err_ix_n1_iy_n1_unsupported_
      call psb_errpush(info,name)
      goto 9999
   end if

   nr = desc_a%get_local_rows()
   x_n = size(x,2)
   y_n = size(y,2)
   lda = size(x,1)
   if (nr > 0) then

      if (trans_) then
         allocate(res(x_n,y_n))
         res = dzero
         ldb = size(y,1)
         call dgemm('T','N',x_n,y_n,nr,done,x,lda,y,ldb,dzero,res,x_n)
      else
         allocate(res(lda,y_n))
         res = dzero
         ldb = x_n
         call dgemm('N','N',nr,y_n,x_n,done,x,lda,y,ldb,dzero,res,lda)
      end if

      ! adjust dot_local because overlapped elements are computed more than once
      if (size(desc_a%ovrlap_elem,1)>0) then
         do j=1,x_n
            do i=1,size(desc_a%ovrlap_elem,1)
               idx = desc_a%ovrlap_elem(i,1)
               ndm = desc_a%ovrlap_elem(i,2)
               res(j,:) = res(j,:) - (real(ndm-1)/real(ndm))*(x(idx,:)*y(idx,:))
            end do
         end do
      end if
   else
      if (trans_) then
         allocate(res(x_n,y_n))
      else
         allocate(res(lda,y_n))
      end if
      res = dzero
   end if

   ! compute global sum
   if (global_) call psb_sum(ctxt, res)

   call psb_erractionrestore(err_act)
   return

9999 call psb_error_handler(ctxt,err_act)

   return

end function psb_dprod_m
