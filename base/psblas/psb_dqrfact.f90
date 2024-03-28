!
! Subroutine: psb_dqrfact
!    Computes QR factorization of a multivector
!
! Arguments:
!    x      -  type(psb_d_multivect_type) The input multivector containing the entries of X
!    desc_a -  type(psb_desc_type)        The communication descriptor.
!    info   -  integer                    Return code
!
!  Note: from a functional point of view, X is input, but here
!        it's declared INOUT because of the sync() methods.
!
function psb_dqrfact(x, desc_a, info) result(res)
   use psb_base_mod, psb_protect_name => psb_dqrfact
   implicit none
   real(psb_dpk_), allocatable                :: res(:,:)
   type(psb_d_multivect_type), intent(inout)  :: x
   type(psb_desc_type), intent(in)            :: desc_a
   integer(psb_ipk_), intent(out)             :: info

   ! locals
   type(psb_ctxt_type) :: ctxt
   integer(psb_ipk_) :: np, me, err_act, iix, jjx, i
   integer(psb_lpk_) :: ix, ijx, m, n
   character(len=20) :: name, ch_err
   real(psb_dpk_), allocatable :: temp(:,:)
   type(psb_d_multivect_type) :: qr_temp

   name='psb_dgqrfact'
   if (psb_errstatus_fatal()) return
   info=psb_success_
   call psb_erractionsave(err_act)

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

   ix = ione
   ijx = ione

   m = desc_a%get_global_rows()
   n = x%get_ncols()

   call psb_chkvect(m,n,x%get_nrows(),ix,ijx,desc_a,info,iix,jjx)
   if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_chkvect'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
   end if

   if (iix /= ione) then
      info=psb_err_ix_n1_iy_n1_unsupported_
      call psb_errpush(info,name)
   end if

   call psb_gather(temp,x,desc_a,info,root=psb_root_)

   if (me == psb_root_) then
      call qr_temp%bld(temp)
      res = qr_temp%qr_fact(info)
      temp = qr_temp%get_vect()
      call psb_bcast(ctxt,res)
   else
      allocate(res(n,n))
      call psb_bcast(ctxt,res)
   end if

   call psb_scatter(temp,x,desc_a,info,root=psb_root_)

   call psb_erractionrestore(err_act)
   return

9999 call psb_error_handler(ctxt,err_act)

   return
end function psb_dqrfact
