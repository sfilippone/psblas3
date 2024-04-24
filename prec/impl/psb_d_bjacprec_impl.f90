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

subroutine psb_d_bjac_dump(prec,info,prefix,head)
  use psb_base_mod
  use psb_d_bjacprec, psb_protect_name => psb_d_bjac_dump
  implicit none
  class(psb_d_bjac_prec_type), intent(in) :: prec
  integer(psb_ipk_), intent(out)                    :: info
  character(len=*), intent(in), optional  :: prefix,head
  integer(psb_ipk_) :: i, j, il1, iln, lname, lev
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: iam, np
  character(len=80)  :: prefix_
  character(len=120) :: fname ! len should be at least 20 more than

  !  len of prefix_

  info = 0
  ctxt = prec%get_ctxt()
  call psb_info(ctxt,iam,np)

  if (present(prefix)) then
    prefix_ = trim(prefix(1:min(len(prefix),len(prefix_))))
  else
    prefix_ = "dump_fact_d"
  end if

  lname = len_trim(prefix_)
  fname = trim(prefix_)
  write(fname(lname+1:lname+5),'(a,i3.3)') '_p',iam
  lname = lname + 5
  write(fname(lname+1:),'(a)')'_lower.mtx'
  if (prec%av(psb_l_pr_)%is_asb())  &
       & call prec%av(psb_l_pr_)%print(fname,head=head)
  write(fname(lname+1:),'(a,a)')'_diag.mtx'
  if (allocated(prec%dv)) &
       & call psb_geprt(fname,prec%dv%v%v,head=head)
  write(fname(lname+1:),'(a)')'_upper.mtx'
  if (prec%av(psb_u_pr_)%is_asb()) &
       & call prec%av(psb_u_pr_)%print(fname,head=head)

end subroutine psb_d_bjac_dump

subroutine psb_d_bjac_apply_vect(alpha,prec,x,beta,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_d_bjacprec, psb_protect_name => psb_d_bjac_apply_vect
  implicit none
  type(psb_desc_type),intent(in)    :: desc_data
  class(psb_d_bjac_prec_type), intent(inout)  :: prec
  real(psb_dpk_),intent(in)         :: alpha,beta
  type(psb_d_vect_type),intent(inout)   :: x
  type(psb_d_vect_type),intent(inout)   :: y
  integer(psb_ipk_), intent(out)              :: info
  character(len=1), optional        :: trans
  real(psb_dpk_),intent(inout), optional, target :: work(:)

  ! Local variables
  integer(psb_ipk_) :: n_row,n_col
  real(psb_dpk_), pointer :: ww(:), aux(:)
  type(psb_d_vect_type) :: wv, wv1
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me
  integer(psb_ipk_) :: err_act, ierr(5)
  integer(psb_ipk_) :: debug_level, debug_unit
  logical            :: do_alloc_wrk
  character          :: trans_
  character(len=20)  :: name='d_bjac_prec_apply'
  character(len=20)  :: ch_err

  info = psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ctxt       = desc_data%get_context()
  call psb_info(ctxt, me, np)


  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N','T','C')
    ! Ok
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select


  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()

  if (x%get_nrows() < n_row) then
    info = 36; ierr(1) = 2; ierr(2) = n_row;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (y%get_nrows() < n_row) then
    info = 36; ierr(1) = 3; ierr(2) = n_row;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (.not.allocated(prec%dv)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: D")
    goto 9999
  end if
  if (prec%dv%get_nrows() < n_row) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: D")
    goto 9999
  end if


  if (n_col <= size(work)) then
    ww => work(1:n_col)
    if ((4*n_col+n_col) <= size(work)) then
      aux => work(n_col+1:)
    else
      allocate(aux(4*n_col),stat=info)

    endif
  else
    allocate(ww(n_col),aux(4*n_col),stat=info)
  endif

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999
  end if

  do_alloc_wrk = .not.prec%is_allocated_wrk()
  if (do_alloc_wrk) call prec%allocate_wrk(info,vmold=x%v)

  associate (wv => prec%wrk(1), wv1 => prec%wrk(2))

    select case(prec%iprcparm(psb_f_type_))
    case(psb_f_ilu_n_,psb_f_ilu_k_,psb_f_ilu_t_)

      select case(trans_)
      case('N')
        call psb_spsm(done,prec%av(psb_l_pr_),x,dzero,wv,desc_data,info,&
             & trans=trans_,scale='L',diag=prec%dv,choice=psb_none_,work=aux)
        if(info == psb_success_) call psb_spsm(alpha,prec%av(psb_u_pr_),wv,&
             & beta,y,desc_data,info,&
             & trans=trans_,scale='U',choice=psb_none_, work=aux)

      case('T')
        call psb_spsm(done,prec%av(psb_u_pr_),x,dzero,wv,desc_data,info,&
             & trans=trans_,scale='L',diag=prec%dv,choice=psb_none_, work=aux)
        if(info == psb_success_)  call psb_spsm(alpha,prec%av(psb_l_pr_),wv,&
             & beta,y,desc_data,info,&
             & trans=trans_,scale='U',choice=psb_none_,work=aux)

      case('C')

        call psb_spsm(done,prec%av(psb_u_pr_),x,dzero,wv,desc_data,info,&
             & trans=trans_,scale='U',choice=psb_none_, work=aux)

        call wv1%mlt(done,prec%dv,wv,dzero,info,conjgx=trans_)

        if(info == psb_success_)  call psb_spsm(alpha,prec%av(psb_l_pr_),wv1,&
             & beta,y,desc_data,info,&
             & trans=trans_,scale='U',choice=psb_none_,work=aux)

      end select
      if (info /= psb_success_) then
        ch_err="psb_spsm"
        goto 9999
      end if

    case(psb_f_ainv_,psb_f_invt_,psb_f_invk_)
      ! Application of approximate inverse preconditioner, just some spmm

      select case(trans_)
      case('N')
        call psb_spmm(done,prec%av(psb_l_pr_),x,dzero,wv,desc_data,info,&
             & trans=trans_,work=aux,doswap=.false.)

        if (info == psb_success_) call wv1%mlt(done,prec%dv,wv,dzero,info)
        if(info == psb_success_) &
             & call psb_spmm(alpha,prec%av(psb_u_pr_),wv1,&
             & beta,y,desc_data,info, trans=trans_, work=aux,doswap=.false.)

       case('T','C')
         call psb_spmm(done,prec%av(psb_l_pr_),x,dzero,wv,desc_data,info,&
              & trans=trans_,work=aux,doswap=.false.)
         if (info == psb_success_) call wv1%mlt(done,prec%dv,wv,dzero,info)
         if (info == psb_success_) &
              & call psb_spmm(alpha,prec%av(psb_u_pr_),wv1, &
              & beta,y,desc_data,info,trans=trans_,work=aux,doswap=.false.)

      end select
      if (info /= psb_success_) then
        ch_err="psb_spsm"
        goto 9999
      end if

    case default
      info = psb_err_internal_error_
      call psb_errpush(info,name,a_err='Invalid factorization')
      goto 9999
    end select

  end associate

  call psb_halo(y,desc_data,info,data=psb_comm_mov_)

  if (do_alloc_wrk) call prec%free_wrk(info)
  if (n_col <= size(work)) then
    if ((4*n_col+n_col) <= size(work)) then
    else
      deallocate(aux)
    endif
  else
    deallocate(ww,aux)
  endif


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name,i_err=ierr,a_err=ch_err)
  call psb_error_handler(err_act)
  return

end subroutine psb_d_bjac_apply_vect

subroutine psb_d_bjac_apply(alpha,prec,x,beta,y,desc_data,info,trans,work)
  use psb_base_mod
  use psb_d_bjacprec, psb_protect_name => psb_d_bjac_apply
  implicit none
  type(psb_desc_type),intent(in)    :: desc_data
  class(psb_d_bjac_prec_type), intent(inout)  :: prec
  real(psb_dpk_),intent(in)         :: alpha,beta
  real(psb_dpk_),intent(inout)      :: x(:)
  real(psb_dpk_),intent(inout)      :: y(:)
  integer(psb_ipk_), intent(out)              :: info
  character(len=1), optional        :: trans
  real(psb_dpk_),intent(inout), optional, target :: work(:)

  ! Local variables
  integer(psb_ipk_) :: n_row,n_col
  real(psb_dpk_), pointer :: ww(:), aux(:)
  type(psb_d_vect_type) :: tx,ty
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me
  integer(psb_ipk_) :: err_act, ierr(5)
  integer(psb_ipk_) :: debug_level, debug_unit
  character          :: trans_
  character(len=20)  :: name='d_bjac_prec_apply'
  character(len=20)  :: ch_err

  info = psb_success_
  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()
  ctxt       = desc_data%get_context()
  call psb_info(ctxt, me, np)


  trans_ = psb_toupper(trans)
  select case(trans_)
  case('N','T','C')
    ! Ok
  case default
    call psb_errpush(psb_err_iarg_invalid_i_,name)
    goto 9999
  end select


  n_row = desc_data%get_local_rows()
  n_col = desc_data%get_local_cols()

  if (size(x) < n_row) then
    info = 36; ierr(1) = 2; ierr(2) = n_row;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (size(y) < n_row) then
    info = 36; ierr(1) = 3; ierr(2) = n_row;
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  end if
  if (.not.allocated(prec%dv)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: D")
    goto 9999
  end if
  if (prec%dv%get_nrows() < n_row) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner: D")
    goto 9999
  end if


  if (n_col <= size(work)) then
    ww => work(1:n_col)
    if ((4*n_col+n_col) <= size(work)) then
      aux => work(n_col+1:)
    else
      allocate(aux(4*n_col),stat=info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
        goto 9999
      end if

    endif
  else
    allocate(ww(n_col),aux(4*n_col),stat=info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999
    end if
  endif


  select case(prec%iprcparm(psb_f_type_))
  case(psb_f_ilu_n_,psb_f_ilu_k_,psb_f_ilu_t_)

    select case(trans_)
    case('N')
      call psb_spsm(done,prec%av(psb_l_pr_),x,dzero,ww,desc_data,info,&
           & trans=trans_,scale='L',diag=prec%dv%v%v,choice=psb_none_,work=aux)
      if(info == psb_success_) call psb_spsm(alpha,prec%av(psb_u_pr_),ww,&
           & beta,y,desc_data,info,&
           & trans=trans_,scale='U',choice=psb_none_, work=aux)

    case('T')
      call psb_spsm(done,prec%av(psb_u_pr_),x,dzero,ww,desc_data,info,&
           & trans=trans_,scale='L',diag=prec%dv%v%v,choice=psb_none_, work=aux)
      if(info == psb_success_)  call psb_spsm(alpha,prec%av(psb_l_pr_),ww,&
           & beta,y,desc_data,info,&
           & trans=trans_,scale='U',choice=psb_none_,work=aux)

    case('C')

      call psb_spsm(done,prec%av(psb_u_pr_),x,dzero,ww,desc_data,info,&
           & trans=trans_,scale='U',choice=psb_none_, work=aux)
      ww(1:n_row) = ww(1:n_row)*(prec%dv%v%v(1:n_row))
      if(info == psb_success_)  call psb_spsm(alpha,prec%av(psb_l_pr_),ww,&
           & beta,y,desc_data,info,&
           & trans=trans_,scale='U',choice=psb_none_,work=aux)

    end select
    if (info /= psb_success_) then
      ch_err="psb_spsm"
      goto 9999
    end if

  case(psb_f_ainv_,psb_f_invt_,psb_f_invk_)
    ! Application of approximate inverse preconditioner, just some spmm

    select case(trans_)

      case('N')
        call psb_spmm(done,prec%av(psb_l_pr_),x,dzero,ww,desc_data,info,&
            & trans=trans_,work=aux,doswap=.false.)
        ww(1:n_row) = ww(1:n_row) * prec%dv%v%v(1:n_row)
        if (info == psb_success_) &
            & call psb_spmm(alpha,prec%av(psb_u_pr_),ww,beta,y,desc_data,info,&
            & trans=trans_,work=aux,doswap=.false.)

      case('T')
        call psb_spmm(done,prec%av(psb_u_pr_),x,dzero,ww,desc_data,info,&
            & trans=trans_,work=aux,doswap=.false.)
        ww(1:n_row) = ww(1:n_row) * prec%dv%v%v(1:n_row)
        if (info == psb_success_) &
            & call psb_spmm(alpha,prec%av(psb_l_pr_),ww,beta,y,desc_data,info,&
            & trans=trans_,work=aux,doswap=.false.)

      case('C')
        call psb_spmm(done,prec%av(psb_u_pr_),x,dzero,ww,desc_data,info,&
             & trans=trans_,work=aux,doswap=.false.)
        ww(1:n_row) = ww(1:n_row) * (prec%dv%v%v(1:n_row))
        if (info == psb_success_) &
             & call psb_spmm(alpha,prec%av(psb_l_pr_),ww,beta,y,desc_data,info,&
             & trans=trans_,work=aux,doswap=.false.)

    end select


  case default
    info = psb_err_internal_error_
    call psb_errpush(info,name,a_err='Invalid factorization')
    goto 9999
  end select

  call psb_halo(y,desc_data,info,data=psb_comm_mov_)

  if (n_col <= size(work)) then
    if ((4*n_col+n_col) <= size(work)) then
    else
      deallocate(aux)
    endif
  else
    deallocate(ww,aux)
  endif


  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name,i_err=ierr,a_err=ch_err)
  call psb_error_handler(err_act)
  return

end subroutine psb_d_bjac_apply

subroutine psb_d_bjac_precinit(prec,info)
  use psb_base_mod
  use psb_d_bjacprec, psb_protect_name => psb_d_bjac_precinit
  Implicit None

  class(psb_d_bjac_prec_type),intent(inout) :: prec
  integer(psb_ipk_), intent(out)                     :: info
  integer(psb_ipk_) :: err_act, nrow
  character(len=20)  :: name='d_bjac_precinit'

  call psb_erractionsave(err_act)

  info = psb_success_
  call psb_realloc(psb_ifpsz,prec%iprcparm,info)
  if (info /= psb_success_) then
    info = psb_err_alloc_dealloc_
    call psb_Errpush(info,name)
    goto 9999
  end if
  call psb_realloc(psb_rfpsz,prec%rprcparm,info)
  if (info /= psb_success_) then
    info = psb_err_alloc_dealloc_
    call psb_Errpush(info,name)
    goto 9999
  end if

  prec%iprcparm(:)                = 0
  prec%iprcparm(psb_p_type_)      = psb_bjac_
  prec%iprcparm(psb_f_type_)      = psb_f_ilu_n_
  prec%iprcparm(psb_ilu_fill_in_) = 0
  prec%iprcparm(psb_ilu_ialg_)    = psb_ilu_n_
  prec%iprcparm(psb_ilu_scale_)   = psb_ilu_scale_none_
  prec%iprcparm(psb_inv_fillin_)  = 0
  prec%iprcparm(psb_ainv_alg_)    = psb_ainv_llk_


  prec%rprcparm(:)                = 0
  prec%rprcparm(psb_fact_eps_)    = 1E-1_psb_dpk_
  prec%rprcparm(psb_inv_thresh_)  = 1E-1_psb_dpk_


  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return

end subroutine psb_d_bjac_precinit


subroutine psb_d_bjac_precbld(a,desc_a,prec,info,amold,vmold,imold)

  use psb_base_mod
  use psb_d_ilu_fact_mod
  use psb_d_ainv_fact_mod
  use psb_d_bjacprec, psb_protect_name => psb_d_bjac_precbld
  Implicit None

  type(psb_dspmat_type), intent(inout), target :: a
  type(psb_desc_type), intent(inout), target   :: desc_a
  class(psb_d_bjac_prec_type),intent(inout) :: prec
  integer(psb_ipk_), intent(out)                      :: info
  class(psb_d_base_sparse_mat), intent(in), optional :: amold
  class(psb_d_base_vect_type), intent(in), optional  :: vmold
  class(psb_i_base_vect_type), intent(in), optional  :: imold

  !     .. Local Scalars ..
  integer(psb_ipk_) ::    i, m, ialg, fill_in, iscale, inv_fill, iinvalg
  integer(psb_ipk_) ::    ierr(5)
  character ::        trans, unitd
  type(psb_dspmat_type), allocatable :: lf, uf
  real(psb_dpk_), allocatable :: dd(:)
  real(psb_dpk_) :: fact_eps, inv_thresh
  integer(psb_ipk_) :: nztota,  err_act, n_row, nrow_a,n_col, nhalo
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np,me
  character(len=20)  :: name='d_bjac_precbld'
  character(len=20)  :: ch_err


  info = psb_success_

  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if

  ctxt=desc_a%get_ctxt()
  call prec%set_ctxt(ctxt)
  call psb_info(ctxt, me, np)

  m = a%get_nrows()
  if (m < 0) then
    info = psb_err_iarg_neg_
    ierr(1) = 1;  ierr(2) = m
    call psb_errpush(info,name,i_err=ierr)
    goto 9999
  endif
  trans = 'N'
  unitd = 'U'

  ! We check if all the information contained in the preconditioner structure
  ! are meaningful, otherwise we give an error and get out of the build
  ! procedure
  ialg    = prec%iprcparm(psb_ilu_ialg_)  ! Integer for ILU type algorithm
  iinvalg = prec%iprcparm(psb_ainv_alg_)  ! Integer for AINV type algorithm
  iscale = prec%iprcparm(psb_ilu_scale_)  ! Integer for scaling of matrix
  fact_eps = prec%rprcparm(psb_fact_eps_) ! Drop-tolerance for factorization
  inv_thresh = prec%rprcparm(psb_inv_thresh_) ! Drop-tolerance for inverse
  fill_in = prec%iprcparm(psb_ilu_fill_in_) ! Fill-In for factorization
  inv_fill = prec%iprcparm(psb_inv_fillin_) ! Fill-In for inverse factorization

  ! Check if the type of scaling is known, pay attention that not all the
  ! scalings make sense for all the factorization, if something that does not
  ! make sense is required the factorization routine will fail in an
  ! unnrecoverable way.
  if ((iscale == psb_ilu_scale_none_).or.&
      (iscale == psb_ilu_scale_maxval_).or.&
      (iscale == psb_ilu_scale_diag_).or.&
      (iscale == psb_ilu_scale_arwsum_).or.&
      (iscale == psb_ilu_scale_aclsum_).or.&
      (iscale == psb_ilu_scale_arcsum_)) then
  ! Do nothing: admissible request
  else
    info=psb_err_from_subroutine_
    ch_err='psb_ilu_scale_'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  select case (prec%iprcparm(psb_f_type_))
    case (psb_f_ainv_)
      ! Check if the variant for the AINV is known to the library
      select case (iinvalg)
      case(psb_ainv_llk_,psb_ainv_s_llk_,psb_ainv_s_ft_llk_,psb_ainv_llk_noth_,&
        & psb_ainv_mlk_)
        ! Do nothing these are okay
      case default
        info=psb_err_from_subroutine_
        ch_err='psb_ainv_alg_'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end select ! AINV Variant
      ! Check if the drop-tolerance make sense
      if( inv_thresh > 1) then
        info=psb_err_from_subroutine_
        ch_err='psb_inv_thresh_'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    case (psb_f_ilu_t_)
      if (fact_eps > 1) then
        ! Check if the drop-tolerance make sense
        info=psb_err_from_subroutine_
        ch_err='psb_fact_eps_'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    case (psb_f_invt_)
      ! Check both tolerances
      if (fact_eps > 1) then
        info=psb_err_from_subroutine_
        ch_err='psb_fact_eps_'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
      if( inv_thresh > 1) then
        info=psb_err_from_subroutine_
        ch_err='psb_inv_thresh_'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    case default
  end select


  ! Checks relative to the fill-in parameters
  if (prec%iprcparm(psb_f_type_) == psb_f_ilu_n_) then
    if(fill_in < 0) then
      info=psb_err_from_subroutine_
      ch_err='psb_ilu_fill_in_'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    else if (fill_in == 0) then
      ! If the requested level of fill is equal to zero, we default to the
      ! specialized ILU(0) routine
      prec%iprcparm(psb_f_type_) = psb_f_ilu_n_
    end if
  end if
  ! If no limit on the fill_in is required we allow every fill, this is needed
  ! since this quantity is used to allocate the auxiliary vectors for the
  ! factorization
  if (inv_fill <= 0) inv_fill = m

  ! Select on the type of factorization to be used
  select case(prec%iprcparm(psb_f_type_))

  case(psb_f_ilu_n_)
  ! ILU(0) Factorization: the total number of nonzeros of the factorized matrix
  ! is equal to the one of the input matrix

    if (allocated(prec%av)) then
      if (size(prec%av) < psb_bp_ilu_avsz) then
        do i=1,size(prec%av)
          call prec%av(i)%free()
        enddo
        deallocate(prec%av,stat=info)
      endif
    end if
    if (.not.allocated(prec%av)) then
      allocate(prec%av(psb_max_avsz),stat=info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      end if
    endif

    nrow_a = desc_a%get_local_rows()
    nztota = a%get_nzeros()

    n_col  = desc_a%get_local_cols()
    nhalo  = n_col-nrow_a
    n_row  = nrow_a

    allocate(lf,uf,stat=info)
    if (info == psb_success_) call lf%allocate(n_row,n_row,nztota)
    if (info == psb_success_) call uf%allocate(n_row,n_row,nztota)

    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    allocate(dd(n_row),stat=info)
    if (info == psb_success_) then
      allocate(prec%dv, stat=info)
      if (info == 0) then
        if (present(vmold)) then
          allocate(prec%dv%v,mold=vmold,stat=info)
        else
          allocate(psb_d_base_vect_type :: prec%dv%v,stat=info)
        end if
      end if
    end if

    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999
    endif
    ! This is where we have no renumbering, thus no need
    call psb_ilu0_fact(ialg,a,lf,uf,dd,info)

    if(info == psb_success_) then
      call prec%av(psb_l_pr_)%mv_from(lf%a)
      call prec%av(psb_u_pr_)%mv_from(uf%a)
      call prec%av(psb_l_pr_)%set_asb()
      call prec%av(psb_u_pr_)%set_asb()
      call prec%av(psb_l_pr_)%trim()
      call prec%av(psb_u_pr_)%trim()
      call prec%dv%bld(dd)
      ! call move_alloc(dd,prec%d)
    else
      info=psb_err_from_subroutine_
      ch_err='psb_ilu0_fact'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case(psb_f_ilu_k_)
  ! ILU(N) Incomplete LU-factorization with N levels of fill-in. Depending on
  ! the type of the variant of the algorithm the may be forgotten or added to
  ! the diagonal (MILU)

  if (allocated(prec%av)) then
    if (size(prec%av) < psb_bp_ilu_avsz) then
      do i=1,size(prec%av)
        call prec%av(i)%free()
      enddo
      deallocate(prec%av,stat=info)
    endif
  end if
  if (.not.allocated(prec%av)) then
    allocate(prec%av(psb_max_avsz),stat=info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_alloc_dealloc_,name)
      goto 9999
    end if
  endif

  nrow_a = desc_a%get_local_rows()
  nztota = a%get_nzeros()

  n_col  = desc_a%get_local_cols()
  nhalo  = n_col-nrow_a
  n_row  = nrow_a

  allocate(lf,uf,stat=info)
  if (info == psb_success_) call lf%allocate(n_row,n_row,nztota)
  if (info == psb_success_) call uf%allocate(n_row,n_row,nztota)

  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_sp_all'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  allocate(dd(n_row),stat=info)
  if (info == psb_success_) then
    allocate(prec%dv, stat=info)
    if (info == 0) then
      if (present(vmold)) then
        allocate(prec%dv%v,mold=vmold,stat=info)
      else
        allocate(psb_d_base_vect_type :: prec%dv%v,stat=info)
      end if
    end if
  end if

  if (info /= psb_success_) then
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999
  endif
  ! This is where we have no renumbering, thus no need
  call psb_iluk_fact(fill_in,ialg,a,lf,uf,dd,info)

  if(info == psb_success_) then
    call prec%av(psb_l_pr_)%mv_from(lf%a)
    call prec%av(psb_u_pr_)%mv_from(uf%a)
    call prec%av(psb_l_pr_)%set_asb()
    call prec%av(psb_u_pr_)%set_asb()
    call prec%av(psb_l_pr_)%trim()
    call prec%av(psb_u_pr_)%trim()
    call prec%dv%bld(dd)
    ! call move_alloc(dd,prec%d)
  else
    info=psb_err_from_subroutine_
    ch_err='psb_iluk_fact'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  case(psb_f_ilu_t_)
  ! ILU(N,E) Incomplete LU factorization with thresholding and level of fill

    if (allocated(prec%av)) then
      if (size(prec%av) < psb_bp_ilu_avsz) then
        do i=1,size(prec%av)
          call prec%av(i)%free()
        enddo
        deallocate(prec%av,stat=info)
      endif
    end if
    if (.not.allocated(prec%av)) then
      allocate(prec%av(psb_max_avsz),stat=info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      end if
    endif

    nrow_a = desc_a%get_local_rows()
    nztota = a%get_nzeros()

    n_col  = desc_a%get_local_cols()
    nhalo  = n_col-nrow_a
    n_row  = nrow_a

    allocate(lf,uf,stat=info)
    if (info == psb_success_) call lf%allocate(n_row,n_row,nztota)
    if (info == psb_success_) call uf%allocate(n_row,n_row,nztota)

    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    allocate(dd(n_row),stat=info)
    if (info == psb_success_) then
      allocate(prec%dv, stat=info)
      if (info == 0) then
        if (present(vmold)) then
          allocate(prec%dv%v,mold=vmold,stat=info)
        else
          allocate(psb_d_base_vect_type :: prec%dv%v,stat=info)
        end if
      end if
    end if

    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999
    endif
    ! This is where we have no renumbering, thus no need
    call psb_ilut_fact(fill_in,fact_eps,a,lf,uf,dd,info,iscale=iscale)

    if(info == psb_success_) then
      call prec%av(psb_l_pr_)%mv_from(lf%a)
      call prec%av(psb_u_pr_)%mv_from(uf%a)
      call prec%av(psb_l_pr_)%set_asb()
      call prec%av(psb_u_pr_)%set_asb()
      call prec%av(psb_l_pr_)%trim()
      call prec%av(psb_u_pr_)%trim()
      call prec%dv%bld(dd)
      ! call move_alloc(dd,prec%d)
    else
      info=psb_err_from_subroutine_
      ch_err='psb_ilut_fact'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case(psb_f_ainv_)
    ! Approximate Inverse Factorizations based on variants of the incomplete
    ! biconjugation algorithms
    if (allocated(prec%av)) then
      if (size(prec%av) < psb_bp_ilu_avsz) then
        do i=1,size(prec%av)
          call prec%av(i)%free()
        enddo
        deallocate(prec%av,stat=info)
      endif
    end if
    if (.not.allocated(prec%av)) then
      allocate(prec%av(psb_max_avsz),stat=info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      end if
    endif

    nrow_a = desc_a%get_local_rows()
    nztota = a%get_nzeros()

    n_col  = desc_a%get_local_cols()
    nhalo  = n_col-nrow_a
    n_row  = nrow_a

    allocate(lf,uf,stat=info)
    if (info == psb_success_) call lf%allocate(n_row,n_row,nztota)
    if (info == psb_success_) call uf%allocate(n_row,n_row,nztota)

    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    allocate(dd(n_row),stat=info)
    if (info == psb_success_) then
      allocate(prec%dv, stat=info)
      if (info == 0) then
        if (present(vmold)) then
          allocate(prec%dv%v,mold=vmold,stat=info)
        else
          allocate(psb_d_base_vect_type :: prec%dv%v,stat=info)
        end if
      end if
    end if

    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999
    endif
    ! Computing the factorization
    call psb_ainv_fact(a,iinvalg,inv_fill,inv_thresh,lf,dd,uf,desc_a,info,iscale=iscale)

    if(info == psb_success_) then
      call prec%av(psb_l_pr_)%mv_from(lf%a)
      call prec%av(psb_u_pr_)%mv_from(uf%a)
      call prec%av(psb_l_pr_)%set_asb()
      call prec%av(psb_u_pr_)%set_asb()
      call prec%av(psb_l_pr_)%trim()
      call prec%av(psb_u_pr_)%trim()
      call prec%dv%bld(dd)
      ! call move_alloc(dd,prec%d)
    else
      info=psb_err_from_subroutine_
      ch_err='psb_ilut_fact'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case(psb_f_invk_)
    ! Approximate Inverse Factorizations based on the sparse inversion of
    ! triangular factors of an ILU(n) factorization
    if (allocated(prec%av)) then
      if (size(prec%av) < psb_bp_ilu_avsz) then
        do i=1,size(prec%av)
          call prec%av(i)%free()
        enddo
        deallocate(prec%av,stat=info)
      endif
    end if
    if (.not.allocated(prec%av)) then
      allocate(prec%av(psb_max_avsz),stat=info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      end if
    endif

    nrow_a = desc_a%get_local_rows()
    nztota = a%get_nzeros()

    n_col  = desc_a%get_local_cols()
    nhalo  = n_col-nrow_a
    n_row  = nrow_a

    allocate(lf,uf,stat=info)
    if (info == psb_success_) call lf%allocate(n_row,n_row,nztota)
    if (info == psb_success_) call uf%allocate(n_row,n_row,nztota)

    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    allocate(dd(n_row),stat=info)
    if (info == psb_success_) then
      allocate(prec%dv, stat=info)
      if (info == 0) then
        if (present(vmold)) then
          allocate(prec%dv%v,mold=vmold,stat=info)
        else
          allocate(psb_d_base_vect_type :: prec%dv%v,stat=info)
        end if
      end if
    end if

    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999
    endif
    ! Computing the factorization
    call psb_invk_fact(a,fill_in, inv_fill,lf,dd,uf,desc_a,info)

    if(info == psb_success_) then
      call prec%av(psb_l_pr_)%mv_from(lf%a)
      call prec%av(psb_u_pr_)%mv_from(uf%a)
      call prec%av(psb_l_pr_)%set_asb()
      call prec%av(psb_u_pr_)%set_asb()
      call prec%av(psb_l_pr_)%trim()
      call prec%av(psb_u_pr_)%trim()
      call prec%dv%bld(dd)
      ! call move_alloc(dd,prec%d)
    else
      info=psb_err_from_subroutine_
      ch_err='psb_ilut_fact'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case(psb_f_invt_)
    ! Approximate Inverse Factorizations based on the sparse inversion of
    ! triangular factors of an ILU(eps) factorization
    if (allocated(prec%av)) then
      if (size(prec%av) < psb_bp_ilu_avsz) then
        do i=1,size(prec%av)
          call prec%av(i)%free()
        enddo
        deallocate(prec%av,stat=info)
      endif
    end if
    if (.not.allocated(prec%av)) then
      allocate(prec%av(psb_max_avsz),stat=info)
      if (info /= psb_success_) then
        call psb_errpush(psb_err_alloc_dealloc_,name)
        goto 9999
      end if
    endif

    nrow_a = desc_a%get_local_rows()
    nztota = a%get_nzeros()

    n_col  = desc_a%get_local_cols()
    nhalo  = n_col-nrow_a
    n_row  = nrow_a

    allocate(lf,uf,stat=info)
    if (info == psb_success_) call lf%allocate(n_row,n_row,nztota)
    if (info == psb_success_) call uf%allocate(n_row,n_row,nztota)

    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    allocate(dd(n_row),stat=info)
    if (info == psb_success_) then
      allocate(prec%dv, stat=info)
      if (info == 0) then
        if (present(vmold)) then
          allocate(prec%dv%v,mold=vmold,stat=info)
        else
          allocate(psb_d_base_vect_type :: prec%dv%v,stat=info)
        end if
      end if
    end if

    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
      goto 9999
    endif
    ! Computing the factorization
    call psb_invt_fact(a,fill_in,inv_fill,fact_eps,inv_thresh,lf,dd,uf,&
      & desc_a,info)

    if(info == psb_success_) then
      call prec%av(psb_l_pr_)%mv_from(lf%a)
      call prec%av(psb_u_pr_)%mv_from(uf%a)
      call prec%av(psb_l_pr_)%set_asb()
      call prec%av(psb_u_pr_)%set_asb()
      call prec%av(psb_l_pr_)%trim()
      call prec%av(psb_u_pr_)%trim()
      call prec%dv%bld(dd)
      ! call move_alloc(dd,prec%d)
    else
      info=psb_err_from_subroutine_
      ch_err='psb_ilut_fact'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  case(psb_f_none_)
    info=psb_err_from_subroutine_
    ch_err='Inconsistent prec  psb_f_none_'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999

  case default
    info=psb_err_from_subroutine_
    ch_err='Unknown psb_f_type_'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end select

  if (present(amold)) then
    call prec%av(psb_l_pr_)%cscnv(info,mold=amold)
    call prec%av(psb_u_pr_)%cscnv(info,mold=amold)
  end if

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_d_bjac_precbld

subroutine psb_d_bjac_precseti(prec,what,val,info)

  use psb_base_mod
  use psb_d_bjacprec, psb_protect_name => psb_d_bjac_precseti
  Implicit None

  class(psb_d_bjac_prec_type),intent(inout) :: prec
  integer(psb_ipk_), intent(in)                      :: what
  integer(psb_ipk_), intent(in)                      :: val
  integer(psb_ipk_), intent(out)                     :: info
  integer(psb_ipk_) :: err_act, nrow
  character(len=20)  :: name='d_bjac_precset'

  call psb_erractionsave(err_act)

  info = psb_success_
  if (.not.allocated(prec%iprcparm)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner")
    goto 9999
  end if

  select case(what)
  case (psb_f_type_)
    prec%iprcparm(psb_f_type_)     = val

  case (psb_ilu_fill_in_)
    prec%iprcparm(psb_ilu_fill_in_) = val

  case (psb_ilu_ialg_)
    prec%iprcparm(psb_ilu_ialg_) = val

  case (psb_ilu_scale_)
    prec%iprcparm(psb_ilu_scale_) = val

  case (psb_ainv_alg_)
    prec%iprcparm(psb_ainv_alg_) = val

  case(psb_inv_fillin_)
    prec%iprcparm(psb_inv_fillin_) = val

  case default
    write(psb_err_unit,'(i0," is invalid, ignoring user specification")') what

  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_d_bjac_precseti

subroutine psb_d_bjac_precsetr(prec,what,val,info)

  use psb_base_mod
  use psb_d_bjacprec, psb_protect_name => psb_d_bjac_precsetr
  Implicit None

  class(psb_d_bjac_prec_type),intent(inout) :: prec
  integer(psb_ipk_), intent(in)                      :: what
  real(psb_dpk_), intent(in)                       :: val
  integer(psb_ipk_), intent(out)                     :: info
  integer(psb_ipk_) :: err_act, nrow
  character(len=20)  :: name='d_bjac_precset'

  call psb_erractionsave(err_act)

  info = psb_success_
  if (.not.allocated(prec%rprcparm)) then
    info = 1124
    call psb_errpush(info,name,a_err="preconditioner")
    goto 9999
  end if

  select case(what)

  case (psb_fact_eps_)
    prec%rprcparm(psb_fact_eps_) = val

  case (psb_inv_thresh_)
    prec%rprcparm(psb_inv_thresh_) = val

  case default
    write(psb_err_unit,'(i0," is invalid, ignoring user specification")') what

  end select

  call psb_erractionrestore(err_act)
  return

9999 call psb_error_handler(err_act)
  return
end subroutine psb_d_bjac_precsetr
