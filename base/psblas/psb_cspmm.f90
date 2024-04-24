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
! File: psb_cspmm.f90
!
!
! Subroutine: psb_cspm_vect
!     Performs one of the distributed matrix-vector operations
!
!     Y := alpha * Pr * A * Pc * X  + beta * Y,  or
!
!     Y := alpha * Pr * A' * Pr * X  + beta * Y,
!
!  alpha and beta are scalars, X and Y are distributed
!  vectors and A is a M-by-N distributed matrix.
!
! Arguments:   
!    alpha   -  complex                The scalar alpha.
!    a       -  type(psb_cspmat_type). The sparse matrix containing A.
!    x       -  type(psb_c_vect_type) The input vector containing the entries of ( X ).
!    beta    -  complex                The scalar beta.
!    y       -  type(psb_c_vect_type) The input vector containing the entries of ( Y ).
!    desc_a  -  type(psb_desc_type).   The communication descriptor.
!    info    -  integer.               Return code
!    trans   -  character(optional).   Whether A or A'. Default:  'N' 
!    work(:) -  complex,(optional).    Working area.
!    doswap  -  logical(optional).     Whether to performe halo updates.
! 
subroutine  psb_cspmv_vect(alpha,a,x,beta,y,desc_a,info,&
     & trans, work, doswap)   
  use psb_base_mod, psb_protect_name => psb_cspmv_vect
  use psi_mod
  implicit none

  complex(psb_spk_), intent(in)            :: alpha, beta
  type(psb_c_vect_type), intent(inout)     :: x
  type(psb_c_vect_type), intent(inout)     :: y
  type(psb_cspmat_type), intent(in)        :: a
  type(psb_desc_type), intent(in)          :: desc_a
  integer(psb_ipk_), intent(out)                     :: info
  complex(psb_spk_), optional, target, intent(inout) :: work(:)
  character, intent(in), optional          :: trans
  logical, intent(in), optional            :: doswap

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, iia, jja,  nrow, ncol, lldx, lldy, &
       & liwork, iiy, jjy, ib, ip, idx
  integer(psb_lpk_) :: ix, ijx, iy, ijy, m, n, ia, ja
  integer(psb_ipk_), parameter       :: nb=4
  complex(psb_spk_), pointer :: iwork(:), xp(:), yp(:)
  complex(psb_spk_), allocatable :: xvsave(:)
  character                :: trans_
  character(len=20)        :: name, ch_err
  logical                  :: aliw, doswap_
  integer(psb_ipk_) :: debug_level, debug_unit
  logical, parameter  :: do_timings=.true.
  integer(psb_ipk_), save  :: mv_phase1=-1, mv_phase2=-1, mv_phase3=-1, mv_phase4=-1
  integer(psb_ipk_), save  :: mv_phase11=-1, mv_phase12=-1

  name='psb_cspmv'
  info=psb_success_
  call psb_erractionsave(err_act)
  if  (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt=desc_a%get_context()
  call psb_info(ctxt, me, np)
  if (np == -1) then
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

  if (present(doswap)) then
    doswap_ = doswap
  else
    doswap_ = .true.
  endif

  if (present(trans)) then     
    trans_ = psb_toupper(trans)
  else
    trans_ = 'N'
  endif
  if ( (trans_ == 'N').or.(trans_ == 'T')&
       & .or.(trans_ == 'C')) then
  else
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info,name)
    goto 9999
  end if
  if ((do_timings).and.(mv_phase1==-1))       &
       & mv_phase1 = psb_get_timer_idx("SPMM: and send ")
  if ((do_timings).and.(mv_phase2==-1))       &
       & mv_phase2 = psb_get_timer_idx("SPMM: and cmp ad")
  if ((do_timings).and.(mv_phase3==-1))       &
       & mv_phase3 = psb_get_timer_idx("SPMM: and rcv")
  if ((do_timings).and.(mv_phase4==-1))       &
       & mv_phase4 = psb_get_timer_idx("SPMM: and cmp and")
  if ((do_timings).and.(mv_phase11==-1))       &
       & mv_phase11 = psb_get_timer_idx("SPMM: noand exch ")
  if ((do_timings).and.(mv_phase12==-1))       &
       & mv_phase12 = psb_get_timer_idx("SPMM: noand cmp")


  m    = desc_a%get_global_rows()
  n    = desc_a%get_global_cols()
  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  lldx = x%get_nrows()
  lldy = y%get_nrows()
  if ((info == 0).and.(lldx<ncol)) call x%reall(ncol,info)
  if ((info == 0).and.(lldy<ncol)) call y%reall(ncol,info)

  if (psb_errstatus_fatal()) then 
    info=psb_err_from_subroutine_
    ch_err='reall'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  iwork => null()
  ! check for presence/size of a work area
  liwork= 2*ncol

  if (present(work)) then
    if (size(work) >= liwork) then
      aliw =.false.
    else
      aliw=.true.
    endif
  else
    aliw=.true.
  end if

  if (aliw) then
    allocate(iwork(liwork),stat=info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='Allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  else
    iwork => work
  endif

  if (debug_level >= psb_debug_comp_) &
       & write(debug_unit,*) me,' ',trim(name),' Allocated work ', info

  if (trans_ == 'N') then
    !  Matrix is not transposed
    
    if (allocated(a%ad)) then
      block
        logical, parameter :: do_timings=.true.
        real(psb_dpk_) :: t1, t2, t3, t4, t5
        !if (me==0) write(0,*) 'going for overlap ',a%ad%get_fmt(),' ',a%and%get_fmt()
        if (do_timings) call psb_barrier(ctxt)
        if (do_timings) call psb_tic(mv_phase1)
        if (doswap_) call psi_swapdata(psb_swap_send_,&
             & czero,x%v,desc_a,iwork,info,data=psb_comm_halo_)
        if (do_timings) call psb_toc(mv_phase1)
        if (do_timings) call psb_tic(mv_phase2)          
        call a%ad%spmm(alpha,x%v,beta,y%v,info)
        if (do_timings) call psb_tic(mv_phase3)
        if (doswap_) call psi_swapdata(psb_swap_recv_,&
             & czero,x%v,desc_a,iwork,info,data=psb_comm_halo_)
        if (do_timings) call psb_toc(mv_phase3)
        if (do_timings) call psb_tic(mv_phase4)          
        call a%and%spmm(alpha,x%v,cone,y%v,info,ivshft=a%ad%get_ncols())
        if (do_timings) call psb_toc(mv_phase4)
      end block
      
    else
      block
        logical, parameter :: do_timings=.true.
        real(psb_dpk_) :: t1, t2, t3, t4, t5
        if (do_timings) call psb_barrier(ctxt)
        
        if (do_timings) call psb_tic(mv_phase11)          
        if (doswap_) then
          call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
               & czero,x%v,desc_a,iwork,info,data=psb_comm_halo_)
        end if
        if (do_timings) call psb_toc(mv_phase11)
        if (do_timings) call psb_tic(mv_phase12)          
        call psb_csmm(alpha,a,x,beta,y,info)
        if (do_timings) call psb_toc(mv_phase12)
      end block
    end if
    
    if(info /= psb_success_) then
      info = psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if

  else
    !  Matrix is transposed

    !
    ! Non-empty overlap, need a buffer to hold
    ! the entries updated with average operator.
    ! Why the average? because in this way they will contribute
    ! with a proper scale factor (1/np) to the overall product.
    ! 
    call psi_ovrl_save(x%v,xvsave,desc_a,info)
    if (info == psb_success_) call psi_ovrl_upd(x%v,desc_a,psb_avg_,info)

    if (beta /= czero) call y%set(czero,nrow+1,ncol)
    !  local Matrix-vector product
    if (info == psb_success_) call psb_csmm(alpha,a,x,beta,y,info,trans=trans_)

    if (debug_level >= psb_debug_comp_) &
         & write(debug_unit,*) me,' ',trim(name),' csmm ', info

    if (info == psb_success_) call psi_ovrl_restore(x%v,xvsave,desc_a,info)
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      ch_err='psb_csmm'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if (doswap_) then
      call psi_swaptran(ior(psb_swap_send_,psb_swap_recv_),&
           & cone,y%v,desc_a,iwork,info)
      if (info == psb_success_) call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
           & cone,y%v,desc_a,iwork,info,data=psb_comm_ovr_)

      if (debug_level >= psb_debug_comp_) &
           & write(debug_unit,*) me,' ',trim(name),' swaptran ', info
      if(info /= psb_success_) then
        info = psb_err_from_subroutine_
        ch_err='PSI_SwapTran'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    end if

  end if

  if (aliw) deallocate(iwork,stat=info)
  if (debug_level >= psb_debug_comp_) &
       & write(debug_unit,*) me,' ',trim(name),' deallocat ',aliw, info
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err='Deallocate iwork'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  nullify(iwork)

  call psb_erractionrestore(err_act)
  if (debug_level >= psb_debug_comp_) then 
    call psb_barrier(ctxt)
    write(debug_unit,*) me,' ',trim(name),' Returning '
  endif
  return  

9999 call psb_error_handler(ctxt,err_act)

  return
end subroutine psb_cspmv_vect
!
! Subroutine: psb_cspmm
!     Performs one of the distributed matrix-vector operations
!
!     sub( Y ) := alpha * Pr * A * Pc * sub( X )  + beta * sub( Y ),  or
!
!     sub( Y ) := alpha * Pr * A' * Pr * sub( X )  + beta * sub( Y ),
!
!     where:
!
!        sub( X ) denotes:  X(1:N,JX:JX+K-1),
!
!        sub( Y ) denotes:  Y(1:M,JY:JY+K-1),
!
!  alpha and beta are scalars, and sub( X ) and sub( Y ) are distributed
!  vectors and A is a M-by-N distributed matrix.
!
! Arguments:   
!    alpha   -  complex                The scalar alpha.
!    a       -  type(psb_cspmat_type). The sparse matrix containing A.
!    x(:,:)  -  complex                The input vector containing the entries of ( X ).
!    beta    -  complex                The scalar beta.
!    y(:,:)  -  complex                The input vector containing the entries of ( Y ).
!    desc_a  -  type(psb_desc_type).   The communication descriptor.
!    info    -  integer.               Return code
!    trans   -  character(optional).   Whether A or A'. Default: 'N' 
!    k       -  integer(optional).     The number of right-hand sides.
!    jx      -  integer(optional).     The column offset for ( X ). Default:  1
!    jy      -  integer(optional).     The column offset for ( Y ). Default:  1
!    work(:) -  complex,(optional).    Working area.
!    doswap  -  logical(optional).     Whether to performe halo updates.
! 
subroutine  psb_cspmm(alpha,a,x,beta,y,desc_a,info,&
     & trans, k, jx, jy, work, doswap)   
  use psb_base_mod, psb_protect_name => psb_cspmm
  use psi_mod
  implicit none

  complex(psb_spk_), intent(in)             :: alpha, beta
  complex(psb_spk_), intent(inout), target  :: x(:,:)
  complex(psb_spk_), intent(inout), target  :: y(:,:)
  type(psb_cspmat_type), intent(in)        :: a
  type(psb_desc_type), intent(in)          :: desc_a
  integer(psb_ipk_), intent(out)                     :: info
  complex(psb_spk_), optional, target, intent(inout)  :: work(:)
  character, intent(in), optional          :: trans
  integer(psb_ipk_), intent(in), optional            :: k, jx, jy
  logical, intent(in), optional            :: doswap

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, iia, jja,  nrow, ncol, lldx, lldy, &
       & liwork, iiy, jjy, i, ib, ib1, ip, idx, ik
  integer(psb_lpk_) :: ix, ijx, iy, ijy, m, n, ia, ja, lik
  integer(psb_ipk_), parameter               :: nb=4
  complex(psb_spk_), pointer     :: xp(:,:), yp(:,:), iwork(:)
  complex(psb_spk_), allocatable :: xvsave(:,:)
  character                        :: trans_
  character(len=20)                :: name, ch_err
  logical                          :: aliw, doswap_
  integer(psb_ipk_) :: debug_level, debug_unit

  name='psb_cspmm'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt=desc_a%get_context()

  call psb_info(ctxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ia = 1
  ja = 1

  ix = 1
  if (present(jx)) then
    ijx = jx
  else
    ijx = 1
  endif

  iy = 1
  if (present(jy)) then
    ijy = jy
  else
    ijy = 1
  endif

  if (present(doswap)) then
    doswap_ = doswap
  else
    doswap_ = .true.
  endif

  if (present(k)) then     
    lik = min(k,size(x,2)-ijx+1)
    lik = min(lik,size(y,2)-ijy+1)
  else
    lik = min(size(x,2)-ijx+1,size(y,2)-ijy+1)
  endif

  if (present(trans)) then     
    trans_ = psb_toupper(trans)
  else
    trans_ = 'N'
  endif
  if ( (trans_ == 'N').or.(trans_ == 'T')&
       & .or.(trans_ == 'C')) then
  else
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info,name)
    goto 9999
  end if

  m    = desc_a%get_global_rows()
  n    = desc_a%get_global_cols()
  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  lldx = size(x,1)
  lldy = size(y,1)

  ! check for presence/size of a work area
  liwork= 2*ncol

  if (present(work)) then
    if (size(work) >= liwork) then
      aliw =.false.
    else
      aliw=.true.
    endif
  else
    aliw=.true.
  end if

  if (aliw) then
    allocate(iwork(liwork),stat=info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  else
    iwork => work
  endif

  iwork(1)=czero

  ! checking for matrix correctness
  call psb_chkmat(m,n,ia,ja,desc_a,info,iia,jja)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkmat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if


  if (trans_ == 'N') then
    !  Matrix is not transposed
    if((ja /= ix).or.(ia /= iy)) then
      ! this case is not yet implemented
      info = psb_err_ja_nix_ia_niy_unsupported_
      call psb_errpush(info,name)
      goto 9999
    end if

    ! checking for vectors correctness
    call psb_chkvect(n,lik,lldx,ix,ijx,desc_a,info,iix,jjx,check_halo=.true.)
    if (info == psb_success_) &
         & call psb_chkvect(m,lik,lldy,iy,ijy,desc_a,info,iiy,jjy)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_chkvect'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if((iix /= 1).or.(iiy /= 1)) then
      ! this case is not yet implemented
      info = psb_err_ix_n1_iy_n1_unsupported_
      call psb_errpush(info,name)
      goto 9999
    end if


    if (doswap_.and.(np>1)) then 
      ib1=min(nb,lik)
      xp => x(iix:lldx,jjx:jjx+ib1-1)
      if (doswap_)&
           & call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
           & ib1,czero,xp,desc_a,iwork,info)


      blk: do i=1, lik, nb
        ib=ib1
        ib1 = max(0,min(nb,(lik)-(i-1+ib)))
        xp => x(iix:lldx,jjx+i-1+ib:jjx+i-1+ib+ib1-1)
        if ((ib1 > 0).and.(doswap_)) &
             & call psi_swapdata(psb_swap_send_,ib1,&
             & czero,xp,desc_a,iwork,info)

        if(info /= psb_success_) exit blk

        !  local Matrix-vector product
        call psb_csmm(alpha,a,x(:,jjx+i-1:jjx+i-1+ib-1),&
             & beta,y(:,jjy+i-1:jjy+i-1+ib-1),info,trans=trans_)

        if(info /= psb_success_) exit blk

        if((ib1 > 0).and.(doswap_))&
             & call psi_swapdata(psb_swap_recv_,ib1,&
             & czero,xp,desc_a,iwork,info)

        if(info /= psb_success_) exit blk
      end do blk
    else
      if (doswap_)&
           & call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
           & ib1,czero,x(:,1:lik),desc_a,iwork,info)
      if (info == psb_success_) call psb_csmm(alpha,a,x(:,1:lik),beta,y(:,1:lik),info)
    end if
    if(info /= psb_success_) then
      info = psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if

  else

    !  Matrix is transposed
    if((ja /= iy).or.(ia /= ix)) then
      ! this case is not yet implemented
      info = psb_err_ja_nix_ia_niy_unsupported_
      call psb_errpush(info,name)
      goto 9999
    end if


    ! checking for vectors correctness
    call psb_chkvect(m,lik,lldx,ix,ijx,desc_a,info,iix,jjx,check_halo=.true.)
    if (info == psb_success_) &
         & call psb_chkvect(n,lik,lldy,iy,ijy,desc_a,info,iiy,jjy,check_halo=.true.)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_chkvect'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if((iix /= 1).or.(iiy /= 1)) then
      ! this case is not yet implemented
      info = psb_err_ix_n1_iy_n1_unsupported_
      call psb_errpush(info,name)
      goto 9999
    end if

    !
    ! Non-empty overlap, need a buffer to hold
    ! the entries updated with average operator.
    ! Why the average? because in this way they will contribute
    ! with a proper scale factor (1/np) to the overall product.
    ! 
    call psi_ovrl_save(x(:,1:lik),xvsave,desc_a,info)
    if (info == psb_success_) call psi_ovrl_upd(x,desc_a,psb_avg_,info)
    y(nrow+1:ncol,1:lik)    = czero

    if (info == psb_success_) &
         & call psb_csmm(alpha,a,x(:,1:lik),beta,y(:,1:lik),info,trans=trans_)
    if (debug_level >= psb_debug_comp_) &
         & write(debug_unit,*) me,' ',trim(name),' csmm ', info
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      ch_err='psb_csmm'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    if (info == psb_success_) call psi_ovrl_restore(x,xvsave,desc_a,info)

    if (doswap_)then
      ik = lik ! This should not be an issue, we are expecting the values
      ! to be small, within IPK
      call psi_swaptran(ior(psb_swap_send_,psb_swap_recv_),&
           & ik,cone,y(:,1:ik),desc_a,iwork,info)
      if (info == psb_success_) call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
           & ik,cone,y(:,1:ik),desc_a,iwork,info,data=psb_comm_ovr_)

      if (debug_level >= psb_debug_comp_) &
           & write(debug_unit,*) me,' ',trim(name),' swaptran ', info
      if(info /= psb_success_) then
        info = psb_err_from_subroutine_
        ch_err='PSI_dSwapTran'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    end if

  end if

  if (aliw) deallocate(iwork,stat=info)
  if (debug_level >= psb_debug_comp_) &
       & write(debug_unit,*) me,' ',trim(name),' deallocat ',aliw, info
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err='Deallocate iwork'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 call psb_error_handler(ctxt,err_act)

  return
end subroutine psb_cspmm
!!$ 
!!$              Parallel Sparse BLAS  version 3.5
!!$    (C) Copyright 2006-2018
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari      
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
!
! Subroutine: psb_cspmv
!     Performs one of the distributed matrix-vector operations
!
!     Y := alpha * Pr * A * Pc * X  + beta * Y,  or
!
!     Y := alpha * Pr * A' * Pr * X  + beta * Y,
!
!  alpha and beta are scalars, and X and Y are distributed
!  vectors and A is a M-by-N distributed matrix.
!
! Arguments:   
!    alpha   -  complex                The scalar alpha.
!    a       -  type(psb_cspmat_type). The sparse matrix containing A.
!    x(:)    -  complex                The input vector containing the entries of ( X ).
!    beta    -  complex                The scalar beta.
!    y(:)    -  complex                The input vector containing the entries of ( Y ).
!    desc_a  -  type(psb_desc_type).   The communication descriptor.
!    info    -  integer.               Return code
!    trans   -  character(optional).   Whether A or A'. Default:  'N' 
!    work(:) -  complex,(optional).    Working area.
!    doswap  -  logical(optional).     Whether to performe halo updates.
! 
subroutine  psb_cspmv(alpha,a,x,beta,y,desc_a,info,&
     & trans, work, doswap)   
  use psb_base_mod, psb_protect_name => psb_cspmv
  use psi_mod
  implicit none

  complex(psb_spk_), intent(in)             :: alpha, beta
  complex(psb_spk_), intent(inout), target  :: x(:)
  complex(psb_spk_), intent(inout), target  :: y(:)
  type(psb_cspmat_type), intent(in)        :: a
  type(psb_desc_type), intent(in)          :: desc_a
  integer(psb_ipk_), intent(out)                     :: info
  complex(psb_spk_), optional, target, intent(inout) :: work(:)
  character, intent(in), optional          :: trans
  logical, intent(in), optional            :: doswap

  ! locals
  type(psb_ctxt_type) :: ctxt
  integer(psb_ipk_) :: np, me,&
       & err_act, iix, jjx, iia, jja, nrow, ncol, lldx, lldy, &
       & liwork, iiy, jjy, ib, ip, idx, ik
  integer(psb_lpk_) :: ix, ijx, iy, ijy, m, n, ia, ja, lik, jx, jy
  integer(psb_ipk_), parameter           :: nb=4
  complex(psb_spk_), pointer :: iwork(:), xp(:), yp(:)
  complex(psb_spk_), allocatable :: xvsave(:)  
  character                    :: trans_
  character(len=20)            :: name, ch_err
  logical                      :: aliw, doswap_
  integer(psb_ipk_) :: debug_level, debug_unit

  name='psb_cspmv'
  info=psb_success_
  call psb_erractionsave(err_act)
  if (psb_errstatus_fatal()) then
    info = psb_err_internal_error_ ;    goto 9999
  end if
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  ctxt=desc_a%get_context()
  call psb_info(ctxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ia = 1
  ja = 1
  ix = 1
  jx = 1
  iy = 1
  jy = 1
  ik = 1
  lik = 1
  ib = 1

  if (present(doswap)) then
    doswap_ = doswap
  else
    doswap_ = .true.
  endif

  if (present(trans)) then     
    trans_ = psb_toupper(trans)
  else
    trans_ = 'N'
  endif
  if ( (trans_ == 'N').or.(trans_ == 'T')&
       & .or.(trans_ == 'C')) then
  else
    info = psb_err_iarg_invalid_value_
    call psb_errpush(info,name)
    goto 9999
  end if

  m    = desc_a%get_global_rows()
  n    = desc_a%get_global_cols()
  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  lldx = size(x)
  lldy = size(y)

  iwork => null()
  ! check for presence/size of a work area
  liwork= 2*ncol

  if (present(work)) then
    if (size(work) >= liwork) then
      aliw =.false.
    else
      aliw=.true.
    endif
  else
    aliw=.true.
  end if

  if (aliw) then
    allocate(iwork(liwork),stat=info)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='Allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  else
    iwork => work
  endif

  if (debug_level >= psb_debug_comp_) &
       & write(debug_unit,*) me,' ',trim(name),' Allocated work ', info
  ! checking for matrix correctness
  call psb_chkmat(m,n,ia,ja,desc_a,info,iia,jja)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkmat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (debug_level >= psb_debug_comp_) &
       & write(debug_unit,*) me,' ',trim(name),' Checkmat ', info
  if (trans_ == 'N') then
    !  Matrix is not transposed
    if((ja /= ix).or.(ia /= iy)) then
      ! this case is not yet implemented
      info = psb_err_ja_nix_ia_niy_unsupported_
      call psb_errpush(info,name)
      goto 9999
    end if

    ! checking for vectors correctness
    call psb_chkvect(n,lik,lldx,ix,jx,desc_a,info,iix,jjx,check_halo=.true.)
    if (info == psb_success_) &
         & call psb_chkvect(m,lik,lldy,iy,jy,desc_a,info,iiy,jjy)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_chkvect'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if((iix /= 1).or.(iiy /= 1)) then
      ! this case is not yet implemented
      info = psb_err_ix_n1_iy_n1_unsupported_
      call psb_errpush(info,name)
      goto 9999
    end if

    if (doswap_) then
      call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
           & czero,x,desc_a,iwork,info,data=psb_comm_halo_)
    end if

    call psb_csmm(alpha,a,x,beta,y,info)

    if(info /= psb_success_) then
      info = psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if

  else
    !  Matrix is transposed
    if((ja /= iy).or.(ia /= ix)) then
      ! this case is not yet implemented
      info = psb_err_ja_nix_ia_niy_unsupported_
      call psb_errpush(info,name)
      goto 9999
    end if

    ! checking for vectors correctness
    call psb_chkvect(m,lik,lldx,ix,jx,desc_a,info,iix,jjx,check_halo=.true.)
    if (info == psb_success_)&
         & call psb_chkvect(n,lik,lldy,iy,jy,desc_a,info,iiy,jjy,check_halo=.true.)
    if(info /= psb_success_) then
      info=psb_err_from_subroutine_
      ch_err='psb_chkvect'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if((iix /= 1).or.(iiy /= 1)) then
      ! this case is not yet implemented
      info = psb_err_ix_n1_iy_n1_unsupported_
      call psb_errpush(info,name)
      goto 9999
    end if

    xp => x(1:lldx)
    yp => y(1:lldy)

    !
    ! Non-empty overlap, need a buffer to hold
    ! the entries updated with average operator.
    ! Why the average? because in this way they will contribute
    ! with a proper scale factor (1/np) to the overall product.
    ! 
    call psi_ovrl_save(x,xvsave,desc_a,info)
    if (info == psb_success_) call psi_ovrl_upd(x,desc_a,psb_avg_,info)
    yp(nrow+1:ncol) = czero

    !  local Matrix-vector product
    if (info == psb_success_) call psb_csmm(alpha,a,x,beta,y,info,trans=trans_)

    if (debug_level >= psb_debug_comp_) &
         & write(debug_unit,*) me,' ',trim(name),' csmm ', info

    if (info == psb_success_) call psi_ovrl_restore(x,xvsave,desc_a,info)
    if (info /= psb_success_) then
      info = psb_err_from_subroutine_
      ch_err='psb_csmm'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if (doswap_) then
      call psi_swaptran(ior(psb_swap_send_,psb_swap_recv_),&
           & cone,yp,desc_a,iwork,info)
      if (info == psb_success_) call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
           & cone,yp,desc_a,iwork,info,data=psb_comm_ovr_)

      if (debug_level >= psb_debug_comp_) &
           & write(debug_unit,*) me,' ',trim(name),' swaptran ', info
      if(info /= psb_success_) then
        info = psb_err_from_subroutine_
        ch_err='PSI_dSwapTran'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if
    end if

  end if

  if (aliw) deallocate(iwork,stat=info)
  if (debug_level >= psb_debug_comp_) &
       & write(debug_unit,*) me,' ',trim(name),' deallocat ',aliw, info
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err='Deallocate iwork'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  nullify(iwork)

  call psb_erractionrestore(err_act)
  if (debug_level >= psb_debug_comp_) then 
    call psb_barrier(ctxt)
    write(debug_unit,*) me,' ',trim(name),' Returning '
  endif
  return  

9999 call psb_error_handler(ctxt,err_act)

  return
end subroutine psb_cspmv
