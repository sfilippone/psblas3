!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
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
! File: psb_dspsm.f90
!
! Subroutine: psb_dspsm
!  Performs one of the distributed matrix-vector operations
!
!     sub( Y ) := alpha * Pr * A-1 * Pc *sub( X ) + beta * sub (Y ),   or
!
!     sub( Y ) := alpha * D * Pr * A-1 * Pc * sub( X ) + beta * sub (Y ), or
!
!     sub( Y ) := alpha * Pr * A-1 * Pc * D * sub( X ) + beta * sub (Y ),  or
!
!     sub( Y ) := alpha * Pr * A-T * Pc * sub( X ) + beta * sub (Y ),   or
!
!     sub( Y ) := alpha * D * Pr * A-T * Pc * sub( X ) + beta * sub (Y ), or
!
!     sub( Y ) := alpha * Pr * A-T * Pc * D * sub( X ) + beta * sub (Y ),  or
!
!  where :
!
!        sub( X ) denotes X(1:M,JX:JX+K-1),
!      
!        sub( Y ) denotes Y(1:M,JY:JY+K-1).
!
!  sub( X ) is a distributed
!  vector and T is a M-by-M distributed triangular matrix.
!
! Arguments:   
!    alpha   -  real                   The scalar alpha.
!    a       -  type(psb_zspmat_type). The sparse matrix containing A.
!    x(:,:)  -  real                   The input vector containing the entries of ( X ).
!    beta    -  real                   The scalar beta.
!    y(:,:)  -  real                   The input vector containing the entries of ( Y ).
!    desc_a  -  type(psb_desc_type).   The communication descriptor.
!    info    -  integer.               Return code
!    trans   -  character(optional).   Whether A or A'. If not present 'N' is assumed.
!    scale   -  character(optional).   Specify some type of operation with
!                                      the diagonal matrix D.
!    choice  -  integer(optional).     The kind of update to perform on overlap elements.
!    d(:)    -  real  , optional       Matrix for diagonal scaling.
!    k       -  integer(optional).     The number of right-hand sides.
!    jx      -  integer(optional).     The column offset for ( X ). Default: 1
!    jy      -  integer(optional).     The column offset for ( Y ). Default: 1 
!    work(:) -  real  , optional       Working area.
!
! 
subroutine  psb_dspsm(alpha,a,x,beta,y,desc_a,info,&
     & trans, scale, choice, diag, k, jx, jy, work)   
  use psb_base_mod, psb_protect_name => psb_dspsm
  use psi_mod
  implicit none

  real(psb_dpk_), intent(in)              :: alpha, beta
  real(psb_dpk_), intent(in), target      :: x(:,:)
  real(psb_dpk_), intent(inout), target   :: y(:,:)
  type(psb_dspmat_type), intent(in)        :: a
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)                      :: info
  real(psb_dpk_), intent(in), optional, target      :: diag(:)
  real(psb_dpk_), optional, target, intent(inout)       :: work(:)
  character, intent(in), optional           :: trans, scale
  integer(psb_ipk_), intent(in), optional             :: choice
  integer(psb_ipk_), intent(in), optional             :: k, jx, jy

  ! locals
  integer(psb_ipk_) :: ictxt, np, me,&
       & err_act, iix, jjx, ia, ja, iia, jja, lldx,lldy, choice_,&
       & ix, iy, ik, ijx, ijy, i, lld,&
       & m, nrow, ncol, liwork, llwork, iiy, jjy, idx, ndm

  character                :: lscale
  integer(psb_ipk_), parameter       :: nb=4
  real(psb_dpk_),pointer :: iwork(:), xp(:,:), yp(:,:), id(:)
  character                :: itrans
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_dspsm'
  if (psb_errstatus_fatal()) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! just this case right now
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

  if (present(k)) then
    ik = min(k,size(x,2)-ijx+1)
    ik = min(ik,size(y,2)-ijy+1)
  else
    ik = min(size(x,2)-ijx+1,size(y,2)-ijy+1)
  endif

  if (present(choice)) then     
    choice_ = choice
  else
    choice_ = psb_avg_
  endif

  if (present(scale)) then     
    lscale = psb_toupper(scale)
  else
    lscale = 'U'
  endif

  if (present(trans)) then     
    itrans = psb_toupper(trans)
    if((itrans == 'N').or.(itrans == 'T').or.  (itrans == 'C')) then
      ! OK 
    else
      info = psb_err_iarg_invalid_value_
      call psb_errpush(info,name)
      goto 9999
    end if
  else
    itrans = 'N'
  endif

  m    = desc_a%get_global_rows()
  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  lldx = size(x,1)
  lldy = size(y,1)

  if((lldx < ncol).or.(lldy < ncol)) then
    info=psb_err_lld_case_not_implemented_
    call psb_errpush(info,name)
    goto 9999
  end if

  ! check for presence/size of a work area
  iwork => null()
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

  iwork(1)=0.d0

  if(present(diag)) then
    lld = size(diag)
    id => diag
  else
    lld=1
    allocate(id(1))
    id=1.d0
  end if

  ! checking for matrix correctness
  call psb_chkmat(m,m,ia,ja,desc_a,info,iia,jja)
  ! checking for vectors correctness
  if (info == psb_success_) &
       & call psb_chkvect(m,ik,size(x,1),ix,ijx,desc_a,info,iix,jjx)
  if (info == psb_success_) &
       & call psb_chkvect(m,ik,size(y,1),iy,ijy,desc_a,info,iiy,jjy)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect/mat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if(ja /= ix) then
    ! this case is not yet implemented
    info = psb_err_ja_nix_ia_niy_unsupported_
  end if

  if((iix /= 1).or.(iiy /= 1)) then
    ! this case is not yet implemented
    info = psb_err_ix_n1_iy_n1_unsupported_
  end if

  if(info /= psb_success_) then
    call psb_errpush(info,name)
    goto 9999
  end if

  ! Perform local triangular system solve
  xp => x(iix:lldx,jjx:jjx+ik-1)
  yp => y(iiy:lldy,jjy:jjy+ik-1)
  call psb_cssm(alpha,a,xp,beta,yp,info,scale=scale,d=diag,trans=trans)

  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err='cssm'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  ! update overlap elements
  if (choice_ > 0) then

    call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),ik,&
         & done,yp,desc_a,iwork,info,data=psb_comm_ovr_)

    if (info == psb_success_) call psi_ovrl_upd(yp,desc_a,choice_,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Inner updates')
      goto 9999
    end if
  end if

  if(aliw) deallocate(iwork)
  if(.not.present(diag)) deallocate(id)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_dspsm
     
!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
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
!
! Subroutine: psb_dspsv
!  Performs one of the distributed matrix-vector operations
!
!     Y := alpha * Pr * A-1 * Pc * X + beta * Y,   or
!
!     Y := alpha * D * Pr * A-1 * Pc * X + beta * Y, or
!
!     Y := alpha * Pr * A-1 * Pc * D * X + beta * Y,  or
!
!     Y := alpha * Pr * A-T * Pc * X + beta * Y,   or
!
!     Y := alpha * D * Pr * A-T * Pc * X + beta * Y, or
!
!     Y := alpha * Pr * A-T * Pc * D * X + beta * Y,  or
!
!  X is a distributed
!  vector and T is a M-by-M distributed triangular matrix.
!
!
! Arguments:   
!    alpha   -  real                   The scalar alpha.
!    a       -  type(psb_zspmat_type). The sparse matrix containing A.
!    x(:)    -  real                   The input vector containing the entries of ( X ).
!    beta    -  real                   The scalar beta.
!    y(:)    -  real                   The input vector containing the entries of ( Y ).
!    desc_a  -  type(psb_desc_type).   The communication descriptor.
!    info    -  integer.               Return code
!    trans   -  character(optional).   Whether A or A'. If not present 'N' is assumed.
!    scale   -  character(optional).   Specify some type of operation with
!                                      the diagonal matrix D.
!    choice  -  integer(optional).     The kind of update to perform on overlap elements.
!    d(:)    -  real  , optional       Matrix for diagonal scaling.
!    work(:) -  real  , optional       Working area.
! 
subroutine  psb_dspsv(alpha,a,x,beta,y,desc_a,info,&
     & trans, scale, choice, diag, work)   
  use psb_base_mod, psb_protect_name => psb_dspsv
  use psi_mod
  implicit none 

  real(psb_dpk_), intent(in)              :: alpha, beta
  real(psb_dpk_), intent(in), target      :: x(:)
  real(psb_dpk_), intent(inout), target   :: y(:)
  type(psb_dspmat_type), intent(in)         :: a
  type(psb_desc_type), intent(in)           :: desc_a
  integer(psb_ipk_), intent(out)                      :: info
  real(psb_dpk_), intent(in), optional, target    :: diag(:)
  real(psb_dpk_), optional, target, intent(inout)        :: work(:)
  character, intent(in), optional           :: trans, scale
  integer(psb_ipk_), intent(in), optional             :: choice

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, &
       & err_act, iix, jjx, ia, ja, iia, jja, lldx,lldy, choice_,&
       & ix, iy, ik, jx, jy, i, lld,&
       & m, nrow, ncol, liwork, llwork, iiy, jjy, idx, ndm

  character                :: lscale
  integer(psb_ipk_), parameter       :: nb=4
  real(psb_dpk_),pointer :: iwork(:), xp(:), yp(:), id(:)
  character                :: itrans
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_dspsv'
  if (psb_errstatus_fatal()) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = psb_err_context_error_
    call psb_errpush(info,name)
    goto 9999
  endif

  ! just this case right now
  ia = 1
  ja = 1
  ix = 1
  iy = 1
  ik = 1
  jx= 1
  jy= 1

  if (present(choice)) then     
    choice_ = choice
  else
    choice_ = psb_avg_
  endif

  if (present(scale)) then     
    lscale = psb_toupper(scale)
  else
    lscale = 'U'
  endif

  if (present(trans)) then     
    itrans = psb_toupper(trans)
    if((itrans == 'N').or.(itrans == 'T').or.(itrans == 'C')) then
      ! Ok
    else
      info = psb_err_iarg_invalid_value_
      call psb_errpush(info,name)
      goto 9999
    end if
  else
    itrans = 'N'
  endif

  m    = desc_a%get_global_rows()
  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  lldx = size(x)
  lldy = size(y)

  if((lldx < ncol).or.(lldy < ncol)) then
    info=psb_err_lld_case_not_implemented_
    call psb_errpush(info,name)
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
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  else
    iwork => work
  endif

  iwork(1)=0.d0

  if (present(diag)) then
    lld = size(diag)
    id => diag
  else
    lld=1
    allocate(id(1))
    id=1.d0
  end if

  ! checking for matrix correctness
  call psb_chkmat(m,m,ia,ja,desc_a,info,iia,jja)
  ! checking for vectors correctness
  if (info == psb_success_) &
       & call psb_chkvect(m,ik,size(x),ix,jx,desc_a,info,iix,jjx)
  if (info == psb_success_)&
       & call psb_chkvect(m,ik,size(y),iy,jy,desc_a,info,iiy,jjy)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect/mat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if(ja /= ix) then
    ! this case is not yet implemented
    info = psb_err_ja_nix_ia_niy_unsupported_
  end if

  if((iix /= 1).or.(iiy /= 1)) then
    ! this case is not yet implemented
    info = psb_err_ix_n1_iy_n1_unsupported_
  end if

  if(info /= psb_success_) then
    call psb_errpush(info,name)
    goto 9999
  end if

  ! Perform local triangular system solve
  xp => x(iix:lldx)
  yp => y(iiy:lldy)
  call psb_cssm(alpha,a,xp,beta,yp,info,scale=scale,d=diag,trans=trans)

  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err='dcssm'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  ! update overlap elements
  if (choice_ > 0) then
    call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
         & done,yp,desc_a,iwork,info,data=psb_comm_ovr_)


    if (info == psb_success_) call psi_ovrl_upd(yp,desc_a,choice_,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Inner updates')
      goto 9999
    end if
  end if

  if (aliw) deallocate(iwork)
  if(.not.present(diag)) deallocate(id)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_dspsv
     
subroutine  psb_dspsv_vect(alpha,a,x,beta,y,desc_a,info,&
     & trans, scale, choice, diag, work)   
  use psb_base_mod, psb_protect_name => psb_dspsv_vect
  use psi_mod
  implicit none 

  real(psb_dpk_), intent(in)              :: alpha, beta
  type(psb_d_vect_type), intent(inout)    :: x
  type(psb_d_vect_type), intent(inout)    :: y
  type(psb_dspmat_type), intent(inout)    :: a
  type(psb_desc_type), intent(in)         :: desc_a
  integer(psb_ipk_), intent(out)                    :: info
  type(psb_d_vect_type), intent(inout), optional  :: diag
  real(psb_dpk_), optional, target, intent(inout)        :: work(:)
  character, intent(in), optional           :: trans, scale
  integer(psb_ipk_), intent(in), optional             :: choice

  ! locals
  integer(psb_ipk_) :: ictxt, np, me, &
       & err_act, iix, jjx, ia, ja, iia, jja, lldx,lldy, choice_,&
       & ix, iy, ik, jx, jy, i, lld,&
       & m, nrow, ncol, liwork, llwork, iiy, jjy, idx, ndm

  character                :: lscale
  integer(psb_ipk_), parameter       :: nb=4
  real(psb_dpk_),pointer :: iwork(:), xp(:), yp(:)
  character                :: itrans
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_dspsv'
  if (psb_errstatus_fatal()) return 
  info=psb_success_
  call psb_erractionsave(err_act)

  ictxt=desc_a%get_context()

  call psb_info(ictxt, me, np)
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


  ! just this case right now
  ia = 1
  ja = 1
  ix = 1
  iy = 1
  ik = 1
  jx = 1
  jy = 1

  if (present(choice)) then     
    choice_ = choice
  else
    choice_ = psb_avg_
  endif

  if (present(scale)) then     
    lscale = psb_toupper(scale)
  else
    lscale = 'U'
  endif

  if (present(trans)) then     
    itrans = psb_toupper(trans)
    if((itrans == 'N').or.(itrans == 'T').or.(itrans == 'C')) then
      ! Ok
    else
      info = psb_err_iarg_invalid_value_
      call psb_errpush(info,name)
      goto 9999
    end if
  else
    itrans = 'N'
  endif

  m    = desc_a%get_global_rows()
  nrow = desc_a%get_local_rows()
  ncol = desc_a%get_local_cols()
  lldx = x%get_nrows()
  lldy = y%get_nrows()

  if((lldx < ncol).or.(lldy < ncol)) then
    info=psb_err_lld_case_not_implemented_
    call psb_errpush(info,name)
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
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  else
    iwork => work
  endif

  iwork(1)=0.d0

  ! checking for matrix correctness
  call psb_chkmat(m,m,ia,ja,desc_a,info,iia,jja)
  ! checking for vectors correctness
  if (info == psb_success_) &
       & call psb_chkvect(m,ik,x%get_nrows(),ix,jx,desc_a,info,iix,jjx)
  if (info == psb_success_)&
       & call psb_chkvect(m,ik,y%get_nrows(),iy,jy,desc_a,info,iiy,jjy)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    ch_err='psb_chkvect/mat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if(ja /= ix) then
    ! this case is not yet implemented
    info = psb_err_ja_nix_ia_niy_unsupported_
  end if

  if((iix /= 1).or.(iiy /= 1)) then
    ! this case is not yet implemented
    info = psb_err_ix_n1_iy_n1_unsupported_
  end if

  if(info /= psb_success_) then
    call psb_errpush(info,name)
    goto 9999
  end if

  ! Perform local triangular system solve
!!$  if (present(diag)) then 
!!$    call psb_cssm(alpha,a,x,beta,y,info,scale=scale,d=diag,trans=trans)    
!!$  else
!!$    call psb_cssm(alpha,a,x,beta,y,info,scale=scale,trans=trans)    
!!$  end if
  if (present(diag)) then 
    call a%cssm(alpha,x,beta,y,info,scale=scale,d=diag,trans=trans)    
  else
    call a%cssm(alpha,x,beta,y,info,scale=scale,trans=trans)    
  end if
  if(info /= psb_success_) then
    info = psb_err_from_subroutine_
    ch_err='dcssm'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  ! update overlap elements
  if (choice_ > 0) then
    call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
         & done,y%v,desc_a,iwork,info,data=psb_comm_ovr_)


    if (info == psb_success_) call psi_ovrl_upd(y%v,desc_a,choice_,info)
    if (info /= psb_success_) then
      call psb_errpush(psb_err_from_subroutine_,name,a_err='Inner updates')
      goto 9999
    end if
  end if

  if (aliw) deallocate(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_dspsv_vect
     
