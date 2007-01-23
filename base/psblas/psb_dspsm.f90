!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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
! Parameters:   
!    alpha  -  real.                        The scalar alpha.
!    a      -  type(<psb_dspmat_type>).     The sparse matrix containing A.
!    x      -  real,dimension(:,:).         The input vector containing the entries of sub( X ).
!    beta   -  real.                        The scalar beta.
!    y      -  real,dimension(:,:).         The input vector containing the entries of sub( Y ).
!    desc_a -  type(<psb_desc_type>).       The communication descriptor.
!    info   -  integer.                     Eventually returns an error code.
!    trans  -  character(optional).         Whether A or A'. If not present 'N' is assumed.
!    unitd  -  character(optional).         Specify some type of operation with the diagonal matrix D.
!    choice -  integer(optional).           The kind of update to perform on overlap elements.
!    d      -  real,dimension(:)(optional). Matrix for diagonal scaling.
!    k      -  integer(optional).           The number of right-hand sides.
!    jx     -  integer(optional).           The column offset for sub( X ). If not present 1 is assumed.
!    jy     -  integer(optional).           The column offset for sub( Y ). If not present 1 is assumed.
!    work   -  real,dimension(:)(optional). Working area.
! 
subroutine  psb_dspsm(alpha,a,x,beta,y,desc_a,info,&
     & trans, unitd, choice, d, k, jx, jy, work)   

  use psb_spmat_type
  use psb_serial_mod
  use psb_descriptor_type
  use psb_comm_mod
  use psi_mod
  use psb_check_mod
  use psb_error_mod
  use psb_string_mod
  use psb_penv_mod
  implicit none

  real(kind(1.D0)), intent(in)              :: alpha, beta
  real(kind(1.d0)), intent(in), target      :: x(:,:)
  real(kind(1.d0)), intent(inout), target   :: y(:,:)
  type (psb_dspmat_type), intent(in)        :: a
  type(psb_desc_type), intent(in)           :: desc_a
  integer, intent(out)                      :: info
  real(kind(1.d0)), intent(in), optional, target      :: d(:)
  real(kind(1.d0)), optional, target       :: work(:)
  character, intent(in), optional           :: trans, unitd
  integer, intent(in), optional             :: choice
  integer, intent(in), optional             :: k, jx, jy

  ! locals
  integer                  :: int_err(5), ictxt, np, me,&
       & err_act, n, iix, jjx, ia, ja, iia, jja, lldx,lldy, lchoice,&
       & ix, iy, ik, ijx, ijy, i, lld,&
       & m, nrow, ncol, liwork, llwork, iiy, jjy

  character                :: lunitd
  integer, parameter       :: nb=4
  real(kind(1.d0)),pointer :: iwork(:), xp(:,:), yp(:,:), id(:)
  character                :: itrans
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_dspsm'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
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
    lchoice = choice
  else
    lchoice = psb_avg_
  endif

  if (present(unitd)) then     
    lunitd = toupper(unitd)
  else
    lunitd = 'U'
  endif

  if (present(trans)) then     
    itrans = toupper(trans)
    if((itrans.eq.'N').or.(itrans.eq.'T').or.  (itrans.eq.'C')) then
      ! OK 
    else
      info = 70
      call psb_errpush(info,name)
      goto 9999
    end if
  else
    itrans = 'N'
  endif

  m    = psb_cd_get_global_rows(desc_a)
  nrow = psb_cd_get_local_rows(desc_a)
  ncol = psb_cd_get_local_cols(desc_a)

  lldx = size(x,1)
  lldy = size(y,1)

  if((lldx.lt.ncol).or.(lldy.lt.ncol)) then
    info=3010
    call psb_errpush(info,name)
    goto 9999
  end if

  ! check for presence/size of a work area
  iwork => null()
  liwork= 2*ncol
  if (a%pr(1) /= 0) llwork = liwork + m * ik
  if (a%pl(1) /= 0) llwork = llwork + m * ik
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
    if(info.ne.0) then
      info=4010
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  else
    iwork => work
  endif

  iwork(1)=0.d0

  if(present(d)) then
    lld = size(d)
    id => d
  else
    lld=1
    allocate(id(1))
    id=1.d0
  end if

  ! checking for matrix correctness
  call psb_chkmat(m,m,ia,ja,desc_a,info,iia,jja)
  ! checking for vectors correctness
  if (info == 0) &
       & call psb_chkvect(m,ik,size(x,1),ix,ijx,desc_a,info,iix,jjx)
  if (info == 0) &
       & call psb_chkvect(m,ik,size(y,1),iy,ijy,desc_a,info,iiy,jjy)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chkvect/mat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if(ja.ne.ix) then
    ! this case is not yet implemented
    info = 3030
  end if

  if((iix.ne.1).or.(iiy.ne.1)) then
    ! this case is not yet implemented
    info = 3040
  end if

  if(info.ne.0) then
    call psb_errpush(info,name)
    goto 9999
  end if

  ! Perform local triangular system solve
  xp => x(iix:lldx,jjx:jjx+ik-1)
  yp => y(iiy:lldy,jjy:jjy+ik-1)
  call psb_cssm(alpha,a,xp,beta,yp,info,unitd=lunitd,d=id,trans=itrans)

  if(info.ne.0) then
    info = 4010
    ch_err='dcssm'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  ! update overlap elements
  if(lchoice.gt.0) then

    call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),ik,&
         & done,yp,desc_a,iwork,info)

    i=0
    ! switch on update type
    select case (lchoice)
    case(psb_square_root_)
      do while(desc_a%ovrlap_elem(i).ne.-ione)
        y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:) =&
             & y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:)/&
             & sqrt(real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_)))
        i = i+2
      end do
    case(psb_avg_)
      do while(desc_a%ovrlap_elem(i).ne.-ione)
        y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:) =&
             & y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_),:)/&
             & real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_))
        i = i+2
      end do
    case(psb_sum_)
      ! do nothing
    case default 
      ! wrong value for choice argument
      info = 70
      int_err=(/10,lchoice,0,0,0/)
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end select
  end if

  if(aliw) deallocate(iwork)
  if(.not.present(d)) deallocate(id)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_dspsm
     
!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
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
! Subroutine: psb_dspsmv
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
! Parameters:   
!    alpha  -  real.                        The scalar alpha.
!    a      -  type(<psb_dspmat_type>).     The sparse matrix containing A.
!    x      -  real,dimension(:).           The input vector containing the entries of X.
!    beta   -  real.                        The scalar beta.
!    y      -  real,dimension(:).           The input vector containing the entries of Y.
!    desc_a -  type(<psb_desc_type>).       The communication descriptor.
!    info   -  integer.                     Eventually returns an error code.
!    trans  -  character(optional).         Whether A or A'. If not present 'N' is assumed.
!    unitd  -  character(optional).         Specify some type of operation with the diagonal matrix D.
!    choice -  integer(optional).           The kind of update to perform on overlap elements.
!    d      -  real,dimension(:)(optional). Matrix for diagonal scaling.
!    work   -  real,dimension(:)(optional). Working area.
! 
subroutine  psb_dspsv(alpha,a,x,beta,y,desc_a,info,&
     & trans, unitd, choice, d, work)   
  use psb_spmat_type
  use psb_serial_mod
  use psb_descriptor_type
  use psb_comm_mod
  use psi_mod
  use psb_check_mod
  use psb_error_mod
  use psb_string_mod
  use psb_penv_mod
  implicit none 

  real(kind(1.D0)), intent(in)              :: alpha, beta
  real(kind(1.d0)), intent(in), target      :: x(:)
  real(kind(1.d0)), intent(inout), target   :: y(:)
  type(psb_dspmat_type), intent(in)         :: a
  type(psb_desc_type), intent(in)           :: desc_a
  integer, intent(out)                      :: info
  real(kind(1.d0)), intent(in), optional, target    :: d(:)
  real(kind(1.d0)), optional, target        :: work(:)
  character, intent(in), optional           :: trans, unitd
  integer, intent(in), optional             :: choice

  ! locals
  integer                  :: int_err(5), ictxt, np, me,&
       & err_act, n, iix, jjx, ia, ja, iia, jja, lldx,lldy, lchoice,&
       & ix, iy, ik, jx, jy, i, lld,&
       & m, nrow, ncol, liwork, llwork, iiy, jjy

  character                :: lunitd
  integer, parameter       :: nb=4
  real(kind(1.d0)),pointer :: iwork(:), xp(:), yp(:), id(:)
  character                :: itrans
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_dspsv'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
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
    lchoice = choice
  else
    lchoice = psb_avg_
  endif

  if (present(unitd)) then     
    lunitd = toupper(unitd)
  else
    lunitd = 'U'
  endif

  if (present(trans)) then     
    itrans = toupper(trans)
    if((itrans.eq.'N').or.(itrans.eq.'T')) then
      ! Ok
    else if (itrans.eq.'C') then
      info = 3020
      call psb_errpush(info,name)
      goto 9999
    else
      info = 70
      call psb_errpush(info,name)
      goto 9999
    end if
  else
    itrans = 'N'
  endif

  m    = psb_cd_get_global_rows(desc_a)
  nrow = psb_cd_get_local_rows(desc_a)
  ncol = psb_cd_get_local_cols(desc_a)

  lldx = size(x)
  lldy = size(y)

  if((lldx.lt.ncol).or.(lldy.lt.ncol)) then
    info=3010
    call psb_errpush(info,name)
    goto 9999
  end if

  iwork => null()
  ! check for presence/size of a work area
  liwork= 2*ncol
  if (a%pr(1) /= 0) llwork = liwork + m * ik
  if (a%pl(1) /= 0) llwork = llwork + m * ik

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
    if(info.ne.0) then
      info=4010
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  else
    iwork => work
  endif

  iwork(1)=0.d0

  if (present(d)) then
    lld = size(d)
    id => d
  else
    lld=1
    allocate(id(1))
    id=1.d0
  end if

  ! checking for matrix correctness
  call psb_chkmat(m,m,ia,ja,desc_a,info,iia,jja)
  ! checking for vectors correctness
  if (info == 0) &
       & call psb_chkvect(m,ik,size(x),ix,jx,desc_a,info,iix,jjx)
  if (info == 0)&
       & call psb_chkvect(m,ik,size(y),iy,jy,desc_a,info,iiy,jjy)
  if(info.ne.0) then
    info=4010
    ch_err='psb_chkvect/mat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if(ja.ne.ix) then
    ! this case is not yet implemented
    info = 3030
  end if

  if((iix.ne.1).or.(iiy.ne.1)) then
    ! this case is not yet implemented
    info = 3040
  end if

  if(info.ne.0) then
    call psb_errpush(info,name)
    goto 9999
  end if

  ! Perform local triangular system solve
  xp => x(iix:lldx)
  yp => y(iiy:lldy)
  call psb_cssm(alpha,a,xp,beta,yp,info,unitd=lunitd,d=id,trans=itrans)

  if(info.ne.0) then
    info = 4010
    ch_err='dcssm'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  ! update overlap elements
  if(lchoice.gt.0) then
    call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
         & done,yp,desc_a,iwork,info)

    i=0
    ! switch on update type
    select case (lchoice)
    case(psb_square_root_)
      do while(desc_a%ovrlap_elem(i).ne.-ione)
        y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_)) =&
             & y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_))/&
             & sqrt(real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_)))
        i = i+2
      end do
    case(psb_avg_)
      do while(desc_a%ovrlap_elem(i).ne.-ione)
        y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_)) =&
             & y(desc_a%ovrlap_elem(i+psb_ovrlp_elem_))/&
             & real(desc_a%ovrlap_elem(i+psb_n_dom_ovr_))
        i = i+2
      end do
    case(psb_sum_)
      ! do nothing
    case default 
      ! wrong value for choice argument
      info = 70
      int_err=(/10,lchoice,0,0,0/)
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end select
  end if

  if (aliw) deallocate(iwork)
  if(.not.present(d)) deallocate(id)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_dspsv
     
