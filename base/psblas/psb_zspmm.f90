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
! File: psb_zspmm.f90
!
! Subroutine: psb_zspmm
!     Performs one of the distributed matrix-vector operations
!
!     sub( Y ) := alpha * Pr * A * Pc * sub( X )  + beta * sub( Y ),  or
!
!     sub( Y ) := alpha * Pr * A' * Pr * sub( X )  + beta * sub( Y ),
!
!     where:
!
!        sub( X ) denotes *if* TRANS = 'N',
!
!                       X(1:N,JX:JX+K-1),
!
!                     *else*
!
!                       X(1:M,JX:JX+K-1).
!
!                     *end if*
!
!        sub( Y ) denotes *if* trans = 'N',
!
!                       Y(1:M,JY:JY+K-1),
!
!                     *else*
!
!                       Y(1:N,JY:JY+K-1)
!
!                     *end* *if*
!
!  alpha and beta are scalars, and sub( X ) and sub( Y ) are distributed
!  vectors and A is a M-by-N distributed matrix.
!
! Parameters:   
!    alpha  -  real.                        The scalar alpha.
!    a      -  type(<psb_zspmat_type>).     The sparse matrix containing A.
!    x      -  real,dimension(:,:).         The input vector containing the entries of sub( X ).
!    beta   -  real.                        The scalar beta.
!    y      -  real,dimension(:,:).         The input vector containing the entries of sub( Y ).
!    desc_a -  type(<psb_desc_type>).       The communication descriptor.
!    info   -  integer.                     Eventually returns an error code.
!    trans  -  character(optional).         Whether A or A'. If not present 'N' is assumed.
!    k      -  integer(optional).           The number of right-hand sides.
!    jx     -  integer(optional).           The column offset for sub( X ). If not present 1 is assumed.
!    jy     -  integer(optional).           The column offset for sub( Y ). If not present 1 is assumed.
!    work   -  real,dimension(:)(optional). Working area.
!    doswap -  integer(optional).           Whether to performe halo updates.
! 
subroutine  psb_zspmm(alpha,a,x,beta,y,desc_a,info,&
     & trans, k, jx, jy, work, doswap)   

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

  complex(kind(1.D0)), intent(in)             :: alpha, beta
  complex(kind(1.d0)), intent(inout), target  :: x(:,:)
  complex(kind(1.d0)), intent(inout), target  :: y(:,:)
  type(psb_zspmat_type), intent(in)        :: a
  type(psb_desc_type), intent(in)          :: desc_a
  integer, intent(out)                     :: info
  complex(kind(1.d0)), optional, target      :: work(:)
  character, intent(in), optional          :: trans
  integer, intent(in), optional            :: k, jx, jy,doswap

  ! locals
  integer                  :: ictxt, np, me,&
       & err_act, n, iix, jjx, ia, ja, iia, jja, ix, iy, ik, ijx, ijy,&
       & idoswap, m, nrow, ncol, lldx, lldy, liwork, iiy, jjy,&
       & i, ib, ib1
  integer, parameter          :: nb=4
  complex(kind(1.d0)),pointer :: xp(:,:), yp(:,:), iwork(:)
  character                :: itrans
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_zspmm'
  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
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
    idoswap = doswap
  else
    idoswap = 1
  endif

  if (present(k)) then     
    ik = min(k,size(x,2)-ijx+1)
    ik = min(ik,size(y,2)-ijy+1)
  else
    ik = min(size(x,2)-ijx+1,size(y,2)-ijy+1)
  endif

  if (present(trans)) then     
    if ( (toupper(trans) == 'N').or.(toupper(trans) == 'T').or. (toupper(trans) == 'C')) then
      itrans = toupper(trans)
    else
      info = 70
      call psb_errpush(info,name)
      goto 9999
    end if
  else
    itrans = 'N'
  endif

  m    = psb_cd_get_global_rows(desc_a)
  n    = psb_cd_get_global_cols(desc_a)
  nrow = psb_cd_get_local_rows(desc_a)
  ncol = psb_cd_get_local_cols(desc_a)
  lldx = size(x,1)
  lldy = size(y,1)

  ! check for presence/size of a work area
  liwork= 2*ncol
  if (a%pr(1) /= 0) liwork = liwork + n * ik
  if (a%pl(1) /= 0) liwork = liwork + m * ik
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
    if(info /= 0) then
      info=4010
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  else
    iwork => work
  endif

  iwork(1)=zzero

  ! checking for matrix correctness
  call psb_chkmat(m,n,ia,ja,desc_a,info,iia,jja)
  if(info /= 0) then
    info=4010
    ch_err='psb_chkmat'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if


  if (itrans == 'N') then
    !  Matrix is not transposed
    if((ja /= ix).or.(ia /= iy)) then
      ! this case is not yet implemented
      info = 3030
      call psb_errpush(info,name)
      goto 9999
    end if

    ! checking for vectors correctness
    call psb_chkvect(n,ik,size(x,1),ix,ijx,desc_a,info,iix,jjx)
    if (info == 0)&
         & call psb_chkvect(m,ik,size(y,1),iy,ijy,desc_a,info,iiy,jjy)
    if(info /= 0) then
      info=4010
      ch_err='psb_chkvect'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if((iix /= 1).or.(iiy /= 1)) then
      ! this case is not yet implemented
      info = 3040
      call psb_errpush(info,name)
      goto 9999
    end if


    ib1=min(nb,ik)
    xp => x(iix:lldx,jjx:jjx+ib1-1)
    if(idoswap > 0)&
         & call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
         & ib1,zzero,xp,desc_a,iwork,info)


    blk: do i=1, ik, nb
      ib=ib1
      ib1 = max(0,min(nb,(ik)-(i-1+ib)))
      xp => x(iix:lldx,jjx+i+ib-1:jjx+i+ib+ib1-2)
      if((ib1 > 0).and.(idoswap > 0))&
           & call psi_swapdata(psb_swap_send_,ib1,&
           & zzero,xp,desc_a,iwork,info)

      if(info /= 0) exit blk

      !  local Matrix-vector product
      call psb_csmm(alpha,a,x(iix:lldx,jjx:jjx+ib-1),&
           & beta,y(iiy:lldy,jjy:jjy+ib-1),info,trans=itrans)

      if(info /= 0) exit blk

      if((ib1 > 0).and.(idoswap > 0))&
           & call psi_swapdata(psb_swap_send_,ib1,&
           & zzero,xp,desc_a,iwork,info)

      if(info /= 0) exit blk
    end do blk

    if(info /= 0) then
      info = 4011
      call psb_errpush(info,name)
      goto 9999
    end if

  else
    !  Matrix is transposed
    if((ja /= iy).or.(ia /= ix)) then
      ! this case is not yet implemented
      info = 3030
      call psb_errpush(info,name)
      goto 9999
    end if

    if(desc_a%ovrlap_elem(1) /= -1) then
      info = 3070
      call psb_errpush(info,name)
      goto 9999
    end if

    ! checking for vectors correctness
    call psb_chkvect(m,ik,size(x,1),ix,ijx,desc_a,info,iix,jjx)
    if (info == 0) &
         & call psb_chkvect(n,ik,size(y,1),iy,ijy,desc_a,info,iiy,jjy)
    if(info /= 0) then
      info=4010
      ch_err='psb_chkvect'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    if((iix /= 1).or.(iiy /= 1)) then
      ! this case is not yet implemented
      info = 3040
      call psb_errpush(info,name)
      goto 9999
    end if

    y(iiy+nrow+1-1:iiy+ncol,1:ik)=zzero

    !  local Matrix-vector product

    call psb_csmm(alpha,a,x(iix:lldx,jjx:jjx+ik-1),&
         & beta,y(iiy:lldy,jjy:jjy+ik-1),info,trans=itrans)

    if(info /= 0) then
      info = 4010
      ch_err='csmm'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    yp => y(iiy:lldy,jjy:jjy+ik-1)
    if (idoswap/=0)&
         & call psi_swaptran(ior(psb_swap_send_,psb_swap_recv_),&
         & ik,zone,yp,desc_a,iwork,info)

    if(info /= 0) then
      info = 4010
      ch_err='PSI_dSwapTran'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

  end if

  if(aliw) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return
end subroutine psb_zspmm




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
! Subroutine: psb_zspmmv
!     Performs one of the distributed matrix-vector operations
!
!     Y := alpha * Pr * A * Pc * X  + beta * Y,  or
!
!     Y := alpha * Pr * A' * Pr * X  + beta * Y,
!
!  alpha and beta are scalars, and X and Y are distributed
!  vectors and A is a M-by-N distributed matrix.
!
! Parameters:   
!    alpha  -  real.                        The scalar alpha.
!    a      -  type(<psb_zspmat_type>).     The sparse matrix containing A.
!    x      -  real,dimension(:).           The input vector containing the entries of X.
!    beta   -  real.                        The scalar beta.
!    y      -  real,dimension(:.         The input vector containing the entries of Y.
!    desc_a -  type(<psb_desc_type>).       The communication descriptor.
!    info   -  integer.                     Eventually returns an error code.
!    trans  -  character(optional).         Whether A or A'. If not present 'N' is assumed.
!    work   -  real,dimension(:)(optional). Working area.
!    doswap -  integer(optional).           Whether to performe halo updates.
! 
subroutine  psb_zspmv(alpha,a,x,beta,y,desc_a,info,&
     & trans, work, doswap)   

  use psb_spmat_type
  use psb_serial_mod
  use psb_descriptor_type
  use psb_comm_mod
  use psb_const_mod
  use psi_mod
  use psb_check_mod
  use psb_error_mod
  use psb_string_mod
  use psb_penv_mod
  implicit none

  complex(kind(1.D0)), intent(in)             :: alpha, beta
  complex(kind(1.d0)), intent(inout), target  :: x(:)
  complex(kind(1.d0)), intent(inout), target  :: y(:)
  type(psb_zspmat_type), intent(in)        :: a
  type(psb_desc_type), intent(in)          :: desc_a
  integer, intent(out)                     :: info
  complex(kind(1.d0)), optional, target      :: work(:)
  character, intent(in), optional          :: trans
  integer, intent(in), optional            :: doswap

  ! locals
  integer                  :: ictxt, np, me,&
       & err_act, n, iix, jjx, ia, ja, iia, jja, ix, iy, ik, ijx, ijy,&
       & idoswap, m, nrow, ncol, lldx, lldy, liwork, jx, jy, iiy, jjy,&
       & i, ib, ib1
  integer, parameter       :: nb=4
  complex(kind(1.d0)),pointer :: iwork(:), xp(:), yp(:)
  character                :: itrans
  character(len=20)        :: name, ch_err
  logical                  :: aliw

  name='psb_zspmv'
  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
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
  ib = 1

  if (present(doswap)) then
     idoswap = doswap
  else
     idoswap = 1
  endif

  if (present(trans)) then     
    if ( (toupper(trans) == 'N').or.(toupper(trans) == 'T') .or.(toupper(trans) == 'C')) then
      itrans = toupper(trans)
    else
      info = 70
      call psb_errpush(info,name)
      goto 9999
    end if
  else
     itrans = 'N'
  endif

  m    = psb_cd_get_global_rows(desc_a)
  n    = psb_cd_get_global_cols(desc_a)
  nrow = psb_cd_get_local_rows(desc_a)
  ncol = psb_cd_get_local_cols(desc_a)
  lldx = size(x)
  lldy = size(y)

  iwork => null()
  ! check for presence/size of a work area
  liwork= 2*ncol
  if (a%pr(1) /= 0) liwork = liwork + n * ik
  if (a%pl(1) /= 0) liwork = liwork + m * ik
  !  write(0,*)'---->>>',work(1)
 if (present(work)) then
    if (size(work) >= liwork) then
        aliw =.false.
    else
        aliw=.true.
    endif
  else
        aliw=.true.
  end if

        aliw=.true.
  if (aliw) then
    allocate(iwork(liwork),stat=info)
    if(info /= 0) then
      info=4010
      ch_err='Allocate'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
  else
     iwork => work
  endif


  ! checking for matrix correctness
  call psb_chkmat(m,n,ia,ja,desc_a,info,iia,jja)
  if(info /= 0) then
     info=4010
     ch_err='psb_chkmat'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if


  if (itrans == 'N') then
     !  Matrix is not transposed
     if((ja /= ix).or.(ia /= iy)) then
        ! this case is not yet implemented
        info = 3030
        call psb_errpush(info,name)
        goto 9999
     end if

     ! checking for vectors correctness
     call psb_chkvect(n,ik,size(x),ix,jx,desc_a,info,iix,jjx)
     if (info == 0) &
          & call psb_chkvect(m,ik,size(y),iy,jy,desc_a,info,iiy,jjy)
     if(info /= 0) then
        info=4010
        ch_err='psb_chkvect'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

     if((iix /= 1).or.(iiy /= 1)) then
        ! this case is not yet implemented
        info = 3040
        call psb_errpush(info,name)
        goto 9999
     end if

     if (idoswap == 0) then
        x(nrow+1:ncol)=zzero
     else
        call psi_swapdata(ior(psb_swap_send_,psb_swap_recv_),&
             & zzero,x,desc_a,iwork,info,data=psb_comm_halo_)
     end if

     !  local Matrix-vector product
      call psb_csmm(alpha,a,x(iix:lldx),beta,y(iiy:lldy),info)

     if(info /= 0) then
        info = 4011
        call psb_errpush(info,name)
        goto 9999
     end if

  else
     !  Matrix is transposed
     if((ja /= iy).or.(ia /= ix)) then
        ! this case is not yet implemented
        info = 3030
        call psb_errpush(info,name)
        goto 9999
     end if

     if(desc_a%ovrlap_elem(1) /= -1) then
        info = 3070
        call psb_errpush(info,name)
        goto 9999
     end if

     ! checking for vectors correctness
     call psb_chkvect(m,ik,size(x),ix,jx,desc_a,info,iix,jjx)
     if (info == 0)&
          & call psb_chkvect(n,ik,size(y),iy,jy,desc_a,info,iiy,jjy)
     if(info /= 0) then
        info=4010
        ch_err='psb_chkvect'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

     if((iix /= 1).or.(iiy /= 1)) then
        ! this case is not yet implemented
        info = 3040
        call psb_errpush(info,name)
        goto 9999
     end if

     xp => x(iix:lldx)
     yp => y(iiy:lldy)

     yp(nrow+1:ncol)=zzero

     !  local Matrix-vector product
     call psb_csmm(alpha,a,xp,beta,yp,info,trans=itrans)

     if(info /= 0) then
        info = 4010
        ch_err='dcsmm'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

     if(idoswap /= 0)&
          & call psi_swaptran(ior(psb_swap_send_,psb_swap_recv_),&
          & zone,yp,desc_a,iwork,info)
     if(info /= 0) then
        info = 4010
        ch_err='PSI_dSwapTran'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

  end if

  if(aliw) deallocate(iwork)
  nullify(iwork)

  call psb_erractionrestore(err_act)

  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_abort_) then
     call psb_error(ictxt)
     return
  end if
  return
end subroutine psb_zspmv
