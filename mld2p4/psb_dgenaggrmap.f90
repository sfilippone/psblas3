!!$ 
!!$ 
!!$                    MD2P4
!!$    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
!!$                      for 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$                       Daniela di Serafino    Second University of Naples
!!$                       Pasqua D'Ambra         ICAR-CNR                      
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the MD2P4 group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MD2P4 GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
subroutine psb_dgenaggrmap(aggr_type,a,desc_a,nlaggr,ilaggr,info)
  use psb_base_mod
  use psb_prec_type
  implicit none
  integer, intent(in)               :: aggr_type
  type(psb_dspmat_type), intent(in) :: a
  type(psb_desc_type), intent(in)   :: desc_a
  integer, allocatable              :: ilaggr(:),nlaggr(:)
  integer, intent(out)              :: info
  ! Locals 
  integer, allocatable  :: ils(:), neigh(:)
  integer :: icnt,nlp,k,n,ia,isz,nr, naggr,i,j,m

  logical :: recovery
  logical, parameter :: debug=.false.
  integer ::ictxt,np,me,err_act
  integer :: nrow, ncol, n_ne
  integer, parameter :: one=1, two=2
  character(len=20)   :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name = 'psb_bldaggrmat'
  call psb_erractionsave(err_act)
  !
  ! Note. At the time being we are ignoring aggr_type 
  ! so that we only have local decoupled aggregation. This might 
  ! change in the future. 
  !
  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt,me,np)
  nrow  = psb_cd_get_local_rows(desc_a)
  ncol  = psb_cd_get_local_cols(desc_a)

  nr = a%m
  allocate(ilaggr(nr),neigh(nr),stat=info)
  if(info.ne.0) then
    info=4000
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  do i=1, nr
    ilaggr(i) = -(nr+1)
  end do
  ! Note: -(nr+1)  Untouched as yet
  !       -i    1<=i<=nr  Adjacent to aggregate i
  !        i    1<=i<=nr  Belonging to aggregate i

  !
  ! Phase one: group nodes together. 
  ! Very simple minded strategy. 
  ! 
  naggr = 0
  nlp   = 0
  do
    icnt = 0
    do i=1, nr 
      if (ilaggr(i) == -(nr+1)) then 
        !
        ! 1. Untouched nodes are marked >0 together 
        !    with their neighbours
        !
        icnt = icnt + 1 
        naggr = naggr + 1 
        ilaggr(i) = naggr

        call psb_neigh(a,i,neigh,n_ne,info,lev=one)
        if (info/=0) then 
          info=4010
          ch_err='psb_neigh'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if
        do k=1, n_ne
          j = neigh(k)
          if ((1<=j).and.(j<=nr)) then 
            ilaggr(j) = naggr
!!$              if (ilaggr(j) < 0) ilaggr(j) = naggr
!!$              if (ilaggr(j) == -(nr+1)) ilaggr(j) = naggr
          endif
        enddo
        !
        ! 2. Untouched neighbours of these nodes are marked <0.
        !
        call psb_neigh(a,i,neigh,n_ne,info,lev=two)
        if (info/=0) then 
          info=4010
          ch_err='psb_neigh'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
        end if

        do n = 1, n_ne
          m = neigh(n)
          if ((1<=m).and.(m<=nr)) then
            if (ilaggr(m) == -(nr+1)) ilaggr(m) = -naggr
          endif
        enddo
      endif
    enddo
    nlp = nlp + 1 
    if (icnt == 0) exit 
  enddo
  if (debug) then 
    write(0,*) 'Check 1:',count(ilaggr == -(nr+1)),(a%ia1(i),i=a%ia2(1),a%ia2(2)-1)
  end if

  !
  ! Phase two: sweep over leftovers. 
  !
  allocate(ils(naggr+10),stat=info) 
  if(info.ne.0) then
    info=4000
    call psb_errpush(info,name)
    goto 9999
  end if

  do i=1, size(ils)
    ils(i) = 0
  end do
  do i=1, nr 
    n = ilaggr(i)
    if (n>0) then 
      if (n>naggr) then 
        write(0,*) 'loc_Aggregate: n > naggr 1 ? ',n,naggr
      else
        ils(n) = ils(n) + 1 
      end if

    end if
  end do
  if (debug) then 
    write(0,*) 'Phase 1: number of aggregates ',naggr
    write(0,*) 'Phase 1: nodes aggregated     ',sum(ils)
  end if

  recovery=.false.
  do i=1, nr
    if (ilaggr(i) < 0) then 
      !
      ! Now some silly rule to break ties:
      ! Group with smallest adjacent aggregate. 
      !
      isz = nr+1
      ia  = -1

      call psb_neigh(a,i,neigh,n_ne,info,lev=one)
      if (info/=0) then 
        info=4010
        ch_err='psb_neigh'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
      end if

      do j=1, n_ne
        k = neigh(j)
        if ((1<=k).and.(k<=nr))  then 
          n = ilaggr(k) 
          if (n>0) then 
            if (n>naggr) then 
              write(0,*) 'loc_Aggregate: n > naggr 2? ',n,naggr
            end if

            if (ils(n) < isz) then 
              ia = n
              isz = ils(n)
            endif
          endif
        endif
      enddo
      if (ia == -1) then 
        if (ilaggr(i) > -(nr+1)) then
          ilaggr(i) = abs(ilaggr(i))
          if (ilaggr(I)>naggr) then 
            write(0,*) 'loc_Aggregate: n > naggr 3? ',ilaggr(i),naggr
          end if
          ils(ilaggr(i)) = ils(ilaggr(i)) + 1
          !
          ! This might happen if the pattern is non symmetric. 
          ! Need a better handling. 
          !
          recovery = .true.
        else
          write(0,*) 'Unrecoverable error !!',ilaggr(i), nr
        endif
      else
        ilaggr(i) = ia
        if (ia>naggr) then 
          write(0,*) 'loc_Aggregate: n > naggr 4? ',ia,naggr
        end if

        ils(ia)  = ils(ia) + 1
      endif
    end if
  enddo
  if (recovery) then 
    write(0,*) 'Had to recover from strange situation in loc_aggregate.'
    write(0,*) 'Perhaps an unsymmetric pattern?'
  endif
  if (debug) then 
    write(0,*) 'Phase 2: number of aggregates ',naggr
    write(0,*) 'Phase 2: nodes aggregated     ',sum(ils) 
    do i=1, naggr 
      write(*,*) 'Size of aggregate ',i,' :',count(ilaggr==i), ils(i)
    enddo
    write(*,*) maxval(ils(1:naggr))
    write(*,*) 'Leftovers ',count(ilaggr<0), '  in ',nlp,' loops'
  end if

!!$    write(0,*) 'desc_a loc_aggr 4 : ', desc_a%matrix_data(m_)
  if (count(ilaggr<0) >0) then 
    write(0,*) 'Fatal error: some leftovers!!!'
  endif

  deallocate(ils,neigh,stat=info)
  if (info/=0) then 
    info=4000
    call psb_errpush(info,name)
    goto 9999
  end if

  if (nrow /= size(ilaggr)) then 
    write(0,*) 'SOmething wrong ilaggr ',nrow,size(ilaggr)
  endif
  call psb_realloc(ncol,ilaggr,info)
  if (info/=0) then 
    info=4010
    ch_err='psb_realloc'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  allocate(nlaggr(np),stat=info)
  if (info/=0) then 
    info=4000
    call psb_errpush(info,name)
    goto 9999
  end if

  nlaggr(:) = 0
  nlaggr(me+1) = naggr
  call psb_sum(ictxt,nlaggr(1:np))

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error()
    return
  end if
  return

end subroutine psb_dgenaggrmap
