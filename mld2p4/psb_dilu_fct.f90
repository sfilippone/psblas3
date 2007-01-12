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
subroutine psb_dilu_fct(a,l,u,d,info,blck)
  
  !
  ! This routine copies and factors "on the fly" from A and BLCK
  ! into L/D/U. 
  !
  !
  use psb_base_mod
  implicit none
  !     .. Scalar Arguments ..
  integer, intent(out)                ::     info
  !     .. Array Arguments ..
  type(psb_dspmat_type),intent(in)    :: a
  type(psb_dspmat_type),intent(inout) :: l,u
  type(psb_dspmat_type),intent(in), optional, target :: blck
  real(kind(1.d0)), intent(inout)     ::  d(:)
  !     .. Local Scalars ..
  real(kind(1.d0)) ::  dia, temp
  integer   ::  i, j, jj, k, kk, l1, l2, ll, low1, low2,m,ma,err_act
  
  type(psb_dspmat_type), pointer  :: blck_
  character(len=20)   :: name, ch_err
  logical, parameter :: debug=.false.
  name='psb_dcsrlu'
  info = 0
  call psb_erractionsave(err_act)
  !     .. Executable Statements ..
  !

  if (present(blck)) then 
    blck_ => blck
  else
    allocate(blck_,stat=info) 
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    call psb_nullify_sp(blck_)  ! Why do we need this? Who knows.... 
    call psb_sp_all(0,0,blck_,1,info)
    if(info.ne.0) then
       info=4010
       ch_err='psb_sp_all'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if

    blck_%m=0
  endif

!!$  write(0,*) 'ilu_fct: ',size(l%ia2),size(u%ia2),a%m,blck_%m
  call psb_dilu_fctint(m,a%m,a,blck_%m,blck_,&
       & d,l%aspk,l%ia1,l%ia2,u%aspk,u%ia1,u%ia2,l1,l2,info)
  if(info.ne.0) then
     info=4010
     ch_err='psb_dilu_fctint'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  
  l%infoa(1) = l1
  l%fida     = 'CSR'
  l%descra   = 'TLU'
  u%infoa(1) = l2
  u%fida     = 'CSR'
  u%descra   = 'TUU'
  l%m = m
  l%k = m
  u%m = m
  u%k = m
  if (present(blck)) then 
    blck_ => null() 
  else
    call psb_sp_free(blck_,info)
    if(info.ne.0) then
       info=4010
       ch_err='psb_sp_free'
       call psb_errpush(info,name,a_err=ch_err)
       goto 9999
    end if
    deallocate(blck_) 
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return

contains
  subroutine psb_dilu_fctint(m,ma,a,mb,b,&
       & d,laspk,lia1,lia2,uaspk,uia1,uia2,l1,l2,info)
    implicit none 

    type(psb_dspmat_type)          :: a,b
    integer                        :: m,ma,mb,l1,l2,info
    integer, dimension(*)          :: lia1,lia2,uia1,uia2
    real(kind(1.d0)), dimension(*) :: laspk,uaspk,d

    integer :: i,j,k,l,low1,low2,kk,jj,ll, irb, ktrw,err_act
    real(kind(1.d0)) :: dia,temp
    integer, parameter :: nrb=16
    logical,parameter  :: debug=.false.
    type(psb_dspmat_type) :: trw
    integer             :: int_err(5) 
    character(len=20)   :: name, ch_err

    name='psb_dilu_fctint'
    if(psb_get_errstatus().ne.0) return 
    info=0
    call psb_erractionsave(err_act)
    call psb_nullify_sp(trw)
    trw%m=0
    trw%k=0
    if(debug) write(0,*)'LUINT Allocating TRW'
    call psb_sp_all(trw,1,info)
    if(info.ne.0) then
      info=4010
      ch_err='psb_sp_all'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    if(debug) write(0,*)'LUINT Done  Allocating TRW'
    lia2(1) = 1
    uia2(1) = 1
    l1=0
    l2=0
    m = ma+mb
    if(debug) write(0,*)'In DCSRLU Begin cycle',m,ma,mb

    do i = 1, ma
      if(debug) write(0,*)'LUINT: Loop index ',i,ma
      d(i) = 0.d0

      !
      ! Here we take a fast shortcut if possible, otherwise 
      ! use spgtblk, slower but able (in principle) to handle 
      ! anything. 
      !
      if (a%fida=='CSR') then 
        do j = a%ia2(i), a%ia2(i+1) - 1
          k = a%ia1(j)
          !           write(0,*)'KKKKK',k
          if ((k < i).and.(k >= 1)) then
            l1 = l1 + 1
            laspk(l1) = a%aspk(j)
            lia1(l1) = k
          else if (k == i) then
            d(i) = a%aspk(j)
          else if ((k > i).and.(k <= m)) then
            l2 = l2 + 1
            uaspk(l2) = a%aspk(j)
            uia1(l2) = k
          end if
        enddo

      else

        if ((mod(i,nrb) == 1).or.(nrb==1)) then 
          irb = min(ma-i+1,nrb)
          call psb_sp_getblk(i,a,trw,info,lrw=i+irb-1)
          if(info.ne.0) then
            info=4010
            ch_err='psb_sp_getblk'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
          ktrw=1
        end if

        do 
          if (ktrw > trw%infoa(psb_nnz_)) exit
          if (trw%ia1(ktrw) > i) exit
          k = trw%ia2(ktrw)
          if ((k < i).and.(k >= 1)) then
            l1 = l1 + 1
            laspk(l1) = trw%aspk(ktrw)
            lia1(l1) = k
          else if (k == i) then
            d(i) = trw%aspk(ktrw)
          else if ((k > i).and.(k <= m)) then
            l2 = l2 + 1
            uaspk(l2) = trw%aspk(ktrw)
            uia1(l2) = k
          end if
          ktrw = ktrw + 1
        enddo

      end if

!!$

      lia2(i+1) = l1 + 1
      uia2(i+1) = l2 + 1

      dia = d(i)
      do kk = lia2(i), lia2(i+1) - 1
        !     
        !     compute element alo(i,k) of incomplete factorization
        !     
        temp = laspk(kk)
        k = lia1(kk)
        laspk(kk) = temp*d(k)
        !     update the rest of row i using alo(i,k)
        low1 = kk + 1
        low2 = uia2(i)
        updateloop: do  jj = uia2(k), uia2(k+1) - 1
          j = uia1(jj)
          !
          if (j < i) then
            !     search alo(i,*) for matching index J
            do  ll = low1, lia2(i+1) - 1
              l = lia1(ll)
              if (l > j) then
                low1 = ll
                exit
              else if (l == j) then
                laspk(ll) = laspk(ll) - temp*uaspk(jj)
                low1 = ll + 1
                cycle updateloop
              end if
            enddo
            !     
          else if (j == i) then
            !    j=i  update diagonal
            !    write(0,*)'aggiorno dia',dia,'temp',temp,'jj',jj,'u%aspk',uaspk(jj)
            dia = dia - temp*uaspk(jj)
            !    write(0,*)'dia',dia,'temp',temp,'jj',jj,'aspk',uaspk(jj)
            cycle updateloop
            !     
          else if (j > i) then
            !     search aup(i,*) for matching index j
            do ll = low2, uia2(i+1) - 1
              l = uia1(ll)
              if (l > j) then
                low2 = ll
                exit
              else if (l == j) then
                uaspk(ll) = uaspk(ll) - temp*uaspk(jj)
                low2 = ll + 1
                cycle updateloop
              end if
            enddo
          end if
          !     
          !     for milu al=1.;  for ilu al=0.
          !     al = 1.d0
          !     dia = dia - al*temp*aup(jj)
        enddo updateloop
      enddo
      !    
      !     
      !     Non singularity
      !     
      if (dabs(dia) < epstol) then
        !
        !     Pivot too small: unstable factorization
        !     
        info = 2
        int_err(1) = i
        write(ch_err,'(g20.10)') dia
        call psb_errpush(info,name,i_err=int_err,a_err=ch_err)
        goto 9999
      else
        dia = 1.d0/dia
      end if
      d(i) = dia
      ! write(6,*)'diag(',i,')=',d(i)
      !     Scale row i of upper triangle
      do  kk = uia2(i), uia2(i+1) - 1
        uaspk(kk) = uaspk(kk)*dia
      enddo
    enddo

    do i = ma+1, m
      d(i) = 0.d0


      if (b%fida=='CSR') then 

        do j = b%ia2(i-ma), b%ia2(i-ma+1) - 1
          k = b%ia1(j)
          !           if (me.eq.2)  write(0,*)'ecco k=',k
          if ((k < i).and.(k >= 1)) then
            l1 = l1 + 1
            laspk(l1) = b%aspk(j)
            lia1(l1) = k
            !              if(me.eq.2) write(0,*)'scrivo l'
          else if (k == i) then
            d(i) = b%aspk(j)
          else if ((k > i).and.(k <= m)) then
            l2 = l2 + 1
            uaspk(l2) = b%aspk(j)
            !              write(0,*)'KKKKK',k
            uia1(l2) = k
          end if
        enddo

      else

        if ((mod((i-ma),nrb) == 1).or.(nrb==1)) then 
          irb = min(m-i+1,nrb)
          call psb_sp_getblk(i-ma,b,trw,info,lrw=i-ma+irb-1)
          if(info.ne.0) then
            info=4010
            ch_err='psb_sp_getblk'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
          ktrw=1
        end if

        do 
          if (ktrw > trw%infoa(psb_nnz_)) exit
          if (trw%ia1(ktrw) > i) exit
          k = trw%ia2(ktrw)
          ! write(0,*)'KKKKK',k
          if ((k < i).and.(k >= 1)) then
            l1 = l1 + 1
            laspk(l1) = trw%aspk(ktrw)
            lia1(l1) = k
          else if (k == i) then
            d(i) = trw%aspk(ktrw)
          else if ((k > i).and.(k <= m)) then
            l2 = l2 + 1
            uaspk(l2) = trw%aspk(ktrw)
            uia1(l2) = k
          end if
          ktrw = ktrw + 1
        enddo

      endif


      lia2(i+1) = l1 + 1
      uia2(i+1) = l2 + 1

      dia = d(i)
      do kk = lia2(i), lia2(i+1) - 1
        !     
        !     compute element alo(i,k) of incomplete factorization
        !     
        temp = laspk(kk)
        k = lia1(kk)
        laspk(kk) = temp*d(k)
        !     update the rest of row i using alo(i,k)
        low1 = kk + 1
        low2 = uia2(i)
        updateloopb: do  jj = uia2(k), uia2(k+1) - 1
          j = uia1(jj)
          !
          if (j < i) then
            !     search alo(i,*) for matching index J
            do  ll = low1, lia2(i+1) - 1
              l = lia1(ll)
              if (l > j) then
                low1 = ll
                exit
              else if (l == j) then
                laspk(ll) = laspk(ll) - temp*uaspk(jj)
                low1 = ll + 1
                cycle updateloopb
              end if
            enddo
            !     
          else if (j == i) then
            !   j=i  update diagonal
            dia = dia - temp*uaspk(jj)
            cycle updateloopb
            !     
          else if (j > i) then
            !    search aup(i,*) for matching index j
            do ll = low2, uia2(i+1) - 1
              l = uia1(ll)
              if (l > j) then
                low2 = ll
                exit
              else if (l == j) then
                uaspk(ll) = uaspk(ll) - temp*uaspk(jj)
                low2 = ll + 1
                cycle updateloopb
              end if
            enddo
          end if
          !     
          !     for milu al=1.;  for ilu al=0.
          !     al = 1.d0
          !     dia = dia - al*temp*aup(jj)
        enddo updateloopb
      enddo
      !     
      !     
      !     Non singularity
      !     
      if (dabs(dia) < epstol) then
        !
        !     Pivot too small: unstable factorization
        !     
        int_err(1) = i
        write(ch_err,'(g20.10)') dia
        info = 2
        call psb_errpush(info,name,i_err=int_err,a_err=ch_err)
        goto 9999
      else
        dia = 1.d0/dia
      end if
      d(i) = dia
      !     Scale row i of upper triangle
      do  kk = uia2(i), uia2(i+1) - 1
        uaspk(kk) = uaspk(kk)*dia
      enddo
    enddo

    call psb_sp_free(trw,info)
    if(info.ne.0) then
      info=4010
      ch_err='psb_sp_free'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    if(debug) write(0,*)'Leaving ilu_fct'

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.act_abort) then
      call psb_error()
      return
    end if
    return
  end subroutine psb_dilu_fctint
end subroutine psb_dilu_fct
