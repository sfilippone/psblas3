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
subroutine psb_zsp_renum(a,desc_a,p,atmp,info)
  use psb_base_mod
  use psb_prec_type
  implicit none

  !     .. array Arguments ..                                                     
  type(psb_zspmat_type), intent(in)      :: a
  type(psb_zspmat_type), intent(inout)   :: atmp
  type(psb_zprec_type), intent(inout) :: p
  type(psb_desc_type), intent(in)        :: desc_a
  integer, intent(out)   :: info


  character(len=20)      :: name, ch_err
  integer   nztota, nztotb, nztmp, nzl, nnr, ir, mglob, mtype, n_row, &
       & nrow_a,n_col, nhalo,lovr,  ind, iind, pi,nr,ns,i,j,jj,k,kk
  integer ::ictxt,np,me, err_act
  integer, allocatable  :: itmp(:), itmp2(:)
  complex(kind(1.d0)), allocatable  :: ztmp(:)
  real(kind(1.d0)) :: t1,t2,t3,t4,t5,t6, t7, t8

  if (psb_get_errstatus().ne.0) return 
  info=0
  name='apply_renum'
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)

!!!!!!!!!!!!!!!! CHANGE FOR NON-CSR A
  !
  ! Renumbering type: 
  !     1. Global column indices
  !     (2. GPS band reduction disabled for the time being)

  if (p%iprcparm(iren_)==renum_glb_) then 
    atmp%m = a%m 
    atmp%k = a%k
    atmp%fida='CSR'
    atmp%descra = 'GUN'

    ! This is the renumbering coherent with global indices..
    mglob = psb_cd_get_global_rows(desc_a)
    !
    !  Remember: we have switched IA1=COLS and IA2=ROWS
    !  Now identify the set of distinct local column indices
    !

    nnr = p%desc_data%matrix_data(psb_n_row_)
    allocate(p%perm(nnr),p%invperm(nnr),itmp2(nnr),stat=info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    do k=1,nnr
      itmp2(k) = p%desc_data%loc_to_glob(k)
    enddo
    !
    !  We want:  NEW(I) = OLD(PERM(I))
    ! 
    call isrx(nnr,itmp2,p%perm)

    do k=1, nnr 
      p%invperm(p%perm(k)) = k
    enddo
    t3 = psb_wtime()

    ! Build  ATMP with new numbering. 
    nztmp=size(atmp%aspk)
    allocate(itmp(max(8,atmp%m+2,nztmp+2)),ztmp(atmp%m),stat=info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    j = 1
    atmp%ia2(1) = 1
    do i=1, atmp%m
      ir = p%perm(i)

      if (ir <= a%m ) then

        nzl = a%ia2(ir+1) - a%ia2(ir)
        if (nzl > size(ztmp)) then
          call psb_realloc(nzl,ztmp,info)
          if(info/=0) then
            info=4010
            ch_err='psb_realloc'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        endif
        jj = a%ia2(ir)
        k=0
        do kk=1, nzl
          if (a%ia1(jj+kk-1)<=atmp%m) then  
            k = k + 1
            ztmp(k) = a%aspk(jj+kk-1)
            atmp%ia1(j+k-1) = p%invperm(a%ia1(jj+kk-1))
          endif
        enddo
        call isrx(k,atmp%ia1(j:j+k-1),itmp2)
        do kk=1,k
          atmp%aspk(j+kk-1) = ztmp(itmp2(kk))
        enddo

      else
        write(0,*) 'Row index error 1 :',i,ir
      endif

      j = j + k
      atmp%ia2(i+1) = j

    enddo

    t4 = psb_wtime()


    deallocate(itmp,itmp2,ztmp)

  else if (p%iprcparm(iren_)==renum_gps_) then 

    atmp%m = a%m
    atmp%k = a%k
    atmp%fida='CSR'
    atmp%descra = 'GUN'
    do i=1, a%m
      atmp%ia2(i) = a%ia2(i)
      do j= a%ia2(i), a%ia2(i+1)-1
        atmp%ia1(j) = a%ia1(j)
      enddo
    enddo
    atmp%ia2(a%m+1) = a%ia2(a%m+1)
    nztmp = atmp%ia2(atmp%m+1) - 1


    ! This is a renumbering with Gibbs-Poole-Stockmeyer 
    ! band reduction. Switched off for now. To be fixed,
    ! gps_reduction should get p%perm. 

    !          write(0,*) me,' Renumbering: realloc perms',atmp%m
    call psb_realloc(atmp%m,p%perm,info)
    if(info/=0) then
      info=4010
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_realloc(atmp%m,p%invperm,info)
    if(info/=0) then
      info=4010
      ch_err='psb_realloc'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    allocate(itmp(max(8,atmp%m+2,nztmp+2)),stat=info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    itmp(1:8) = 0
    !          write(0,*) me,' Renumbering: Calling Metis'

    !          write(0,*) size(p%av(u_pr_)%pl),size(p%av(l_pr_)%pr)
    call  gps_reduction(atmp%m,atmp%ia2,atmp%ia1,p%perm,p%invperm,info)
    if(info/=0) then
      info=4010
      ch_err='gps_reduction'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    !      write(0,*) me,' Renumbering: Done GPS'
    !    call psb_barrier(ictxt)
    do i=1, atmp%m 
      if (p%perm(i) /= i) then 
        write(0,*) me,' permutation is not identity '
        exit
      endif
    enddo



    do k=1, nnr 
      p%invperm(p%perm(k)) = k
    enddo
    t3 = psb_wtime()

    ! Build  ATMP with new numbering. 

    allocate(itmp2(max(8,atmp%m+2,nztmp+2)),ztmp(atmp%m),stat=info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    j = 1
    atmp%ia2(1) = 1
    do i=1, atmp%m
      ir = p%perm(i)

      if (ir <= a%m ) then

        nzl = a%ia2(ir+1) - a%ia2(ir)
        if (nzl > size(ztmp)) then
          call psb_realloc(nzl,ztmp,info)
          if(info/=0) then
            info=4010
            ch_err='psb_realloc'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        endif
        jj = a%ia2(ir)
        k=0
        do kk=1, nzl
          if (a%ia1(jj+kk-1)<=atmp%m) then  
            k = k + 1
            ztmp(k) = a%aspk(jj+kk-1)
            atmp%ia1(j+k-1) = p%invperm(a%ia1(jj+kk-1))
          endif
        enddo
        call isrx(k,atmp%ia1(j:j+k-1),itmp2)
        do kk=1,k
          atmp%aspk(j+kk-1) = ztmp(itmp2(kk))
        enddo

      else
        write(0,*) 'Row index error 1 :',i,ir
      endif

      j = j + k
      atmp%ia2(i+1) = j

    enddo

    t4 = psb_wtime()



    deallocate(itmp,itmp2,ztmp)

  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

contains


  subroutine gps_reduction(m,ia,ja,perm,iperm,info)
    integer i,j,dgConn,Npnt,m
    integer n,idpth,ideg,ibw2,ipf2
    integer,dimension(:) :: perm,iperm,ia,ja
    integer, intent(out) :: info

    integer,dimension(:,:),allocatable::NDstk
    integer,dimension(:),allocatable::iOld,renum,ndeg,lvl,lvls1,lvls2,ccstor
    !--- Per la common area.

    character(len=20)      :: name, ch_err

    if(psb_get_errstatus().ne.0) return 
    info=0
    name='gps_reduction'
    call psb_erractionsave(err_act)


    !--- Calcolo il massimo grado di connettivita'.
    npnt = m
    write(6,*) ' GPS su ',npnt
    dgConn=0
    do i=1,m
      dgconn = max(dgconn,(ia(i+1)-ia(i)))
    enddo
    !--- Il max valore di connettivita' e "dgConn"

    !--- Valori della common
    n=Npnt       !--- Numero di righe
    iDeg=dgConn  !--- Massima connettivita'
    !    iDpth=       !--- Numero di livelli non serve settarlo

    allocate(NDstk(Npnt,dgConn),stat=info)
    if (info/=0) then 
      info=4000
      call psb_errpush(info,name)
      goto 9999
    else
      write(0,*) 'gps_reduction first alloc OK'
    endif
    allocate(iOld(Npnt),renum(Npnt+1),ndeg(Npnt),lvl(Npnt),lvls1(Npnt),&
         &lvls2(Npnt),ccstor(Npnt),stat=info)
    if (info/=0) then 
      info=4000
      call psb_errpush(info,name)
      goto 9999
    else
      write(0,*) 'gps_reduction 2nd alloc OK'
    endif

    !--- Prepariamo il grafo della matrice
    Ndstk(:,:)=0
    do i=1,Npnt
      k=0
      do j = ia(i),ia(i+1) - 1 
        if ((1<=ja(j)).and.( ja( j ) /= i ).and.(ja(j)<=npnt)) then
          k = k+1
          Ndstk(i,k)=ja(j)
        endif
      enddo
      ndeg(i)=k
    enddo

    !--- Numerazione.
    do i=1,Npnt
      iOld(i)=i
    enddo
    write(0,*) 'gps_red : Preparation done'
    !--- 
    !--- Chiamiamo funzione reduce.
    call psb_gps_reduce(Ndstk,Npnt,iOld,renum,ndeg,lvl,lvls1, lvls2,ccstor,&
         & ibw2,ipf2,n,idpth,ideg)
    write(0,*) 'gps_red : Done reduce'
    !--- Permutazione
    perm(1:Npnt)=renum(1:Npnt)
    !--- Inversa permutazione
    do i=1,Npnt
      iperm(perm(i))=i
    enddo
    !--- Puliamo tutto.
    deallocate(NDstk,iOld,renum,ndeg,lvl,lvls1,lvls2,ccstor)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.psb_act_abort_) then
      call psb_error()
      return
    end if
    return

  end subroutine gps_reduction

end subroutine psb_zsp_renum
