!!$ 
!!$ 
!!$                    MD2P4
!!$    Multilevel Domain Decomposition Parallel Preconditioner Package for PSBLAS
!!$                      for 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$                       Daniela Di Serafino    II University of Naples
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
subroutine psb_zsp_renum(a,desc_a,blck,p,atmp,info)
  use psb_serial_mod
  use psb_const_mod
  use psb_prec_type
  use psb_descriptor_type
  use psb_spmat_type
  use psb_tools_mod
  use psb_psblas_mod
  use psb_error_mod
  implicit none

  !     .. array Arguments ..                                                     
  type(psb_zspmat_type), intent(in)      :: a,blck
  type(psb_zspmat_type), intent(inout)   :: atmp
  type(psb_zbaseprc_type), intent(inout) :: p
  type(psb_desc_type), intent(in)        :: desc_a
  integer, intent(out)   :: info


  character(len=20)      :: name, ch_err
  integer   istpb, istpe, ifctb, ifcte, err_act, irank, icomm, nztota, nztotb,&
       & nztmp, nzl, nnr, ir, mglob, mtype, n_row, nrow_a,n_col, nhalo,lovr, &
       & ind, iind, pi,nr,ns,i,j,jj,k,kk
  integer ::ictxt,nprow,npcol,me,mycol
  integer, pointer :: itmp(:), itmp2(:)
  complex(kind(1.d0)), pointer :: ztmp(:)
  real(kind(1.d0)) :: t1,t2,t3,t4,t5,t6,mpi_wtime, t7, t8
  external  mpi_wtime

  if (psb_get_errstatus().ne.0) return 
  info=0
  name='apply_renum'
  call psb_erractionsave(err_act)

  ictxt=desc_a%matrix_data(psb_ctxt_)
  call blacs_gridinfo(ictxt,nprow,npcol,me,mycol)

!!!!!!!!!!!!!!!! CHANGE FOR NON-CSR A
  !
  ! Renumbering type: 
  !     1. Global column indices
  !     (2. GPS band reduction disabled for the time being)

  if (p%iprcparm(iren_)==renum_glb_) then 
    atmp%m = a%m + blck%m
    atmp%k = a%k
    atmp%fida='CSR'
    atmp%descra = 'GUN'

    ! This is the renumbering coherent with global indices..
    mglob = desc_a%matrix_data(psb_m_)
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
    t3 = mpi_wtime()

    ! Build  ATMP with new numbering. 

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

      else if (ir <= atmp%m ) then 

        ir = ir - a%m
        nzl = blck%ia2(ir+1) - blck%ia2(ir)
        if (nzl > size(ztmp)) then
          call psb_realloc(nzl,ztmp,info)
          if(info/=0) then
            info=4010
            ch_err='psb_realloc'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        endif
        jj = blck%ia2(ir)
        k=0
        do kk=1, nzl
          if (blck%ia1(jj+kk-1)<=atmp%m) then  
            k = k + 1
            ztmp(k)         = blck%aspk(jj+kk-1)
            atmp%ia1(j+k-1) = p%invperm(blck%ia1(jj+kk-1))
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

    t4 = mpi_wtime()


    deallocate(itmp,itmp2,ztmp)

  else if (p%iprcparm(iren_)==renum_gps_) then 

    atmp%m = a%m + blck%m
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

    if (blck%m>0) then 
      do i=1, blck%m
        atmp%ia2(a%m+i) = nztota+blck%ia2(i)
        do j= blck%ia2(i), blck%ia2(i+1)-1
          atmp%ia1(nztota+j) = blck%ia1(j)
        enddo
      enddo
      atmp%ia2(atmp%m+1) = nztota+blck%ia2(blck%m+1)
    endif
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
    !        call blacs_barrier(ictxt,'All')

    !          write(0,*) size(p%av(u_pr_)%pl),size(p%av(l_pr_)%pr)
    call  gps_reduction(atmp%m,atmp%ia2,atmp%ia1,p%perm,p%invperm,info)
    if(info/=0) then
      info=4010
      ch_err='gps_reduction'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    !      write(0,*) me,' Renumbering: Done GPS'
    call blacs_barrier(ictxt,'All')
    do i=1, atmp%m 
      if (p%perm(i) /= i) then 
        write(0,*) me,' permutation is not identity '
        exit
      endif
    enddo



    do k=1, nnr 
      p%invperm(p%perm(k)) = k
    enddo
    t3 = mpi_wtime()

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

      else if (ir <= atmp%m ) then 

        ir = ir - a%m
        nzl = blck%ia2(ir+1) - blck%ia2(ir)
        if (nzl > size(ztmp)) then
          call psb_realloc(nzl,ztmp,info)
          if(info/=0) then
            info=4010
            ch_err='psb_realloc'
            call psb_errpush(info,name,a_err=ch_err)
            goto 9999
          end if
        endif
        jj = blck%ia2(ir)
        k=0
        do kk=1, nzl
          if (blck%ia1(jj+kk-1)<=atmp%m) then  
            k = k + 1
            ztmp(k)         = blck%aspk(jj+kk-1)
            atmp%ia1(j+k-1) = p%invperm(blck%ia1(jj+kk-1))
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

    t4 = mpi_wtime()



    deallocate(itmp,itmp2,ztmp)

  end if

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


  subroutine gps_reduction(m,ia,ja,perm,iperm,info)
    integer i,j,dgConn,Npnt,m
    integer n,idpth,ideg,ibw2,ipf2
    integer,dimension(:) :: perm,iperm,ia,ja
    integer, intent(out) :: info

    integer,dimension(:,:),allocatable::NDstk
    integer,dimension(:),allocatable::iOld,renum,ndeg,lvl,lvls1,lvls2,ccstor
    !--- Per la common area.

    common /gra/ n,iDpth,iDeg

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
    call reduce(Ndstk,Npnt,iOld,renum,ndeg,lvl,lvls1, lvls2,ccstor,ibw2,ipf2)
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
    if (err_act.eq.act_abort) then
      call psb_error()
      return
    end if
    return

  end subroutine gps_reduction

end subroutine psb_zsp_renum
