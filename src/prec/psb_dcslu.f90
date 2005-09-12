!*****************************************************************************
!*                                                                           *
!* This is where the action takes place. As you may notice, the only         *
!* piece that's really enabled is that for CSR. This is to be fixed.         *
!* CSRSETUP does the setup: building the prec descriptor plus retrieving     *
!*                           matrix rows if needed                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!* some open code does the renumbering                                       *
!*                                                                           *
!* DSPLU90  does the actual factorization.                                  *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*****************************************************************************
subroutine psb_dcslu(a,desc_a,p,upd,info)
  use psb_serial_mod
  use psb_const_mod
  use psb_prec_type
  use psb_descriptor_type
  use psb_spmat_type
  use psb_tools_mod
  use psb_psblas_mod
  use psb_error_mod
  implicit none
  !                                                                               
  !     .. Scalar Arguments ..                                                    
  integer, intent(out)                      :: info
  !     .. array Arguments ..                                                     
  type(psb_dspmat_type), intent(in), target :: a
  type(psb_dbase_prec), intent(inout)       :: p
  type(psb_desc_type), intent(in)           :: desc_a
  character, intent(in)                     :: upd

  !     .. Local Scalars ..                                                       
  integer  ::    i, j, jj, k, kk, m
  integer  ::    int_err(5)
  character ::        trans, unitd
  type(psb_dspmat_type) :: blck, atmp
  real(kind(1.d0)) :: t1,t2,t3,t4,t5,t6,mpi_wtime, t7, t8
  integer, pointer :: itmp(:), itmp2(:)
  real(kind(1.d0)), pointer :: rtmp(:)
  external  mpi_wtime
  logical, parameter :: debugprt=.false., debug=.false., aggr_dump=.false.
  integer   istpb, istpe, ifctb, ifcte, err_act, irank, icomm, nztota, nztotb,&
       & nztmp, nzl, nnr, ir, mglob, mtype, n_row, nrow_a,n_col, nhalo,lovr
  integer ::icontxt,nprow,npcol,me,mycol
  character(len=20)      :: name, ch_err

  interface
     subroutine psb_dsplu(a,l,u,d,info,blck)
       use psb_spmat_type
       integer, intent(out)                ::     info
       type(psb_dspmat_type),intent(in)    :: a
       type(psb_dspmat_type),intent(inout) :: l,u
       type(psb_dspmat_type),intent(in), optional, target :: blck
       real(kind(1.d0)), intent(inout)     ::  d(:)
     end subroutine psb_dsplu
  end interface

  info=0
  name='psb_dcslu'
  call psb_erractionsave(err_act)

  icontxt=desc_a%matrix_data(psb_ctxt_)
  call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)

  m = a%m
  if (m < 0) then
     info = 10
     int_err(1) = 1
     int_err(2) = m
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  endif
  trans = 'N'
  unitd = 'U'
  if (p%iprcparm(n_ovr_) < 0) then
     info = 11
     int_err(1) = 1
     int_err(2) = p%iprcparm(n_ovr_)
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  endif

  !  call blacs_gridinfo(icontxt,nprow,npcol,me,mycol)


  icontxt=desc_a%matrix_data(psb_ctxt_)
  call psb_nullify_sp(blck)
  t1= mpi_wtime()

  if(debug) write(0,*)me,': calling psb_csrsetup'
  call psb_csrsetup(p%iprcparm(p_type_),p%iprcparm(n_ovr_),a,&
       & blck,desc_a,upd,p%desc_data,info)
  if(info/=0) then
     info=4010
     ch_err='psb_csrsetup'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  t2= mpi_wtime()
  if(debug) write(0,*)me,': out of psb_csrsetup'

  if (associated(p%av)) then 
     if (size(p%av) < bp_ilu_avsz) then 
        do k=1,size(p%av)
           call psb_spfree(p%av(k),info)
        end do
        deallocate(p%av)
        p%av => null()
     endif
  endif

  if (.not.associated(p%av)) then 
     allocate(p%av(bp_ilu_avsz))
  endif
  do k=1,size(p%av)
     call psb_nullify_sp(p%av(k))
  end do
  nrow_a = desc_a%matrix_data(psb_n_row_)
  call psb_spinfo(psb_nztotreq_,a,nztota,info)
  if(info/=0) then
     info=4010
     ch_err='psb_spinfo'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if
  n_col  = desc_a%matrix_data(psb_n_col_)
  nhalo  = n_col-nrow_a
  n_row  = p%desc_data%matrix_data(psb_n_row_)
  lovr   = ((nztota+nrow_a-1)/nrow_a)*nhalo*p%iprcparm(n_ovr_)
  p%av(l_pr_)%m  = n_row
  p%av(l_pr_)%k  = n_row
  p%av(u_pr_)%m  = n_row
  p%av(u_pr_)%k  = n_row
  call psb_spall(n_row,n_row,p%av(l_pr_),nztota+lovr,info)
  call psb_spall(n_row,n_row,p%av(u_pr_),nztota+lovr,info)
  if(info/=0) then
     info=4010
     ch_err='psb_spall'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if (associated(p%d)) then 
     if (size(p%d) < n_row) then 
        deallocate(p%d)
     endif
  endif
  if (.not.associated(p%d)) then 
     allocate(p%d(n_row),stat=info)
  endif


  if (debug) then 
     write(0,*) me,'Done psb_csrsetup'
     call blacs_barrier(icontxt,'All')
  endif

!  if (me==0) write(0,*) 'setup time',t2-t1, blck%fida, p%iprcparm(p_type_),blck%m,upd
!    write(0,'(i3,1x,a,4(1x,g14.5))') me,' setup time',t2-t1

  if (p%iprcparm(iren_) > 0) then 

     !
     ! Here we allocate a full copy to hold local A and received BLK
     !

     call psb_spinfo(psb_nztotreq_,a,nztota,info)
     call psb_spinfo(psb_nztotreq_,blck,nztotb,info)
     call psb_spall(atmp,nztota+nztotb,info)
     if(info/=0) then
        info=4011
        call psb_errpush(info,name)
        goto 9999
     end if


!      write(0,*) 'DCSLU ',nztota,nztotb,a%m


     call apply_renum(info)
     if(info/=0) then
        info=4010
        ch_err='apply_ernum'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if

     t3 = mpi_wtime()
     if (debugprt) then 
        open(40+me) 
        call psb_csprt(40+me,atmp,head='% Local matrix')
        close(40+me)
     endif
     if (debug) write(0,*) me,' Factoring rows ',&
          &atmp%m,a%m,blck%m,atmp%ia2(atmp%m+1)-1

     !
     ! Ok, factor the matrix.  
     !
     t5 = mpi_wtime()
     blck%m=0
     call psb_dsplu(atmp,p%av(l_pr_),p%av(u_pr_),p%d,info,blck=blck)
     if(info/=0) then
        info=4010
        ch_err='psb_dsplu'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if
     
     call psb_spfree(atmp,info) 
     if(info/=0) then
        info=4010
        ch_err='psb_spfree'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if


  else if (p%iprcparm(iren_) == 0) then
     t3 = mpi_wtime()
!    ierr = MPE_Log_event( ifctb, 0, "st SIMPLE" )
     ! This is where we have mo renumbering, thus no need 
     ! for ATMP

     if (debugprt) then 
        open(40+me)
        call psb_csprt(40+me,a,iv=p%desc_data%loc_to_glob,&
             &    head='% Local matrix')
        if (p%iprcparm(p_type_)==asm_) then 
           call psb_csprt(40+me,blck,iv=p%desc_data%loc_to_glob,&
                &  irs=a%m,head='% Received rows')
        endif
        close(40+me)
     endif



     t5= mpi_wtime()
     if (debug) write(0,*) me,' Going for dsplu'
     call psb_dsplu(a,p%av(l_pr_),p%av(u_pr_),p%d,info,blck=blck)
     if(info/=0) then
        info=4010
        ch_err='psb_dsplu'
        call psb_errpush(info,name,a_err=ch_err)
        goto 9999
     end if
     if (debug) write(0,*) me,' Done dsplu'
  endif


  if (debugprt) then 
     !
     ! Print out the factors on file.
     !
     open(80+me)

     call psb_csprt(80+me,p%av(l_pr_),head='% Local L factor')
     write(80+me,*) '% Diagonal: ',p%av(l_pr_)%m
     do i=1,p%av(l_pr_)%m
        write(80+me,*) i,i,p%d(i)
     enddo
     call psb_csprt(80+me,p%av(u_pr_),head='% Local U factor')

     close(80+me)
  endif


!    ierr = MPE_Log_event( ifcte, 0, "st SIMPLE" )
  t6 = mpi_wtime()
!
!    write(0,'(i3,1x,a,3(1x,g18.9))') me,'renum/factor time',t3-t2,t6-t5
!    if (me==0) write(0,'(a,3(1x,g18.9))') 'renum/factor time',t3-t2,t6-t5

  call psb_spfree(blck,info)
  if(info/=0) then
     info=4010
     ch_err='psb_spfree'
     call psb_errpush(info,name,a_err=ch_err)
     goto 9999
  end if

  if (debug) write(0,*) me,'End of cslu'

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

  subroutine apply_renum(info)

    integer, intent(out)   :: info
    character(len=20)      :: name, ch_err

    info=0
    name='apply_renum'
    call psb_erractionsave(err_act)

!!!!!!!!!!!!!!!! CHANGE FOR NON-CSR A
    !
    ! Renumbering type: 
    !     1. Global column indices
    !     (2. GPS band reduction disabled for the time being)

    if (p%iprcparm(iren_)==1) then 
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
       allocate(p%perm(nnr),p%invperm(nnr),itmp2(nnr))
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

       allocate(itmp(max(8,atmp%m+2,nztmp+2)),rtmp(atmp%m))

       j = 1
       atmp%ia2(1) = 1
       do i=1, atmp%m
          ir = p%perm(i)

          if (ir <= a%m ) then

             nzl = a%ia2(ir+1) - a%ia2(ir)
             if (nzl > size(rtmp)) then
                call psb_realloc(nzl,rtmp,info)
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
                   rtmp(k) = a%aspk(jj+kk-1)
                   atmp%ia1(j+k-1) = p%invperm(a%ia1(jj+kk-1))
                endif
             enddo
             call isrx(k,atmp%ia1(j:j+k-1),itmp2)
             do kk=1,k
                atmp%aspk(j+kk-1) = rtmp(itmp2(kk))
             enddo

          else if (ir <= atmp%m ) then 

             ir = ir - a%m
             nzl = blck%ia2(ir+1) - blck%ia2(ir)
             if (nzl > size(rtmp)) then
                call psb_realloc(nzl,rtmp,info)
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
                   rtmp(k)         = blck%aspk(jj+kk-1)
                   atmp%ia1(j+k-1) = p%invperm(blck%ia1(jj+kk-1))
                endif
             enddo
             call isrx(k,atmp%ia1(j:j+k-1),itmp2)
             do kk=1,k
                atmp%aspk(j+kk-1) = rtmp(itmp2(kk))
             enddo

          else
             write(0,*) 'Row index error 1 :',i,ir
          endif

          j = j + k
          atmp%ia2(i+1) = j

       enddo

       t4 = mpi_wtime()


       deallocate(itmp,itmp2,rtmp)

    else if (p%iprcparm(iren_)==2) then 

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

       allocate(itmp(max(8,atmp%m+2,nztmp+2)))
       itmp(1:8) = 0
!          write(0,*) me,' Renumbering: Calling Metis'
!        call blacs_barrier(icontxt,'All')

!          write(0,*) size(p%av(u_pr_)%pl),size(p%av(l_pr_)%pr)
       call  gps_reduction(atmp%m,atmp%ia2,atmp%ia1,p%perm,p%invperm,info)
       if(info/=0) then
          info=4010
          ch_err='gps_reduction'
          call psb_errpush(info,name,a_err=ch_err)
          goto 9999
       end if

!      write(0,*) me,' Renumbering: Done GPS'
       call blacs_barrier(icontxt,'All')
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

       allocate(itmp2(max(8,atmp%m+2,nztmp+2)),rtmp(atmp%m))

       j = 1
       atmp%ia2(1) = 1
       do i=1, atmp%m
          ir = p%perm(i)

          if (ir <= a%m ) then

             nzl = a%ia2(ir+1) - a%ia2(ir)
             if (nzl > size(rtmp)) then
                call psb_realloc(nzl,rtmp,info)
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
                   rtmp(k) = a%aspk(jj+kk-1)
                   atmp%ia1(j+k-1) = p%invperm(a%ia1(jj+kk-1))
                endif
             enddo
             call isrx(k,atmp%ia1(j:j+k-1),itmp2)
             do kk=1,k
                atmp%aspk(j+kk-1) = rtmp(itmp2(kk))
             enddo

          else if (ir <= atmp%m ) then 

             ir = ir - a%m
             nzl = blck%ia2(ir+1) - blck%ia2(ir)
             if (nzl > size(rtmp)) then
                call psb_realloc(nzl,rtmp,info)
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
                   rtmp(k)         = blck%aspk(jj+kk-1)
                   atmp%ia1(j+k-1) = p%invperm(blck%ia1(jj+kk-1))
                endif
             enddo
             call isrx(k,atmp%ia1(j:j+k-1),itmp2)
             do kk=1,k
                atmp%aspk(j+kk-1) = rtmp(itmp2(kk))
             enddo

          else
             write(0,*) 'Row index error 1 :',i,ir
          endif

          j = j + k
          atmp%ia2(i+1) = j

       enddo

       t4 = mpi_wtime()



       deallocate(itmp,itmp2,rtmp)

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

  end subroutine apply_renum


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


end subroutine psb_dcslu


