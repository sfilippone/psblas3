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
!*****************************************************************************
!*                                                                           *
!* This is where the action takes place.                                     *
!* ASMATBLD does the setup: building the prec descriptor plus retrieving     *
!*                           matrix rows if needed                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!* some open code does the renumbering                                       *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!*****************************************************************************
subroutine psb_dilu_bld(a,desc_a,p,upd,info)
  use psb_base_mod
  use psb_prec_type
  implicit none
  !                                                                               
  !     .. Scalar Arguments ..                                                    
  integer, intent(out)                      :: info
  !     .. array Arguments ..                                                     
  type(psb_dspmat_type), intent(in), target :: a
  type(psb_dbaseprc_type), intent(inout)    :: p
  type(psb_desc_type), intent(in)           :: desc_a
  character, intent(in)                     :: upd

  !     .. Local Scalars ..                                                       
  integer  ::    i, j, jj, k, kk, m
  integer  ::    int_err(5)
  character ::        trans, unitd
  type(psb_dspmat_type) :: blck, atmp
  real(kind(1.d0)) :: t1,t2,t3,t4,t5,t6,mpi_wtime, t7, t8
  external  mpi_wtime
  logical, parameter :: debugprt=.false., debug=.false., aggr_dump=.false.
  integer   nztota, nztotb, nztmp, nzl, nnr, ir, err_act,&
       & n_row, nrow_a,n_col, nhalo, ind, iind, i1,i2,ia
  integer :: ictxt,np,me
  character(len=20)      :: name, ch_err

  interface psb_ilu_fct
    subroutine psb_dilu_fct(a,l,u,d,info,blck)
      use psb_base_mod
      integer, intent(out)                ::     info
      type(psb_dspmat_type),intent(in)    :: a
      type(psb_dspmat_type),intent(inout) :: l,u
      type(psb_dspmat_type),intent(in), optional, target :: blck
      real(kind(1.d0)), intent(inout)     ::  d(:)
    end subroutine psb_dilu_fct
  end interface

  interface psb_asmatbld
    Subroutine psb_dasmatbld(ptype,novr,a,blk,desc_data,upd,desc_p,info,outfmt)
      use psb_base_mod
      use psb_prec_type
      integer, intent(in)                  :: ptype,novr
      Type(psb_dspmat_type), Intent(in)    ::  a
      Type(psb_dspmat_type), Intent(inout) ::  blk
      Type(psb_desc_type), Intent(inout)   :: desc_p
      Type(psb_desc_type), Intent(in)      :: desc_data 
      Character, Intent(in)                :: upd
      integer, intent(out)                 :: info
      character(len=5), optional           :: outfmt
    end Subroutine psb_dasmatbld
  end interface

  interface psb_sp_renum
    subroutine psb_dsp_renum(a,desc_a,blck,p,atmp,info)
      use psb_base_mod
      use psb_prec_type
      implicit none
      type(psb_dspmat_type), intent(in)      :: a,blck
      type(psb_dspmat_type), intent(inout)   :: atmp
      type(psb_dbaseprc_type), intent(inout) :: p
      type(psb_desc_type), intent(in)        :: desc_a
      integer, intent(out)   :: info
    end subroutine psb_dsp_renum
  end interface

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='psb_ilu_bld'
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)
  call psb_info(ictxt, me, np)

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

  call psb_nullify_sp(blck)
  call psb_nullify_sp(atmp)

  t1= mpi_wtime()

  if(debug) write(0,*)me,': calling psb_asmatbld',p%iprcparm(p_type_),p%iprcparm(n_ovr_)
  if (debug) call psb_barrier(ictxt)
  call psb_asmatbld(p%iprcparm(p_type_),p%iprcparm(n_ovr_),a,&
       & blck,desc_a,upd,p%desc_data,info)
  if(info/=0) then
    info=4010
    ch_err='psb_asmatbld'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if
  t2= mpi_wtime()
  if (debug) write(0,*)me,': out of psb_asmatbld'
  if (debug) call psb_barrier(ictxt)

  if (allocated(p%av)) then 
    if (size(p%av) < bp_ilu_avsz) then 
      call psb_errpush(4010,name,a_err='Insufficient av size')
      goto 9999      
    endif
  else
    call psb_errpush(4010,name,a_err='AV not associated')
    goto 9999      
  endif
!!$  call psb_csprt(50+me,a,head='% (A)')    

  nrow_a = psb_cd_get_local_rows(desc_a)
  nztota = psb_sp_get_nnzeros(a)
  nztotb = psb_sp_get_nnzeros(blck)
  if (debug) write(0,*)me,': out get_nnzeros',nztota
  if (debug) call psb_barrier(ictxt)

  n_col  = psb_cd_get_local_cols(desc_a)
  nhalo  = n_col-nrow_a
  n_row  = p%desc_data%matrix_data(psb_n_row_)
  p%av(l_pr_)%m  = n_row
  p%av(l_pr_)%k  = n_row
  p%av(u_pr_)%m  = n_row
  p%av(u_pr_)%k  = n_row
  call psb_sp_all(n_row,n_row,p%av(l_pr_),nztota+nztotb,info)
  if (info == 0) call psb_sp_all(n_row,n_row,p%av(u_pr_),nztota+nztotb,info)
  if(info/=0) then
    info=4010
    ch_err='psb_sp_all'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (allocated(p%d)) then 
    if (size(p%d) < n_row) then 
      deallocate(p%d)
    endif
  endif
  if (.not.allocated(p%d)) then 
    allocate(p%d(n_row),stat=info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

  endif


  if (debug) then 
    write(0,*) me,'Done psb_asmatbld'
    call psb_barrier(ictxt)
  endif


  if (p%iprcparm(iren_) > 0) then 

    !
    ! Here we allocate a full copy to hold local A and received BLK
    !

    nztota = psb_sp_get_nnzeros(a)
    nztotb = psb_sp_get_nnzeros(blck)
    call psb_sp_all(atmp,nztota+nztotb,info)
    if(info/=0) then
      info=4011
      call psb_errpush(info,name)
      goto 9999
    end if


    !      write(0,*) 'ILU_BLD ',nztota,nztotb,a%m

    call  psb_sp_renum(a,desc_a,blck,p,atmp,info)

    if(info/=0) then
      info=4010
      ch_err='psb_sp_renum'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    t3 = mpi_wtime()
    if (debugprt) then 
      call psb_barrier(ictxt)
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
    call psb_ilu_fct(atmp,p%av(l_pr_),p%av(u_pr_),p%d,info,blck=blck)
    if(info/=0) then
      info=4010
      ch_err='psb_ilu_fct'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    call psb_sp_free(atmp,info) 
    if(info/=0) then
      info=4010
      ch_err='psb_sp_free'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if


  else if (p%iprcparm(iren_) == 0) then
    t3 = mpi_wtime()
    ! This is where we have mo renumbering, thus no need 
    ! for ATMP

    if (debugprt) then 
      open(40+me)
      call psb_barrier(ictxt)
      call psb_csprt(40+me,a,iv=p%desc_data%loc_to_glob,&
           &    head='% Local matrix')
      if (p%iprcparm(p_type_)==asm_) then 
        call psb_csprt(40+me,blck,iv=p%desc_data%loc_to_glob,&
             &  irs=a%m,head='% Received rows')
      endif
      close(40+me)
    endif

    t5= mpi_wtime()
    if (debug) write(0,*) me,' Going for ilu_fct'
    if (debug) call psb_barrier(ictxt)
    call psb_ilu_fct(a,p%av(l_pr_),p%av(u_pr_),p%d,info,blck=blck)
    if(info/=0) then
      info=4010
      ch_err='psb_ilu_fct'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if
    if (debug) write(0,*) me,' Done dilu_fct'
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

!!$  call psb_csprt(60+me,a,head='% (A)')    


  !    ierr = MPE_Log_event( ifcte, 0, "st SIMPLE" )
  t6 = mpi_wtime()
  !
  !    write(0,'(i3,1x,a,3(1x,g18.9))') me,'renum/factor time',t3-t2,t6-t5
  !    if (me==0) write(0,'(a,3(1x,g18.9))') 'renum/factor time',t3-t2,t6-t5

  call psb_sp_free(blck,info)
  if(info/=0) then
    info=4010
    ch_err='psb_sp_free'
    call psb_errpush(info,name,a_err=ch_err)
    goto 9999
  end if

  if (psb_sp_getifld(psb_upd_,p%av(u_pr_),info) /= psb_upd_perm_) then
    call psb_sp_trimsize(p%av(u_pr_),i1,i2,ia,info)
    if (info == 0) call psb_sp_reall(p%av(u_pr_),i1,i2,ia,info)
  endif

  if (psb_sp_getifld(psb_upd_,p%av(l_pr_),info) /= psb_upd_perm_) then
    call psb_sp_trimsize(p%av(l_pr_),i1,i2,ia,info)
    if (info == 0) call psb_sp_reall(p%av(l_pr_),i1,i2,ia,info)
  endif


  if (debug) write(0,*) me,'End of ilu_bld'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act.eq.act_abort) then
    call psb_error()
    return
  end if
  return


end subroutine psb_dilu_bld


