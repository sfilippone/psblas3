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
subroutine psb_zbldaggrmat(a,desc_a,ac,desc_ac,p,info)
  use psb_base_mod
  use psb_prec_type
  implicit none

  type(psb_zspmat_type), intent(in), target  :: a
  type(psb_zbaseprc_type), intent(inout),target     :: p
  type(psb_zspmat_type), intent(out), target :: ac
  type(psb_desc_type), intent(in)            :: desc_a
  type(psb_desc_type), intent(inout)         :: desc_ac
  integer, intent(out)                       :: info

  logical, parameter :: aggr_dump=.false.
  integer ::ictxt,np,me, err_act
  character(len=20) :: name, ch_err
  name='psb_zbldaggrmat'
  if(psb_get_errstatus().ne.0) return 
  info=0
  call psb_erractionsave(err_act)

  ictxt=psb_cd_get_context(desc_a)

  call psb_info(ictxt, me, np)

  select case (p%iprcparm(smth_kind_))
  case (no_smth_) 

    call raw_aggregate(info)

    if(info /= 0) then
      call psb_errpush(4010,name,a_err='raw_aggregate')
      goto 9999
    end if
    if (aggr_dump) call psb_csprt(90+me,ac,head='% Raw aggregate.')

  case(smth_omg_,smth_biz_) 

    call smooth_aggregate(info)

    if(info /= 0) then
      call psb_errpush(4010,name,a_err='smooth_aggregate')
      goto 9999
    end if
  case default
    call psb_errpush(4010,name,a_err=name)
    goto 9999

  end select

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

  subroutine raw_aggregate(info)
    use psb_base_mod
    use psb_prec_type
    implicit none

    include 'mpif.h'
    integer, intent(out)   :: info
    type(psb_zspmat_type)          :: b, tmp
    integer, pointer :: nzbr(:), idisp(:)
    integer :: ictxt, nrow, nglob, ncol, ntaggr, nzac, ip, ndx,&
         & naggr, np, me, nzt,jl,nzl,nlr,&
         & icomm,naggrm1, i, j, k, err_act

    name='raw_aggregate'
    if(psb_get_errstatus().ne.0) return 
    info=0
    call psb_erractionsave(err_act)

    call psb_nullify_sp(b)

    ictxt = psb_cd_get_context(desc_a)
    call psb_info(ictxt, me, np)
    nglob = psb_cd_get_global_rows(desc_a)
    nrow  = psb_cd_get_local_rows(desc_a)
    ncol  = psb_cd_get_local_cols(desc_a)

    naggr  = p%nlaggr(me+1)
    ntaggr = sum(p%nlaggr)
    allocate(nzbr(np), idisp(np),stat=info)

    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    naggrm1=sum(p%nlaggr(1:me))

    if (p%iprcparm(coarse_mat_) == mat_repl_) then
      do i=1, nrow
        p%mlia(i) = p%mlia(i) + naggrm1
      end do
    end if
    call psb_halo(p%mlia,desc_a,info)

    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_halo')
      goto 9999
    end if

    nzt = psb_sp_get_nnzeros(a)

    call psb_sp_all(b,nzt,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='spall')
      goto 9999
    end if

    call psb_sp_setifld(psb_dupl_ovwrt_,psb_dupl_,b,info)
    call psb_sp_setifld(psb_upd_dflt_,psb_upd_,b,info)
    b%fida = 'COO'
    b%m=a%m
    b%k=a%k
    call psb_csdp(a,b,info)
    if(info /= 0) then
      info=4010
      ch_err='psb_csdp'
      call psb_errpush(info,name,a_err=ch_err)
      goto 9999
    end if

    nzt = psb_sp_get_nnzeros(b)

    do i=1, nzt 
      b%ia1(i) = p%mlia(b%ia1(i))
      b%ia2(i) = p%mlia(b%ia2(i))
    enddo
    call psb_fixcoo(b,info)

    nzt = psb_sp_get_nnzeros(b)

    call psb_sp_reall(b,nzt,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='spreall')
      goto 9999
    end if
    b%m = naggr
    b%k = naggr

    if (p%iprcparm(coarse_mat_) == mat_repl_) then 

      call psb_cdrep(ntaggr,ictxt,desc_ac,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_cdrep')
        goto 9999
      end if

      nzbr(:) = 0
      nzbr(me+1) = nzt
      call psb_sum(ictxt,nzbr(1:np))
      nzac = sum(nzbr)
      call psb_sp_all(ntaggr,ntaggr,ac,nzac,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='spall')
        goto 9999
      end if

      call psb_get_mpicomm(ictxt,icomm)
      do ip=1,np
        idisp(ip) = sum(nzbr(1:ip-1))
      enddo
      ndx = nzbr(me+1) 

      call mpi_allgatherv(b%aspk,ndx,mpi_double_complex,ac%aspk,nzbr,idisp,&
           & mpi_double_complex,icomm,info)
      call mpi_allgatherv(b%ia1,ndx,mpi_integer,ac%ia1,nzbr,idisp,&
           & mpi_integer,icomm,info)
      call mpi_allgatherv(b%ia2,ndx,mpi_integer,ac%ia2,nzbr,idisp,&
           & mpi_integer,icomm,info)
      if(info /= 0) then
        info=-1
        call psb_errpush(info,name)
        goto 9999
      end if

      ac%m = ntaggr
      ac%k = ntaggr
      ac%infoa(psb_nnz_) = nzac
      ac%fida='COO'
      ac%descra='G'
      call psb_fixcoo(ac,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='sp_free')
        goto 9999
      end if

    else if (p%iprcparm(coarse_mat_) == mat_distr_) then 

      call psb_cdall(ictxt,desc_ac,info,nl=naggr)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_cdall')
        goto 9999
      end if
      call psb_cdasb(desc_ac,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_cdasb')
        goto 9999
      end if

      call psb_sp_clone(b,ac,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='spclone')
        goto 9999
      end if
      call psb_sp_free(b,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='sp_free')
        goto 9999
      end if

      !if(.not.associated(p%av(ap_nd_)%aspk)) p%iprcparm(jac_sweeps_) = 1
      !------------------------------------------------------------------
      ! Split AC=M+N  N off-diagonal part
      call psb_sp_all(ac%m,ac%k,p%av(ap_nd_),nzl,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_sp_all')
        goto 9999
      end if
      if(.not.allocated(p%av(ap_nd_)%aspk)) write(0,*) '.not.associated(p%av(ap_nd_)%ia1)'
      if(.not.allocated(p%av(ap_nd_)%ia1)) write(0,*) '.not.associated(p%av(ap_nd_)%ia1)'
      !write(0,*) 'ok line 238'

      k=0
      do i=1,nzl
        if (ac%ia2(i)>ac%m) then 
          k = k + 1
          p%av(ap_nd_)%aspk(k) = ac%aspk(i)
          p%av(ap_nd_)%ia1(k)  = ac%ia1(i)
          p%av(ap_nd_)%ia2(k)  = ac%ia2(i)
        endif
      enddo
      p%av(ap_nd_)%infoa(psb_nnz_) = k

      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_ipcoo2csr')
        goto 9999
      end if
      call psb_sum(ictxt,k)

      if (k == 0) then 
        ! If the off diagonal part is emtpy, there's no point 
        ! in doing multiple  Jacobi sweeps. This is certain 
        ! to happen when running on a single processor.
        p%iprcparm(jac_sweeps_) = 1
      end if
      !write(0,*) 'operations in bldaggrmat are ok !'
      !------------------------------------------------------------------

      call psb_ipcoo2csr(p%av(ap_nd_),info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='ipcoo2csr')
        goto 9999
      end if

    else

      write(0,*) 'Unknown p%iprcparm(coarse_mat) in aggregate_sp',p%iprcparm(coarse_mat_)
    end if

    call psb_ipcoo2csr(ac,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='ipcoo2csr')
      goto 9999
    end if

    deallocate(nzbr,idisp)

    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_erractionrestore(err_act)
    if (err_act.eq.act_abort) then
      call psb_error()
      return
    end if
    return

  end subroutine raw_aggregate



  subroutine smooth_aggregate(info)
    use psb_base_mod
    use psb_prec_type
    use mpi
    implicit none 

    integer, intent(out)   :: info

    type(psb_zspmat_type)          :: b
    integer, pointer :: nzbr(:), idisp(:), ivall(:)
    integer :: ictxt, nrow, nglob, ncol, ntaggr, nzac, ip, ndx,&
         & naggr, np, me, &
         & icomm, naggrm1,naggrp1,i,j,err_act,k,nzl
    type(psb_zspmat_type), pointer  :: am1,am2
    type(psb_zspmat_type) :: am3,am4
    logical       :: ml_global_nmb

    logical, parameter :: test_dump=.false., debug=.false.
    integer, parameter :: ncmax=16
    real(kind(1.d0))   :: omega, anorm, tmp, dg
    character(len=20) :: name


    name='smooth_aggregate'
    if(psb_get_errstatus().ne.0) return 
    info=0
    call psb_erractionsave(err_act)

    ictxt = psb_cd_get_context(desc_a)

    call psb_info(ictxt, me, np)

    call psb_nullify_sp(b)
    call psb_nullify_sp(am3)
    call psb_nullify_sp(am4)

    am2 => p%av(sm_pr_t_)
    am1 => p%av(sm_pr_)

    nglob = psb_cd_get_global_rows(desc_a)
    nrow  = psb_cd_get_local_rows(desc_a)
    ncol  = psb_cd_get_local_cols(desc_a)

    naggr  = p%nlaggr(me+1)
    ntaggr = sum(p%nlaggr)

    allocate(nzbr(np), idisp(np),stat=info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if


    naggrm1 = sum(p%nlaggr(1:me))
    naggrp1 = sum(p%nlaggr(1:me+1))

    ml_global_nmb = ( (p%iprcparm(smth_kind_) == smth_omg_).or.&
         & ( (p%iprcparm(smth_kind_) == smth_biz_).and.&
         &    (p%iprcparm(coarse_mat_) == mat_repl_)) ) 


    if (ml_global_nmb) then 
      p%mlia(1:nrow) = p%mlia(1:nrow) + naggrm1
      call psb_halo(p%mlia,desc_a,info)

      if(info /= 0) then
        call psb_errpush(4010,name,a_err='f90_pshalo')
        goto 9999
      end if
    end if

    if (aggr_dump) then 
      open(30+me)
      write(30+me,*) '% Aggregation map'
      do i=1,ncol
        write(30+me,*) i,p%mlia(i)
      end do
      close(30+me)
    end if

    ! naggr: number of local aggregates
    ! nrow: local rows. 
    ! 
    allocate(p%dorig(nrow),stat=info)
    if (info /= 0) then 
      call psb_errpush(4010,name,a_err='Allocate')
      goto 9999      
    end if

    ! Get diagonal D
    call psb_sp_getdiag(a,p%dorig,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='sp_getdiag')
      goto 9999
    end if

    do i=1,size(p%dorig)
      if (p%dorig(i) /= zzero) then
        p%dorig(i) = zone / p%dorig(i)
      else
        p%dorig(i) = zone
      end if
    end do

    !     where (p%dorig /= dzero) 
    !       p%dorig = done / p%dorig
    !     elsewhere
    !       p%dorig = done
    !     end where


    ! 1. Allocate Ptilde in sparse matrix form 
    am4%fida='COO'
    am4%m=ncol
    if (ml_global_nmb) then 
      am4%k=ntaggr
      call psb_sp_all(ncol,ntaggr,am4,ncol,info)
    else 
      am4%k=naggr
      call psb_sp_all(ncol,naggr,am4,ncol,info)
    endif
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='spall')
      goto 9999
    end if

    if (ml_global_nmb) then 
      do i=1,ncol
        am4%aspk(i) = zone
        am4%ia1(i)  = i
        am4%ia2(i)  = p%mlia(i)  
      end do
      am4%infoa(psb_nnz_) = ncol
    else
      do i=1,nrow
        am4%aspk(i) = zone
        am4%ia1(i)  = i
        am4%ia2(i)  = p%mlia(i)  
      end do
      am4%infoa(psb_nnz_) = nrow
    endif




    call psb_ipcoo2csr(am4,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='ipcoo2csr')
      goto 9999
    end if

    call psb_sp_clone(a,am3,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='spclone')
      goto 9999
    end if

    !
    ! WARNING: the cycles below assume that AM3 does have 
    ! its diagonal elements stored explicitly!!! 
    ! Should we switch to something safer? 
    !
    call psb_sp_scal(am3,p%dorig,info)
    if(info /= 0) goto 9999

    if (p%iprcparm(om_choice_) == lib_choice_) then 

      if (p%iprcparm(smth_kind_) == smth_biz_) then 

        ! 
        ! This only works with CSR.
        !
        anorm = dzero
        dg    = done
        do i=1,am3%m
          tmp = dzero
          do j=am3%ia2(i),am3%ia2(i+1)-1
            if (am3%ia1(j) <= am3%m) then 
              tmp = tmp + abs(am3%aspk(j))
            endif
            if (am3%ia1(j) == i ) then 
              dg = abs(am3%aspk(j))
            end if
          end do
          anorm = max(anorm,tmp/dg) 
        enddo

        call psb_amx(ictxt,anorm)     
      else
        anorm = psb_spnrmi(am3,desc_a,info)
      endif
      omega = 4.d0/(3.d0*anorm)
      p%dprcparm(smooth_omega_) = omega 

    else if (p%iprcparm(om_choice_) == user_choice_) then 

      omega = p%dprcparm(smooth_omega_) 

    else if (p%iprcparm(om_choice_) /= user_choice_) then 
      write(0,*) me,'Error: invalid choice for OMEGA in blaggrmat?? ',&
           &   p%iprcparm(om_choice_)    
    end if


    if (am3%fida=='CSR') then 
      do i=1,am3%m
        do j=am3%ia2(i),am3%ia2(i+1)-1
          if (am3%ia1(j) == i) then 
            am3%aspk(j) = done - omega*am3%aspk(j) 
          else
            am3%aspk(j) = - omega*am3%aspk(j) 
          end if
        end do
      end do
    else  if (am3%fida=='COO') then 
      do j=1,am3%infoa(psb_nnz_) 
        if (am3%ia1(j) /= am3%ia2(j)) then 
          am3%aspk(j) = - omega*am3%aspk(j) 
        else
          am3%aspk(j) = done - omega*am3%aspk(j) 
        endif
      end do
      call psb_ipcoo2csr(am3,info)      
    else
      write(0,*) 'Missing implementation of I sum' 
      call psb_errpush(4010,name)
      goto 9999
    end if

    if (test_dump) then 
      open(30+me)
      write(30+me,*) 'OMEGA: ',omega
      do i=1,size(p%dorig)
        write(30+me,*) p%dorig(i) 
      end do
      close(30+me)
    end if

    if (test_dump)  call &
         & psb_csprt(20+me,am4,head='% Operator Ptilde.',ivr=desc_a%loc_to_glob)
    if (test_dump) call psb_csprt(40+me,am3,head='% (I-wDA)',ivr=desc_a%loc_to_glob,&
         & ivc=desc_a%loc_to_glob)    
    if (debug) write(0,*) me,'Done gather, going for SYMBMM 1'
    !
    ! Symbmm90 does the allocation for its result.
    ! 
    ! am1 = (i-wDA)Ptilde
    ! Doing it this way means to consider diag(Ai)
    ! 
    !
    call psb_symbmm(am3,am4,am1,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='symbmm 1')
      goto 9999
    end if

    call psb_numbmm(am3,am4,am1)

    if (debug) write(0,*) me,'Done NUMBMM 1'

    call psb_sp_free(am4,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='sp_free')
      goto 9999
    end if

    if (ml_global_nmb) then 
      !
      ! Now we have to gather the halo of am1, and add it to itself
      ! to multiply it by A,
      !
      call psb_sphalo(am1,desc_a,am4,info,clcnv=.false.)

      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_sphalo')
        goto 9999
      end if

      call psb_rwextd(ncol,am1,info,b=am4)      
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_rwextd')
        goto 9999
      end if

      call psb_sp_free(am4,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_sp_free')
        goto 9999
      end if

    else 

      call psb_rwextd(ncol,am1,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='rwextd')
        goto 9999
      end if
    endif

    if (test_dump) &
         & call psb_csprt(60+me,am1,head='% (I-wDA)Pt',ivr=desc_a%loc_to_glob)    

    call psb_symbmm(a,am1,am3,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='symbmm 2')
      goto 9999
    end if

    call psb_numbmm(a,am1,am3)
    if (debug) write(0,*) me,'Done NUMBMM 2'

    if  (p%iprcparm(smth_kind_) == smth_omg_) then 
      call psb_transc(am1,am2,fmt='COO')
      nzl = am2%infoa(psb_nnz_)
      i=0
      !
      ! Now we have to fix this.  The only rows of B that are correct 
      ! are those corresponding to "local" aggregates, i.e. indices in p%mlia(:)
      !
      do k=1, nzl
        if ((naggrm1 < am2%ia1(k)) .and.(am2%ia1(k) <= naggrp1)) then
          i = i+1
          am2%aspk(i) = am2%aspk(k)
          am2%ia1(i)  = am2%ia1(k)
          am2%ia2(i)  = am2%ia2(k)
        end if
      end do

      am2%infoa(psb_nnz_) = i
      call psb_ipcoo2csr(am2,info)
    else
      call psb_transc(am1,am2)
    endif
    if (debug) write(0,*) me,'starting sphalo/ rwxtd'

    if  (p%iprcparm(smth_kind_) == smth_omg_) then 
      ! am2 = ((i-wDA)Ptilde)^T
      call psb_sphalo(am3,desc_a,am4,info,clcnv=.false.)

      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_sphalo')
        goto 9999
      end if
      call psb_rwextd(ncol,am3,info,b=am4)      
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_rwextd')
        goto 9999
      end if
      call psb_sp_free(am4,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_sp_free')
        goto 9999
      end if

    else if  (p%iprcparm(smth_kind_) == smth_biz_) then 

      call psb_rwextd(ncol,am3,info)
      if(info /= 0) then
        call psb_errpush(4010,name,a_err='psb_rwextd')
        goto 9999
      end if
    endif

    if (debug) write(0,*) me,'starting symbmm 3'
    call psb_symbmm(am2,am3,b,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='symbmm 3')
      goto 9999
    end if

    if (debug) write(0,*) me,'starting numbmm 3'
    call psb_numbmm(am2,am3,b)
    if (debug) write(0,*) me,'Done NUMBMM 3'

!!$    if (aggr_dump) call csprt(50+me,am1,head='% Operator PTrans.')
    call psb_sp_free(am3,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_sp_free')
      goto 9999
    end if

    call psb_ipcsr2coo(b,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='ipcsr2coo')
      goto 9999
    end if

    call psb_fixcoo(b,info)      
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='fixcoo')
      goto 9999
    end if


    if (test_dump) call psb_csprt(80+me,b,head='% Smoothed aggregate AC.')    

    select case(p%iprcparm(smth_kind_))

    case(smth_omg_) 

      select case(p%iprcparm(coarse_mat_))

      case(mat_distr_) 

        call psb_sp_clone(b,ac,info)
        if(info /= 0) goto 9999
        nzac = ac%infoa(psb_nnz_) 
        nzl =  ac%infoa(psb_nnz_) 

        allocate(ivall(ntaggr),stat=info)
        if (info /= 0) then 
          call psb_errpush(4010,name,a_err='Allocate')
          goto 9999      
        end if

        i = 1
        do ip=1,np
          do k=1, p%nlaggr(ip)
            ivall(i) = ip
            i = i + 1
          end do
        end do

        call psb_cdall(ictxt,desc_ac,info,vg=ivall(1:ntaggr),flag=1)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_cdall')
          goto 9999
        end if


        call psb_cdins(nzl,ac%ia1,ac%ia2,desc_ac,info)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_cdins')
          goto 9999
        end if

        if (debug) write(0,*) me,'Created aux descr. distr.'
        call psb_cdasb(desc_ac,info)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_cdasb')
          goto 9999
        end if


        if (debug) write(0,*) me,'Asmbld aux descr. distr.'

        call psb_glob_to_loc(ac%ia1(1:nzl),desc_ac,info,iact='I')
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psglob_to_loc')
          goto 9999
        end if


        call psb_glob_to_loc(ac%ia2(1:nzl),desc_ac,info,iact='I')
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psglob_to_loc')
          goto 9999
        end if


        ac%m=desc_ac%matrix_data(psb_n_row_)
        ac%k=desc_ac%matrix_data(psb_n_col_)
        ac%fida='COO'
        ac%descra='G'

        call psb_sp_free(b,info)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_sp_free')
          goto 9999
        end if


        deallocate(ivall,nzbr,idisp)

        ! Split AC=M+N  N off-diagonal part
        call psb_sp_all(ac%m,ac%k,p%av(ap_nd_),nzl,info)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_sp_all')
          goto 9999
        end if

        k=0
        do i=1,nzl
          if (ac%ia2(i)>ac%m) then 
            k = k + 1
            p%av(ap_nd_)%aspk(k) = ac%aspk(i)
            p%av(ap_nd_)%ia1(k)  = ac%ia1(i)
            p%av(ap_nd_)%ia2(k)  = ac%ia2(i)
          endif
        enddo
        p%av(ap_nd_)%infoa(psb_nnz_) = k
        call psb_ipcoo2csr(p%av(ap_nd_),info)

        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_ipcoo2csr')
          goto 9999
        end if
        call psb_sum(ictxt,k)

        if (k == 0) then 
          ! If the off diagonal part is emtpy, there's no point 
          ! in doing multiple  Jacobi sweeps. This is certain 
          ! to happen when running on a single processor.
          p%iprcparm(jac_sweeps_) = 1
        end if


        if (np>1) then 
          nzl = psb_sp_get_nnzeros(am1)
          call psb_glob_to_loc(am1%ia1(1:nzl),desc_ac,info,'I')

          if(info /= 0) then
            call psb_errpush(4010,name,a_err='psb_glob_to_loc')
            goto 9999
          end if
        endif
        am1%k=desc_ac%matrix_data(psb_n_col_)

        if (np>1) then 
          call psb_ipcsr2coo(am2,info)
          if(info /= 0) then
            call psb_errpush(4010,name,a_err='psb_ipcsr2coo')
            goto 9999
          end if

          nzl = am2%infoa(psb_nnz_) 
          call psb_glob_to_loc(am2%ia1(1:nzl),desc_ac,info,'I')
          if(info /= 0) then
            call psb_errpush(4010,name,a_err='psb_glob_to_loc')
            goto 9999
          end if

          call psb_ipcoo2csr(am2,info)
          if(info /= 0) then
            call psb_errpush(4010,name,a_err='psb_ipcoo2csr')
            goto 9999
          end if
        end if
        am2%m=desc_ac%matrix_data(psb_n_col_)

      case(mat_repl_) 
        !
        !
        call psb_cdrep(ntaggr,ictxt,desc_ac,info)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_cdrep')
          goto 9999
        end if

        nzbr(:) = 0
        nzbr(me+1) = b%infoa(psb_nnz_)

        call psb_sum(ictxt,nzbr(1:np))
        nzac = sum(nzbr)
        call psb_sp_all(ntaggr,ntaggr,ac,nzac,info)
        if(info /= 0) goto 9999


        call psb_get_mpicomm(ictxt,icomm)
        do ip=1,np
          idisp(ip) = sum(nzbr(1:ip-1))
        enddo
        ndx = nzbr(me+1) 

        call mpi_allgatherv(b%aspk,ndx,mpi_double_complex,ac%aspk,nzbr,idisp,&
             & mpi_double_complex,icomm,info)
        call mpi_allgatherv(b%ia1,ndx,mpi_integer,ac%ia1,nzbr,idisp,&
             & mpi_integer,icomm,info)
        call mpi_allgatherv(b%ia2,ndx,mpi_integer,ac%ia2,nzbr,idisp,&
             & mpi_integer,icomm,info)
        if(info /= 0) goto 9999


        ac%m = ntaggr
        ac%k = ntaggr
        ac%infoa(psb_nnz_) = nzac
        ac%fida='COO'
        ac%descra='G'
        call psb_fixcoo(ac,info)
        if(info /= 0) goto 9999
        call psb_sp_free(b,info)
        if(info /= 0) goto 9999
        if (me==0) then 
          if (test_dump) call psb_csprt(80+me,ac,head='% Smoothed aggregate AC.')    
        endif

        deallocate(nzbr,idisp)

      case default 
        write(0,*) 'Inconsistent input in smooth_new_aggregate'
      end select


    case(smth_biz_) 

      select case(p%iprcparm(coarse_mat_))

      case(mat_distr_) 

        call psb_sp_clone(b,ac,info)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='spclone')
          goto 9999
        end if
        call psb_cdall(ictxt,desc_ac,info,nl=naggr)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_cdall')
          goto 9999
        end if
        call psb_cdasb(desc_ac,info)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_cdasb')
          goto 9999
        end if
        call psb_sp_free(b,info)
        if(info /=  0) then
          call psb_errpush(4010,name,a_err='sp_free')
          goto 9999
        end if


      case(mat_repl_) 
        !
        !

        call psb_cdrep(ntaggr,ictxt,desc_ac,info)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_cdrep')
          goto 9999
        end if

        nzbr(:) = 0
        nzbr(me+1) = b%infoa(psb_nnz_)
        call psb_sum(ictxt,nzbr(1:np))
        nzac = sum(nzbr)
        call psb_sp_all(ntaggr,ntaggr,ac,nzac,info)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_sp_all')
          goto 9999
        end if

        call psb_get_mpicomm(ictxt,icomm)
        do ip=1,np
          idisp(ip) = sum(nzbr(1:ip-1))
        enddo
        ndx = nzbr(me+1) 

        call mpi_allgatherv(b%aspk,ndx,mpi_double_complex,ac%aspk,nzbr,idisp,&
             & mpi_double_complex,icomm,info)
        call mpi_allgatherv(b%ia1,ndx,mpi_integer,ac%ia1,nzbr,idisp,&
             & mpi_integer,icomm,info)
        call mpi_allgatherv(b%ia2,ndx,mpi_integer,ac%ia2,nzbr,idisp,&
             & mpi_integer,icomm,info)
        if(info /= 0) then
          info=-1
          call psb_errpush(info,name)
          goto 9999
        end if


        ac%m = ntaggr
        ac%k = ntaggr
        ac%infoa(psb_nnz_) = nzac
        ac%fida='COO'
        ac%descra='G'
        call psb_fixcoo(ac,info)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_fixcoo')
          goto 9999
        end if
        call psb_sp_free(b,info)
        if(info /= 0) then
          call psb_errpush(4010,name,a_err='psb_sp_free')
          goto 9999
        end if

      end select
      deallocate(nzbr,idisp)

    end select

    call psb_ipcoo2csr(ac,info)
    if(info /= 0) then
      call psb_errpush(4010,name,a_err='psb_ipcoo2csr')
      goto 9999
    end if

    if (debug) write(0,*) me,'Done smooth_aggregate '
    call psb_erractionrestore(err_act)
    return

9999 continue
    call psb_errpush(info,name)
    call psb_erractionrestore(err_act)
    if (err_act.eq.act_abort) then
      call psb_error()
      return
    end if
    return


  end subroutine smooth_aggregate



end subroutine psb_zbldaggrmat
