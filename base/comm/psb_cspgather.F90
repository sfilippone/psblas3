subroutine  psb_csp_allgather(globa, loca, desc_a, info, root, dupl,keepnum,keeploc)
  use psb_descriptor_type
  use psb_error_mod
  use psb_mat_mod
  use psb_penv_mod
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif
  type(psb_cspmat_type), intent(inout) :: loca
  type(psb_cspmat_type), intent(inout) :: globa
  type(psb_desc_type), intent(in) :: desc_a
  integer, intent(out)            :: info
  integer, intent(in), optional   :: root, dupl
  logical, intent(in), optional   :: keepnum,keeploc

  type(psb_c_coo_sparse_mat)      :: loc_coo, glob_coo
  integer :: ictxt,np,me, err_act, icomm, dupl_, nrg, ncg, nzg
  integer :: ip, ndx,naggrm1,naggrp1, i, j, k
  logical :: keepnum_, keeploc_
  integer, allocatable :: nzbr(:), idisp(:)
  character(len=20) :: name
  integer            :: debug_level, debug_unit

  name='psb_gather'
  if (psb_get_errstatus().ne.0) return 
  info=psb_success_

  call psb_erractionsave(err_act)
  ictxt = psb_cd_get_context(desc_a)
  icomm = psb_cd_get_mpic(desc_a)
  call psb_info(ictxt, me, np)

  if (present(keepnum)) then 
    keepnum_ = keepnum
  else
    keepnum_ = .true.
  end if
  if (present(keeploc)) then 
    keeploc_ = keeploc
  else
    keeploc_ = .true.
  end if
  call globa%free()

  if (keepnum_) then 
    nrg = psb_cd_get_global_rows(desc_a)
    ncg = psb_cd_get_global_rows(desc_a)

    allocate(nzbr(np), idisp(np),stat=info)
    if (info /= psb_success_) then 
      info=psb_err_alloc_request_
      call psb_errpush(info,name,i_err=(/2*np,0,0,0,0/),&
           & a_err='integer')
      goto 9999      
    end if
    call loca%mv_to(loc_coo)
    nzbr(:) = 0
    nzbr(me+1) = loc_coo%get_nzeros()
    call psb_sum(ictxt,nzbr(1:np))
    nzg = sum(nzbr)
    if (info == psb_success_) call glob_coo%allocate(nrg,ncg,nzg)
    if (info /= psb_success_) goto 9999
    do ip=1,np
      idisp(ip) = sum(nzbr(1:ip-1))
    enddo
    ndx = nzbr(me+1) 
    call mpi_allgatherv(loc_coo%val,ndx,mpi_complex,&
         & glob_coo%val,nzbr,idisp,&
         & mpi_double_precision,icomm,info)
    if (info == psb_success_) call mpi_allgatherv(loc_coo%ia,ndx,mpi_integer,&
         & glob_coo%ia,nzbr,idisp,&
         & mpi_integer,icomm,info)
    if (info == psb_success_) call mpi_allgatherv(loc_coo%ja,ndx,mpi_integer,&
         & glob_coo%ja,nzbr,idisp,&
         & mpi_integer,icomm,info)
    
    if (info /= psb_success_) then 
      call psb_errpush(psb_err_internal_error_,name,a_err=' from mpi_allgatherv')
      goto 9999
    end if
    
    if (keeploc_) then
      call loca%mv_from(loc_coo)
    else
      call loc_coo%free()
    end if
    call glob_coo%set_nzeros(nzg)
    if (present(dupl)) call glob_coo%set_dupl(dupl)
    call globa%mv_from(glob_coo)

  else
    write(psb_err_unit,*) 'SP_ALLGATHER: Not implemented yet with keepnum ',keepnum_
  end if



  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_errpush(info,name)
  call psb_erractionrestore(err_act)
  if (err_act.eq.psb_act_abort_) then
    call psb_error()
    return
  end if
  return

 
end subroutine psb_csp_allgather