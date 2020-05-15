
subroutine psi_build_mtpart(n,ja,irp,nparts,graph_vect,weights)
  use psb_base_mod
  use iso_c_binding
  implicit none 
  integer(psb_ipk_), intent(in) :: n, nparts
  integer(psb_ipk_), intent(in) :: ja(:), irp(:)
  integer(psb_ipk_), allocatable, intent(inout) :: graph_vect(:)
  real(psb_spk_),optional, intent(in) :: weights(:)
  ! local variables
  integer(psb_ipk_) :: i,numflag, nedc,wgflag
  integer(psb_ipk_) :: iopt(10),idummy(2),jdummy(2), info
  integer(psb_ipk_) :: nl,nptl
  integer(psb_ipk_), allocatable :: irpl(:),jal(:),gvl(:)
  real(psb_spk_),allocatable  :: wgh_(:)

#if defined(HAVE_METIS) && defined(IPK4)
  interface 
    ! subroutine METIS_PartGraphKway(n,ixadj,iadj,ivwg,iajw,&
    !     & wgflag,numflag,nparts,weights,iopt,nedc,part) bind(c)
    !   use iso_c_binding
    !   integer(c_int) :: n,wgflag,numflag,nparts,nedc
    !   integer(c_int) :: ixadj(*),iadj(*),ivwg(*),iajw(*),iopt(*),part(*)
    !   real(c_float)  :: weights(*)
    !   !integer(psb_ipk_) :: n,wgflag,numflag,nparts,nedc
    !   !integer(psb_ipk_) :: ixadj(*),iadj(*),ivwg(*),iajw(*),iopt(*),part(*)
    ! end subroutine METIS_PartGraphKway

    function METIS_PartGraphKway(n,ixadj,iadj,ivwg,iajw,&
         & nparts,weights,part) bind(c,name="metis_PartGraphKway_C") result(res)
      use iso_c_binding
      integer(c_int) :: res
      integer(c_int) :: n,nparts
      integer(c_int) :: ixadj(*),iadj(*),ivwg(*),iajw(*),part(*)
      real(c_float)  :: weights(*)
      !integer(psb_ipk_) :: n,wgflag,numflag,nparts,nedc
      !integer(psb_ipk_) :: ixadj(*),iadj(*),ivwg(*),iajw(*),iopt(*),part(*)
    end function METIS_PartGraphKway
  end interface

  call psb_realloc(n,graph_vect,info)
  if (info == psb_success_) allocate(gvl(n),wgh_(nparts),stat=info)

  if (info /= psb_success_) then
    write(psb_err_unit,*) 'Fatal error in BUILD_MTPART: memory allocation ',&
         & ' failure.'
    return
  endif
  if (nparts > 1) then
    iopt(1) = 0
    numflag  = 1
    wgflag   = 0

!!$        write(*,*) 'Before allocation',nparts

    irpl=irp
    jal = ja
    nl = n
    nptl = nparts
    wgh_ = -1.0
    if(present(weights)) then
      if (size(weights) == nptl) then 
!!$            write(*,*) 'weights present',weights
        ! call METIS_PartGraphKway(n,irp,ja,idummy,jdummy,&
        !      & wgflag,numflag,nparts,weights,iopt,nedc,graph_vect)
        info = METIS_PartGraphKway(nl,irpl,jal,idummy,jdummy,&
             & nptl,weights,gvl)

      else
!!$            write(*,*) 'weights absent',wgh_
        info = METIS_PartGraphKway(nl,irpl,jal,idummy,jdummy,&
             & nptl,wgh_,gvl)
      end if
    else
!!$          write(*,*) 'weights absent',wgh_
      info = METIS_PartGraphKway(nl,irpl,jal,idummy,jdummy,&
           & nptl,wgh_,gvl)
    endif
!!$        write(*,*) 'after allocation',info

    do i=1, n
      graph_vect(i) = gvl(i) - 1 
    enddo
  else
    do i=1, n
      graph_vect(i) = 0
    enddo
  endif
#else
  write(psb_err_unit,*) 'Warning: METIS was not configured at PSBLAS compile time !'
#endif

  return

end subroutine psi_build_mtpart
