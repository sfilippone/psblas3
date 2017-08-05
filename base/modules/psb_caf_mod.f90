module psb_caf_mod
use psb_const_mod

logical, parameter :: if_caf=.true.
logical, parameter :: if_caf2=.true.

interface caf_scatterv
  module procedure caf_iscatterv
  module procedure caf_sscatterv
  module procedure caf_dscatterv
  module procedure caf_dscatterv_s
  module procedure caf_cscatterv
  module procedure caf_zscatterv
end interface caf_scatterv

interface caf_gatherv
  module procedure caf_igatherv
  module procedure caf_sgatherv
  module procedure caf_dgatherv
  module procedure caf_cgatherv
  module procedure caf_zgatherv
end interface caf_gatherv

interface caf_allgatherv
  module procedure caf_iallgatherv
  module procedure caf_sallgatherv
  module procedure caf_dallgatherv
  module procedure caf_callgatherv
  module procedure caf_zallgatherv
end interface caf_allgatherv

interface caf_gather
   module procedure caf_igather_s
   module procedure caf_dgather_s
   module procedure caf_sgather_s
   module procedure caf_cgather_s
   module procedure caf_zgather_s
end interface caf_gather

interface caf_alltoallv
   module procedure caf_ialltoallv
   module procedure caf_dalltoallv
   module procedure caf_salltoallv
   module procedure caf_calltoallv
   module procedure caf_zalltoallv
end interface caf_alltoallv

interface caf_amx_reduce
   module procedure caf_camx_reduces
   module procedure caf_zamx_reduces
   module procedure caf_camx_reducev
   module procedure caf_zamx_reducev
   module procedure caf_camx_reducem
   module procedure caf_zamx_reducem
end interface caf_amx_reduce

interface caf_amn_reduce
   module procedure caf_camn_reduces
   module procedure caf_zamn_reduces
   module procedure caf_camn_reducev
   module procedure caf_zamn_reducev
   module procedure caf_camn_reducem
   module procedure caf_zamn_reducem
end interface caf_amn_reduce

contains
  subroutine caf_barrier()
    sync all
  end subroutine caf_barrier


  subroutine caf_iscatterv(snd, scount, rcv, sdispls,rcount,root, info)
    implicit none
    integer(psb_ipk_), intent(in) :: scount(:), snd(:), rcount, sdispls(:), root
    integer(psb_ipk_), allocatable, intent(inout) :: rcv(:)
    integer(psb_ipk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    integer(psb_ipk_), save :: max_rcount[*]
    integer(psb_ipk_), allocatable :: rcv_buf(:)[:]

    np = num_images()
    me = this_image()
    if (size(scount,1) /= np) then
      info = -3
      return
    endif
    if (size(sdispls,1) /= np) then
      info = -4
      return
    endif

    max_rcount=rcount
    call co_max(max_rcount)
    if (allocated(rcv_buf)) deallocate(rcv_buf)
    allocate(rcv_buf(max_rcount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rcount )) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rcount))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif



    if (me == root) then
      do i=1, np
        start = sdispls(i) + 1
        finish= start + scount(i) - 1
        rcv_buf(:)[i]=snd(start:finish)
      enddo
    endif
    sync all
    rcv(1:rcount)=rcv_buf(1:rcount)

    if (allocated(rcv_buf)) deallocate(rcv_buf)
    sync all
  end subroutine caf_iscatterv

  subroutine caf_sscatterv(snd, scount, rcv, sdispls,rcount,root, info)
    implicit none
    integer(psb_ipk_), intent(in) :: scount(:), rcount, sdispls(:), root
    real(psb_spk_) ::  snd(:)
    real(psb_spk_), allocatable, intent(inout) :: rcv(:)
    real(psb_spk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    integer(psb_ipk_), save :: max_rcount[*]
    real(psb_spk_), allocatable :: rcv_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(scount,1) /= np) then
      info = -3
      return
    endif
    if (size(sdispls,1) /= np) then
      info = -4
      return
    endif

    max_rcount=rcount
    call co_max(max_rcount)
    if (allocated(rcv_buf)) deallocate(rcv_buf)
    allocate(rcv_buf(max_rcount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rcount )) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rcount))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif



    if (me == root) then
      do i=1, np
        start = sdispls(i) + 1
        finish= start + scount(i) - 1
        rcv_buf(:)[i]=snd(start:finish)
      enddo
    endif
    sync all
    rcv(1:rcount)=rcv_buf(1:rcount)

    if (allocated(rcv_buf)) deallocate(rcv_buf)
    sync all
  end subroutine caf_sscatterv

  subroutine caf_dscatterv(snd, scount, rcv, sdispls,rcount,root, info)
    implicit none
    integer(psb_ipk_), intent(in) :: scount(:), rcount, sdispls(:), root
    real(psb_dpk_) ::  snd(:)
    real(psb_dpk_), allocatable, intent(inout) :: rcv(:)
    real(psb_dpk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    integer(psb_ipk_), save :: max_rcount[*]
    real(psb_dpk_), allocatable :: rcv_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(scount,1) /= np) then
      info = -3
      return
    endif
    if (size(sdispls,1) /= np) then
      info = -4
      return
    endif

    max_rcount=rcount
    call co_max(max_rcount)
    if (allocated(rcv_buf)) deallocate(rcv_buf)
    allocate(rcv_buf(max_rcount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rcount )) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rcount))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif


    if (me == root) then
      do i=1, np
        start = sdispls(i) + 1
        finish= start + scount(i) - 1
        rcv_buf(:)[i]=snd(start:finish)
      enddo
    endif
    sync all
    rcv(1:rcount)=rcv_buf(1:rcount)

    if (allocated(rcv_buf)) deallocate(rcv_buf)
    sync all
  end subroutine caf_dscatterv

  subroutine caf_dscatterv_s(snd, scount, rcv, sdispls,rcount,root, info)
    implicit none
    integer(psb_ipk_), intent(in) :: scount(:), rcount, sdispls(:), root
    real(psb_dpk_) ::  snd(:)
    real(psb_dpk_), intent(inout) :: rcv
    real(psb_dpk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    integer(psb_ipk_), save :: max_rcount[*]
    real(psb_dpk_), allocatable :: rcv_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(scount,1) /= np) then
      info = -3
      return
    endif
    if (size(sdispls,1) /= np) then
      info = -4
      return
    endif

    max_rcount=rcount
    call co_max(max_rcount)
    if (allocated(rcv_buf)) deallocate(rcv_buf)
    allocate(rcv_buf(max_rcount)[*], STAT=info)
    if (info/=0) return


    if (me == root) then
      do i=1, np
        start = sdispls(i) + 1
        finish= start + scount(i) - 1
        rcv_buf(:)[i]=snd(start:finish)
      enddo
    endif
    sync all
    rcv=rcv_buf(1)

    if (allocated(rcv_buf)) deallocate(rcv_buf)
    sync all
  end subroutine caf_dscatterv_s


  subroutine caf_cscatterv(snd, scount, rcv, sdispls,rcount,root, info)
    implicit none
    integer(psb_ipk_), intent(in) :: scount(:), rcount, sdispls(:), root
    complex(psb_spk_) ::  snd(:)
    complex(psb_spk_), allocatable, intent(inout) :: rcv(:)
    complex(psb_spk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    integer(psb_ipk_), save :: max_rcount[*]
    complex(psb_spk_), allocatable :: rcv_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(scount,1) /= np) then
      info = -3
      return
    endif
    if (size(sdispls,1) /= np) then
      info = -4
      return
    endif

    max_rcount=rcount
    call co_max(max_rcount)
    if (allocated(rcv_buf)) deallocate(rcv_buf)
    allocate(rcv_buf(max_rcount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rcount )) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rcount))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif


    if (me == root) then
      do i=1, np
        start = sdispls(i) + 1
        finish= start + scount(i) - 1
        rcv_buf(:)[i]=snd(start:finish)
      enddo
    endif
    sync all
    rcv(1:rcount)=rcv_buf(1:rcount)

    if (allocated(rcv_buf)) deallocate(rcv_buf)
    sync all
  end subroutine caf_cscatterv


  subroutine caf_zscatterv(snd, scount, rcv, sdispls,rcount,root, info)
    implicit none
    integer(psb_ipk_), intent(in) :: scount(:), rcount, sdispls(:), root
    complex(psb_dpk_) ::  snd(:)
    complex(psb_dpk_), allocatable, intent(inout) :: rcv(:)
    complex(psb_dpk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    integer(psb_ipk_), save :: max_rcount[*]
    complex(psb_dpk_), allocatable :: rcv_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(scount,1) /= np) then
      info = -3
      return
    endif
    if (size(sdispls,1) /= np) then
      info = -4
      return
    endif

    max_rcount=rcount
    call co_max(max_rcount)
    if (allocated(rcv_buf)) deallocate(rcv_buf)
    allocate(rcv_buf(max_rcount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rcount )) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rcount))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif


    if (me == root) then
      do i=1, np
        start = sdispls(i) + 1
        finish= start + scount(i) - 1
        rcv_buf(:)[i]=snd(start:finish)
      enddo
    endif
    sync all
    rcv(1:rcount)=rcv_buf(1:rcount)

    if (allocated(rcv_buf)) deallocate(rcv_buf)
    sync all
  end subroutine caf_zscatterv

  !In this gather subroutine, each image sends just one element (s stands for scalar)
  subroutine caf_igather_s(snd,rcv,root, info)
    integer(psb_ipk_), intent(in) :: snd
    integer(psb_ipk_), intent(inout) :: rcv(:)
    integer(psb_ipk_), intent(inout):: info
    integer(psb_ipk_), intent(in) :: root
    !Local
    integer(psb_ipk_), allocatable :: rcv_buf(:)[:]
    integer(psb_ipk_) :: me, np
    np=num_images()
    me = this_image()
    allocate (rcv_buf(np)[*], STAT=info)
    if (info /= 0) then
       print*,'allocation error', info
       return
    endif
    if (me == root) then
      rcv_buf(1:np)=rcv(1:np)
    endif
    sync all
    rcv_buf(me)[root] = snd
    sync all
    if (me == root) then
      rcv(1:np) = rcv_buf(1:np)
    endif
  end subroutine caf_igather_s

  !In this gather subroutine, each image sends just one element (s stands for scalar)
  subroutine caf_dgather_s(snd,rcv,root, info)
    real(psb_dpk_), intent(in) :: snd
    real(psb_dpk_), intent(inout) :: rcv(:)
    integer(psb_ipk_), intent(inout):: info
    integer(psb_ipk_), intent(in) :: root
    !Local
    real(psb_dpk_), allocatable :: rcv_buf(:)[:]
    integer(psb_ipk_) :: me, np
    np=num_images()
    me = this_image()
    allocate (rcv_buf(np)[*], STAT=info)
    if (info /= 0) then
       print*,'allocation error', info
       return
    endif
    if (me == root) then
      rcv_buf(1:np)=rcv(1:np)
    endif
    sync all
    rcv_buf(me)[root] = snd
    sync all
    if (me == root) then
      rcv(1:np) = rcv_buf(1:np)
    endif
  end subroutine caf_dgather_s

  !In this gather subroutine, each image sends just one element (s stands for scalar)
  subroutine caf_sgather_s(snd,rcv,root, info)
    real(psb_spk_), intent(in) :: snd
    real(psb_spk_), intent(inout) :: rcv(:)
    integer(psb_ipk_), intent(inout):: info
    integer(psb_ipk_), intent(in) :: root
    !Local
    real(psb_spk_), allocatable :: rcv_buf(:)[:]
    integer(psb_ipk_) :: me, np

    np=num_images()
    me = this_image()
    allocate (rcv_buf(np)[*], STAT=info)
    if (info /= 0) then
       print*,'allocation error', info
       return
    endif
    if (me == root) then
      rcv_buf(1:np)=rcv(1:np)
    endif
    sync all
    rcv_buf(me)[root] = snd
    sync all
    if (me == root) then
      rcv(1:np) = rcv_buf(1:np)
    endif
  end subroutine caf_sgather_s

  !In this gather subroutine, each image sends just one element (s stands for scalar)
  subroutine caf_zgather_s(snd,rcv,root, info)
    complex(psb_dpk_), intent(in) :: snd
    complex(psb_dpk_), intent(inout) :: rcv(:)
    integer(psb_ipk_), intent(inout):: info
    integer(psb_ipk_), intent(in) :: root
    !Local
    complex(psb_dpk_), allocatable :: rcv_buf(:)[:]
    integer(psb_ipk_) :: me, np

    np=num_images()
    me = this_image()
    allocate (rcv_buf(np)[*], STAT=info)
    if (info /= 0) then
       print*,'allocation error', info
       return
    endif
    if (me == root) then
      rcv_buf(1:np)=rcv(1:np)
    endif
    sync all
    rcv_buf(me)[root] = snd
    sync all
    if (me == root) then
      rcv(1:np) = rcv_buf(1:np)
    endif
  end subroutine caf_zgather_s


  !In this gather subroutine, each image sends just one element (s stands for scalar)
  subroutine caf_cgather_s(snd,rcv,root, info)
    complex(psb_spk_), intent(in) :: snd
    complex(psb_spk_), intent(inout) :: rcv(:)
    integer(psb_ipk_), intent(inout):: info
    integer(psb_ipk_), intent(in) :: root
    !Local
    complex(psb_spk_), allocatable :: rcv_buf(:)[:]
    integer(psb_ipk_) :: me, np

    np=num_images()
    me = this_image()
    allocate (rcv_buf(np)[*], STAT=info)
    if (info /= 0) then
       print*,'allocation error', info
       return
    endif
    if (me == root) then
      rcv_buf(1:np)=rcv(1:np)
    endif
    sync all
    rcv_buf(me)[root] = snd
    sync all
    if (me == root) then
      rcv(1:np) = rcv_buf(1:np)
    endif
  end subroutine caf_cgather_s


  subroutine caf_alltoall(snd,rcv, m, info)
    use mpi
    implicit none
    integer(psb_ipk_), intent(in) :: m
    integer(psb_ipk_), intent(in) :: snd(:)
    integer(psb_ipk_), intent(out):: rcv(:), info
    !Local
    integer(psb_ipk_) :: me, np,i, snd_start, snd_finish, snd_tot
    integer(psb_ipk_), allocatable :: buffer(:)[:], rcv_start, rcv_finish
    double precision :: t1, t2
t1 = mpi_wtime()
    if ( m < 0) then
      print*,'Error, m must be greater or equal to zero'
      info = -1
    endif
    me = this_image()
    np = num_images()
    snd_tot=m*np
    if (allocated(buffer)) deallocate(buffer)
    allocate(buffer(snd_tot)[*], STAT = info)
    if (info /= 0) then
       print*,'allocation error', info
       return
    endif
    rcv_start = (me-1)*m +1
    rcv_finish = rcv_start + m - 1
    do i=1,np
      snd_start = (i-1)*m +1 
      snd_finish = snd_start + m - 1
      if (rcv_finish > snd_tot) then
        print*,'Error, rcv_finish > snd_tot'
        info = -2 
	return
      endif
      if (snd_finish > size(snd,1)) then
        print*,'Error, snd_finish > size(snd,1)'
        info = -3 
      endif
      buffer(rcv_start:rcv_finish)[i]=snd(snd_start:snd_finish)
    enddo
    !Copy buffer
    sync all
    rcv(1:snd_tot)=buffer(1:snd_tot)
    if (allocated(buffer)) deallocate(buffer)
    sync all
t2 = mpi_wtime() - t1
  end subroutine caf_alltoall

  subroutine caf_ialltoallv(snd, scount, sdispl, rcv, rcount, rdispl, info)
    use mpi
    implicit none
    integer(psb_ipk_), intent(in) :: snd(:)
    integer(psb_ipk_), intent(out):: rcv(:)
    integer(psb_ipk_), intent(in) :: scount(:), sdispl(:) 
    integer(psb_ipk_), intent(in) :: rcount(:), rdispl(:)
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), allocatable :: buffer(:)[:], buffer_rcount(:)[:], buffer_rdispl(:)[:]
    integer(psb_ipk_) :: np, i, size_, rcv_start, rcv_finish, snd_start, snd_finish, me
    type(event_type), allocatable :: done(:)[:], alltoall(:)[:]
    double precision :: t1, t2

t1 = mpi_wtime()
    np=num_images()
    me=this_image()
    size_=size(rcv,1)

    if (allocated(buffer_rdispl)) deallocate(buffer_rdispl)
    if (allocated(buffer_rcount)) deallocate(buffer_rcount)
    if (allocated(buffer)) deallocate(buffer)
    if (allocated(done)) deallocate(done)
    if (allocated(alltoall)) deallocate(alltoall)

    allocate(buffer_rdispl(np)[*],buffer_rcount(np)[*],buffer(size_)[*],done(np)[*],alltoall(np)[*], STAT=info)

    if (info /= 0) then
       print*,'Allocation error'
       return
    endif

    ! All to all rcount e rdispl
    rcv_start = me
    do i=1,np
      snd_start = i
      buffer_rdispl(rcv_start)[i]=rdispl(snd_start)
      buffer_rcount(rcv_start)[i]=rcount(snd_start)
      event post(alltoall(me)[i])
    enddo


    if (info /= 0) then
       return
    endif

    do i=1,np
      snd_start = sdispl(i) + 1
      snd_finish = snd_start + scount(i) - 1
      event wait(alltoall(i))
      rcv_start = buffer_rdispl(i) + 1
      rcv_finish = rcv_start + buffer_rcount(i) - 1
      buffer(rcv_start:rcv_finish)[i]=snd(snd_start:snd_finish)
      event post(done(me)[i])
    enddo

    do i=1,np
      rcv_start = rdispl(i) + 1
      rcv_finish = rcv_start + rcount(i) - 1
      event wait(done(i))
      rcv(rcv_start:rcv_finish)=buffer(rcv_start:rcv_finish)
    enddo

    if (allocated(buffer_rdispl)) deallocate(buffer_rdispl)
    if (allocated(buffer_rcount)) deallocate(buffer_rcount)
    if (allocated(buffer)) deallocate(buffer)
    if (allocated(done)) deallocate(done)
    if (allocated(alltoall)) deallocate(alltoall)
    sync all
t2 = mpi_wtime() - t1
  end subroutine caf_ialltoallv

  subroutine caf_dalltoallv(snd, scount, sdispl, rcv, rcount, rdispl, info)
    implicit none
    real(psb_dpk_), intent(in) :: snd(:)
    real(psb_dpk_), intent(out):: rcv(:)
    integer(psb_ipk_), intent(in) :: scount(:), sdispl(:) 
    integer(psb_ipk_), intent(in) :: rcount(:), rdispl(:)
    integer(psb_ipk_), intent(out) :: info
    real(psb_dpk_), allocatable :: buffer(:)[:]
    integer(psb_ipk_), allocatable :: buffer_rcount(:)[:], buffer_rdispl(:)[:]
    integer(psb_ipk_) :: np, i, size_, rcv_start, rcv_finish, snd_start, snd_finish, me
    type(event_type), allocatable :: done(:)[:], alltoall(:)[:]

    np=num_images()
    me=this_image()
    size_=size(rcv,1)

    if (allocated(buffer_rdispl)) deallocate(buffer_rdispl)
    if (allocated(buffer_rcount)) deallocate(buffer_rcount)
    if (allocated(buffer)) deallocate(buffer)
    if (allocated(done)) deallocate(done)
    if (allocated(alltoall)) deallocate(alltoall)

    allocate(buffer_rdispl(np)[*],buffer_rcount(np)[*],buffer(size_)[*],done(np)[*],alltoall(np)[*], STAT=info)

    if (info /= 0) then
       print*,'Allocation error'
       return
    endif

    ! All to all rcount e rdispl
    rcv_start = me
    do i=1,np
      snd_start = i
      buffer_rdispl(rcv_start)[i]=rdispl(snd_start)
      buffer_rcount(rcv_start)[i]=rcount(snd_start)
      event post(alltoall(me)[i])
    enddo


    if (info /= 0) then
       return
    endif

    do i=1,np
      snd_start = sdispl(i) + 1
      snd_finish = snd_start + scount(i) - 1
      event wait(alltoall(i))
      rcv_start = buffer_rdispl(i) + 1
      rcv_finish = rcv_start + buffer_rcount(i) - 1
      buffer(rcv_start:rcv_finish)[i]=snd(snd_start:snd_finish)
      event post(done(me)[i])
    enddo

    do i=1,np
      rcv_start = rdispl(i) + 1
      rcv_finish = rcv_start + rcount(i) - 1
      event wait(done(i))
      rcv(rcv_start:rcv_finish)=buffer(rcv_start:rcv_finish)
    enddo
    if (allocated(buffer_rdispl)) deallocate(buffer_rdispl)
    if (allocated(buffer_rcount)) deallocate(buffer_rcount)
    if (allocated(buffer)) deallocate(buffer)
    if (allocated(done)) deallocate(done)
    if (allocated(alltoall)) deallocate(alltoall)
    sync all
  end subroutine caf_dalltoallv

  subroutine caf_salltoallv(snd, scount, sdispl, rcv, rcount, rdispl, info)
    implicit none
    real(psb_spk_), intent(in) :: snd(:)
    real(psb_spk_), intent(out):: rcv(:)
    integer(psb_ipk_), intent(in) :: scount(:), sdispl(:) 
    integer(psb_ipk_), intent(in) :: rcount(:), rdispl(:)
    integer(psb_ipk_), intent(out) :: info
    real(psb_spk_), allocatable :: buffer(:)[:]
    integer(psb_ipk_), allocatable :: buffer_rcount(:)[:], buffer_rdispl(:)[:]
    integer(psb_ipk_) :: np, i, size_, rcv_start, rcv_finish, snd_start, snd_finish, me
    type(event_type), allocatable :: done(:)[:], alltoall(:)[:]

    np=num_images()
    me=this_image()
    size_=size(rcv,1)

    if (allocated(buffer_rdispl)) deallocate(buffer_rdispl)
    if (allocated(buffer_rcount)) deallocate(buffer_rcount)
    if (allocated(buffer)) deallocate(buffer)
    if (allocated(done)) deallocate(done)
    if (allocated(alltoall)) deallocate(alltoall)

    allocate(buffer_rdispl(np)[*],buffer_rcount(np)[*],buffer(size_)[*],done(np)[*],alltoall(np)[*], STAT=info)

    if (info /= 0) then
       print*,'Allocation error'
       return
    endif

    ! All to all rcount e rdispl
    rcv_start = me
    do i=1,np
      snd_start = i
      buffer_rdispl(rcv_start)[i]=rdispl(snd_start)
      buffer_rcount(rcv_start)[i]=rcount(snd_start)
      event post(alltoall(me)[i])
    enddo


    if (info /= 0) then
       return
    endif

    do i=1,np
      snd_start = sdispl(i) + 1
      snd_finish = snd_start + scount(i) - 1
      event wait(alltoall(i))
      rcv_start = buffer_rdispl(i) + 1
      rcv_finish = rcv_start + buffer_rcount(i) - 1
      buffer(rcv_start:rcv_finish)[i]=snd(snd_start:snd_finish)
      event post(done(me)[i])
    enddo

    do i=1,np
      rcv_start = rdispl(i) + 1
      rcv_finish = rcv_start + rcount(i) - 1
      event wait(done(i))
      rcv(rcv_start:rcv_finish)=buffer(rcv_start:rcv_finish)
    enddo
    if (allocated(buffer_rdispl)) deallocate(buffer_rdispl)
    if (allocated(buffer_rcount)) deallocate(buffer_rcount)
    if (allocated(buffer)) deallocate(buffer)
    if (allocated(done)) deallocate(done)
    if (allocated(alltoall)) deallocate(alltoall)
    sync all
  end subroutine caf_salltoallv

  subroutine caf_calltoallv(snd, scount, sdispl, rcv, rcount, rdispl, info)
    implicit none
    complex(psb_spk_), intent(in) :: snd(:)
    complex(psb_spk_), intent(out):: rcv(:)
    integer(psb_ipk_), intent(in) :: scount(:), sdispl(:) 
    integer(psb_ipk_), intent(in) :: rcount(:), rdispl(:)
    integer(psb_ipk_), intent(out) :: info
    complex(psb_spk_), allocatable :: buffer(:)[:]
    integer(psb_ipk_), allocatable :: buffer_rcount(:)[:], buffer_rdispl(:)[:]
    integer(psb_ipk_) :: np, i, size_, rcv_start, rcv_finish, snd_start, snd_finish, me
    type(event_type), allocatable :: done(:)[:], alltoall(:)[:]

    np=num_images()
    me=this_image()
    size_=size(rcv,1)

    if (allocated(buffer_rdispl)) deallocate(buffer_rdispl)
    if (allocated(buffer_rcount)) deallocate(buffer_rcount)
    if (allocated(buffer)) deallocate(buffer)
    if (allocated(done)) deallocate(done)
    if (allocated(alltoall)) deallocate(alltoall)

    allocate(buffer_rdispl(np)[*],buffer_rcount(np)[*],buffer(size_)[*],done(np)[*],alltoall(np)[*], STAT=info)

    if (info /= 0) then
       print*,'Allocation error'
       return
    endif

    ! All to all rcount e rdispl
    rcv_start = me
    do i=1,np
      snd_start = i
      buffer_rdispl(rcv_start)[i]=rdispl(snd_start)
      buffer_rcount(rcv_start)[i]=rcount(snd_start)
      event post(alltoall(me)[i])
    enddo


    if (info /= 0) then
       return
    endif

    do i=1,np
      snd_start = sdispl(i) + 1
      snd_finish = snd_start + scount(i) - 1
      event wait(alltoall(i))
      rcv_start = buffer_rdispl(i) + 1
      rcv_finish = rcv_start + buffer_rcount(i) - 1
      buffer(rcv_start:rcv_finish)[i]=snd(snd_start:snd_finish)
      event post(done(me)[i])
    enddo

    do i=1,np
      rcv_start = rdispl(i) + 1
      rcv_finish = rcv_start + rcount(i) - 1
      event wait(done(i))
      rcv(rcv_start:rcv_finish)=buffer(rcv_start:rcv_finish)
    enddo
    if (allocated(buffer_rdispl)) deallocate(buffer_rdispl)
    if (allocated(buffer_rcount)) deallocate(buffer_rcount)
    if (allocated(buffer)) deallocate(buffer)
    if (allocated(done)) deallocate(done)
    if (allocated(alltoall)) deallocate(alltoall)
    sync all
  end subroutine caf_calltoallv

  subroutine caf_zalltoallv(snd, scount, sdispl, rcv, rcount, rdispl, info)
    implicit none
    complex(psb_dpk_), intent(in) :: snd(:)
    complex(psb_dpk_), intent(out):: rcv(:)
    integer(psb_ipk_), intent(in) :: scount(:), sdispl(:) 
    integer(psb_ipk_), intent(in) :: rcount(:), rdispl(:)
    integer(psb_ipk_), intent(out) :: info
    complex(psb_dpk_), allocatable :: buffer(:)[:]
    integer(psb_ipk_), allocatable :: buffer_rcount(:)[:], buffer_rdispl(:)[:]
    integer(psb_ipk_) :: np, i, size_, rcv_start, rcv_finish, snd_start, snd_finish, me
    type(event_type), allocatable :: done(:)[:], alltoall(:)[:]

    np=num_images()
    me=this_image()
    size_=size(rcv,1)

    if (allocated(buffer_rdispl)) deallocate(buffer_rdispl)
    if (allocated(buffer_rcount)) deallocate(buffer_rcount)
    if (allocated(buffer)) deallocate(buffer)
    if (allocated(done)) deallocate(done)
    if (allocated(alltoall)) deallocate(alltoall)

    allocate(buffer_rdispl(np)[*],buffer_rcount(np)[*],buffer(size_)[*],done(np)[*],alltoall(np)[*], STAT=info)

    if (info /= 0) then
       print*,'Allocation error'
       return
    endif

    ! All to all rcount e rdispl
    rcv_start = me
    do i=1,np
      snd_start = i
      buffer_rdispl(rcv_start)[i]=rdispl(snd_start)
      buffer_rcount(rcv_start)[i]=rcount(snd_start)
      event post(alltoall(me)[i])
    enddo


    if (info /= 0) then
       return
    endif

    do i=1,np
      snd_start = sdispl(i) + 1
      snd_finish = snd_start + scount(i) - 1
      event wait(alltoall(i))
      rcv_start = buffer_rdispl(i) + 1
      rcv_finish = rcv_start + buffer_rcount(i) - 1
      buffer(rcv_start:rcv_finish)[i]=snd(snd_start:snd_finish)
      event post(done(me)[i])
    enddo

    do i=1,np
      rcv_start = rdispl(i) + 1
      rcv_finish = rcv_start + rcount(i) - 1
      event wait(done(i))
      rcv(rcv_start:rcv_finish)=buffer(rcv_start:rcv_finish)
    enddo
    if (allocated(buffer_rdispl)) deallocate(buffer_rdispl)
    if (allocated(buffer_rcount)) deallocate(buffer_rcount)
    if (allocated(buffer)) deallocate(buffer)
    if (allocated(done)) deallocate(done)
    if (allocated(alltoall)) deallocate(alltoall)
    sync all
  end subroutine caf_zalltoallv

  subroutine caf_igatherv(snd, scount, rcv, rcount, rdispls, root, info)
    implicit none
    integer(psb_ipk_), intent(in) :: scount, snd(:), rcount(:), rdispls(:), root
    integer(psb_ipk_), allocatable, intent(inout) :: rcv(:)
    integer(psb_ipk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    integer(psb_ipk_), allocatable :: snd_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(rcount,1) /= np) then
      info = -3
      return
    endif
    if (size(rdispls,1) /= np) then
      info = -4
      return
    endif

    if (allocated(snd_buf)) deallocate(snd_buf)
    allocate(snd_buf(scount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rdispls(np) + rcount(np))) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rdispls(np)+rcount(np)))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif
    snd_buf(1:scount)=snd(1:scount)

    sync all

    if (me == root) then
      do i=1, np
        start = rdispls(i) + 1
        finish= start + rcount(i) - 1
        rcv(start:finish)=snd_buf(1:rcount(i))[i]
      enddo
    endif

    if (allocated(snd_buf)) deallocate(snd_buf)
    sync all
  end subroutine caf_igatherv

  subroutine caf_sgatherv(snd, scount, rcv, rcount, rdispls, root, info)
    implicit none
    real(psb_spk_), intent(in) :: snd(:)
    integer(psb_ipk_), intent(in) :: scount, rcount(:), rdispls(:)
    real(psb_spk_), allocatable, intent(inout) :: rcv(:)
    real(psb_spk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), intent(in) :: root
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    real(psb_spk_), allocatable :: snd_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(rcount,1) /= np) then
      info = -3
      return
    endif
    if (size(rdispls,1) /= np) then
      info = -4
      return
    endif

    if (allocated(snd_buf)) deallocate(snd_buf)
    allocate(snd_buf(scount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rdispls(np) + rcount(np))) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rdispls(np)+rcount(np)))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif
    snd_buf(1:scount)=snd(1:scount)

    sync all

    if (me == root) then
      do i=1, np
        start = rdispls(i) + 1
        finish= start + rcount(i) - 1
        rcv(start:finish)=snd_buf(1:rcount(i))[i]
      enddo
    endif

    if (allocated(snd_buf)) deallocate(snd_buf)
    sync all
  end subroutine caf_sgatherv

  subroutine caf_dgatherv(snd, scount, rcv, rcount, rdispls, root, info)
    implicit none
    real(psb_dpk_), intent(in) :: snd(:)
    integer(psb_ipk_), intent(in) :: scount, rcount(:), rdispls(:)
    real(psb_dpk_), allocatable, intent(inout) :: rcv(:)
    real(psb_dpk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), intent(in) :: root
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    real(psb_dpk_), allocatable :: snd_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(rcount,1) /= np) then
      info = -3
      return
    endif
    if (size(rdispls,1) /= np) then
      info = -4
      return
    endif

    if (allocated(snd_buf)) deallocate(snd_buf)
    allocate(snd_buf(scount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rdispls(np) + rcount(np))) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rdispls(np)+rcount(np)))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif
    snd_buf(1:scount)=snd(1:scount)

    sync all

    if (me == root) then
      do i=1, np
        start = rdispls(i) + 1
        finish= start + rcount(i) - 1
        rcv(start:finish)=snd_buf(1:rcount(i))[i]
      enddo
    endif

    if (allocated(snd_buf)) deallocate(snd_buf)
    sync all
  end subroutine caf_dgatherv

  subroutine caf_cgatherv(snd, scount, rcv, rcount, rdispls, root, info)
    implicit none
    complex(psb_spk_), intent(in) :: snd(:)
    integer(psb_ipk_), intent(in) :: scount, rcount(:), rdispls(:)
    complex(psb_spk_), allocatable, intent(inout) :: rcv(:)
    complex(psb_spk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), intent(in) :: root
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    complex(psb_spk_), allocatable :: snd_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(rcount,1) /= np) then
      info = -3
      return
    endif
    if (size(rdispls,1) /= np) then
      info = -4
      return
    endif

    if (allocated(snd_buf)) deallocate(snd_buf)
    allocate(snd_buf(scount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rdispls(np) + rcount(np))) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rdispls(np)+rcount(np)))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif
    snd_buf(1:scount)=snd(1:scount)

    sync all

    if (me == root) then
      do i=1, np
        start = rdispls(i) + 1
        finish= start + rcount(i) - 1
        rcv(start:finish)=snd_buf(1:rcount(i))[i]
      enddo
    endif

    if (allocated(snd_buf)) deallocate(snd_buf)
    sync all
  end subroutine caf_cgatherv

  subroutine caf_zgatherv(snd, scount, rcv, rcount, rdispls, root, info)
    implicit none
    complex(psb_dpk_), intent(in) :: snd(:)
    integer(psb_ipk_), intent(in) :: scount, rcount(:), rdispls(:)
    complex(psb_dpk_), allocatable, intent(inout) :: rcv(:)
    complex(psb_dpk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    integer(psb_ipk_), intent(in) :: root   
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    complex(psb_dpk_), allocatable :: snd_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(rcount,1) /= np) then
      info = -3
      return
    endif
    if (size(rdispls,1) /= np) then
      info = -4
      return
    endif

    if (allocated(snd_buf)) deallocate(snd_buf)
    allocate(snd_buf(scount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rdispls(np) + rcount(np))) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rdispls(np)+rcount(np)))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif
    snd_buf(1:scount)=snd(1:scount)

    sync all

    if (me == root) then
      do i=1, np
        start = rdispls(i) + 1
        finish= start + rcount(i) - 1
        rcv(start:finish)=snd_buf(1:rcount(i))[i]
      enddo
    endif
    if (allocated(snd_buf)) deallocate(snd_buf)
    sync all
  end subroutine caf_zgatherv

  subroutine caf_iallgatherv(snd, scount, rcv, rcount, rdispls, info)
    use psb_realloc_mod
    implicit none
    integer(psb_ipk_), intent(in) :: scount, snd(:), rcount(:), rdispls(:)
    integer(psb_ipk_), allocatable, intent(inout) :: rcv(:)
    integer(psb_ipk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    integer(psb_ipk_), allocatable :: snd_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(rcount,1) /= np) then
      info = -3
      return
    endif
    if (size(rdispls,1) /= np) then
      info = -4
      return
    endif
!   write(0,*) 'caf_allgatherv',me,np
!   call flush(0)
!   sync all 
    if (allocated(snd_buf)) deallocate(snd_buf)
    allocate(snd_buf(scount)[*], STAT=info)
    if (info/=0)then
      write(0,*) me, 'Info on allocate ',info
      return
    end if

    if (allocated(rcv).and.(size(rcv,1) < rdispls(np) + rcount(np))) then
      call psb_realloc(rdispls(np)+rcount(np),rcv,info)
!!$      call move_alloc(rcv,rcv_tmp)
!!$      allocate(rcv(rdispls(np)+rcount(np)))
!!$      rcv(1:size(rcv_tmp,1))=rcv_tmp
!!$      deallocate(rcv_tmp)
    endif
    snd_buf(1:scount)=snd(1:scount)

    sync all

    do i=1, np
      start = rdispls(i) + 1
      finish= start + rcount(i) - 1
      rcv(start:finish)=snd_buf(1:rcount(i))[i]
    enddo

    if (allocated(snd_buf)) deallocate(snd_buf)
    sync all
  end subroutine caf_iallgatherv

  subroutine caf_sallgatherv(snd, scount, rcv, rcount, rdispls, info)
    implicit none
    real(psb_spk_), intent(in) :: snd(:)
    integer(psb_ipk_), intent(in) :: scount, rcount(:), rdispls(:)
    real(psb_spk_), allocatable, intent(inout) :: rcv(:)
    real(psb_spk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    real(psb_spk_), allocatable :: snd_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(rcount,1) /= np) then
      info = -3
      return
    endif
    if (size(rdispls,1) /= np) then
      info = -4
      return
    endif

    if (allocated(snd_buf)) deallocate(snd_buf)
    allocate(snd_buf(scount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rdispls(np) + rcount(np))) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rdispls(np)+rcount(np)))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif
    snd_buf(1:scount)=snd(1:scount)

    sync all

    do i=1, np
      start = rdispls(i) + 1
      finish= start + rcount(i) - 1
      rcv(start:finish)=snd_buf(1:rcount(i))[i]
    enddo

    if (allocated(snd_buf)) deallocate(snd_buf)
    sync all
  end subroutine caf_sallgatherv

  subroutine caf_dallgatherv(snd, scount, rcv, rcount, rdispls, info)
    implicit none
    real(psb_dpk_), intent(in) :: snd(:)
    integer(psb_ipk_), intent(in) :: scount, rcount(:), rdispls(:)
    real(psb_dpk_), allocatable, intent(inout) :: rcv(:)
    real(psb_dpk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    real(psb_dpk_), allocatable :: snd_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(rcount,1) /= np) then
      info = -3
      return
    endif
    if (size(rdispls,1) /= np) then
      info = -4
      return
    endif

    if (allocated(snd_buf)) deallocate(snd_buf)
    allocate(snd_buf(scount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rdispls(np) + rcount(np))) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rdispls(np)+rcount(np)))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif
    snd_buf(1:scount)=snd(1:scount)

    sync all

    do i=1, np
      start = rdispls(i) + 1
      finish= start + rcount(i) - 1
      rcv(start:finish)=snd_buf(1:rcount(i))[i]
    enddo

    if (allocated(snd_buf)) deallocate(snd_buf)
    sync all
  end subroutine caf_dallgatherv

  subroutine caf_callgatherv(snd, scount, rcv, rcount, rdispls, info)
    implicit none
    complex(psb_spk_), intent(in) :: snd(:)
    integer(psb_ipk_), intent(in) :: scount, rcount(:), rdispls(:)
    complex(psb_spk_), allocatable, intent(inout) :: rcv(:)
    complex(psb_spk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    complex(psb_spk_), allocatable :: snd_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(rcount,1) /= np) then
      info = -3
      return
    endif
    if (size(rdispls,1) /= np) then
      info = -4
      return
    endif

    if (allocated(snd_buf)) deallocate(snd_buf)
    allocate(snd_buf(scount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rdispls(np) + rcount(np))) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rdispls(np)+rcount(np)))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif
    snd_buf(1:scount)=snd(1:scount)

    sync all

    do i=1, np
      start = rdispls(i) + 1
      finish= start + rcount(i) - 1
      rcv(start:finish)=snd_buf(1:rcount(i))[i]
    enddo

    if (allocated(snd_buf)) deallocate(snd_buf)
    sync all
  end subroutine caf_callgatherv

  subroutine caf_zallgatherv(snd, scount, rcv, rcount, rdispls, info)
    implicit none
    complex(psb_dpk_), intent(in) :: snd(:)
    integer(psb_ipk_), intent(in) :: scount, rcount(:), rdispls(:)
    complex(psb_dpk_), allocatable, intent(inout) :: rcv(:)
    complex(psb_dpk_), allocatable :: rcv_tmp(:)
    integer(psb_ipk_), intent(out) :: info
    ! ---- local variables ---
    integer(psb_ipk_) :: i, np, me, start, finish, snd_tot
    complex(psb_dpk_), allocatable :: snd_buf(:)[:]
    np = num_images()
    me = this_image()
    if (size(rcount,1) /= np) then
      info = -3
      return
    endif
    if (size(rdispls,1) /= np) then
      info = -4
      return
    endif

    if (allocated(snd_buf)) deallocate(snd_buf)
    allocate(snd_buf(scount)[*], STAT=info)
    if (info/=0) return

    if (allocated(rcv).and.(size(rcv,1) < rdispls(np) + rcount(np))) then
      call move_alloc(rcv,rcv_tmp)
      allocate(rcv(rdispls(np)+rcount(np)))
      rcv(1:size(rcv_tmp,1))=rcv_tmp
      deallocate(rcv_tmp)
    endif
    snd_buf(1:scount)=snd(1:scount)

    sync all

    do i=1, np
      start = rdispls(i) + 1
      finish= start + rcount(i) - 1
      rcv(start:finish)=snd_buf(1:rcount(i))[i]
    enddo

    if (allocated(snd_buf)) deallocate(snd_buf)
    sync all
  end subroutine caf_zallgatherv

  subroutine caf_allgather(snd, scount, rcv, info)
    use mpi
    implicit none
    integer(psb_ipk_), intent(in) :: scount, snd(:)
    integer(psb_ipk_), intent(inout) :: rcv(:,:)
    integer(psb_ipk_), intent(out) :: info
    ! ---- local variables ---
    integer(psb_ipk_) :: img,np, me, i
    integer(psb_ipk_), allocatable :: snd_buf(:)[:]
    type(event_type), allocatable :: snd_copied(:)[:]
    double precision :: t1, t2
    t1 = mpi_wtime()
    np = num_images()
    me = this_image()
    info = 0
    if (size(rcv,1) < scount) then
      info = -3
      print*,'error', info, size(rcv,1), scount
      return
    endif
    if (size(rcv,2) < np) then
      info = -4
      print*,'error', info, size(rcv,2), np
      return
    endif
!   write(*,*) 'Hello from ',me,' of:', np,size(snd,1)
!   sync all 
    if (.true.) then 
!!$    if (allocated(snd_buf)) deallocate(snd_buf)
!!$    if (allocated(snd_copied)) deallocate(snd_copied)
    allocate(snd_buf(size(snd,1))[*],stat=info)
    if (info /= 0) then
      write(*,*) 'Error on allocating snd_buf ',me
      stop
    end if
    allocate(snd_copied(np)[*],stat=info) 
    if (info /= 0) then
      write(*,*) 'Error on allocating snd_copied ',me
      stop
    end if
    !allocate(snd_buf(size(snd,1))[*]) 
    snd_buf(:)=snd(:)
!   write(*,*) 'Sending  from ',me,':', snd(:)
    do img=1,np
      event post(snd_copied(me)[img])
    enddo
    !sync all
    if(.false.) then 
      do img=1,np
        event wait(snd_copied(img))
        rcv(:,img)=snd_buf(:)[img]
      enddo
    else
      do img=me+1,np
        event wait(snd_copied(img))
        rcv(:,img)=snd_buf(:)[img]
      enddo
      do img=1,me
        event wait(snd_copied(img))
        rcv(:,img)=snd_buf(:)[img]
      enddo

    end if
    !if (allocated(snd_buf))
    deallocate(snd_buf)
    !if (allocated(snd_copied))
    deallocate(snd_copied)
    !Not sure this is necessary...
    !sync all
  else if(.false.) then 
    do img=1,np
      if (me == img) rcv(:,img) = snd(:)
      call co_broadcast(rcv(:,img),img)
    end do
  else if(.false.) then
    rcv(:,:) = 0
    rcv(:,me) = snd(:)
    call co_sum(rcv)
  end if
    t2 = mpi_wtime() - t1
  end subroutine caf_allgather


  pure integer(psb_ipk_) function caf_iamx(a, b)
    implicit none
    integer(psb_ipk_), intent(in)    :: a
    integer(psb_ipk_), intent(in) :: b
    integer(psb_ipk_) :: i
    integer(psb_ipk_) :: w, z
    w = abs( a )
    z = abs( b )
    if ( w>z ) then
      caf_iamx = a
    else
      caf_iamx = b
    end if
  end function caf_iamx

  pure real(psb_spk_) function caf_samx(a, b)
    implicit none
    real(psb_spk_), intent(in)    :: a
    real(psb_spk_), intent(in) :: b
    real(psb_spk_) :: w, z
    w = abs( a )
    z = abs( b )
    if ( w>z ) then
      caf_samx = a
    else
      caf_samx = b
    end if
  end function caf_samx

  pure real(psb_dpk_) function caf_damx(a, b)
    implicit none
    real(psb_dpk_), intent(in)    :: a
    real(psb_dpk_), intent(in) :: b
    real(psb_dpk_) :: w, z
    w = abs( a )
    z = abs( b )
    if ( w>z ) then
      caf_damx = a
    else
      caf_damx = b
    end if
  end function caf_damx


  pure complex(psb_spk_) function caf_camx(a, b)
    implicit none
    complex(psb_spk_), intent(in)    :: a
    complex(psb_spk_), intent(in) :: b
    real(psb_spk_) :: w, z
    w = abs( a )
    z = abs( b )
    if ( w>z ) then
      caf_camx = a
    else
      caf_camx = b
    end if
  end function caf_camx

  pure complex(psb_dpk_) function caf_zamx(a, b)
    implicit none
    complex(psb_dpk_), intent(in)    :: a
    complex(psb_dpk_), intent(in) :: b
    real(psb_dpk_) :: w, z
    w = abs( a )
    z = abs( b )
    if ( w>z ) then
      caf_zamx = a
    else
      caf_zamx = b
    end if
  end function caf_zamx

  pure integer(psb_ipk_) function caf_iamn(a, b)
    implicit none
    integer(psb_ipk_), intent(in)    :: a
    integer(psb_ipk_), intent(in) :: b
    integer(psb_ipk_) :: i
    integer(psb_ipk_) :: w, z
    w = abs( a )
    z = abs( b )
    if ( w<z ) then
      caf_iamn = a
    else
      caf_iamn = b
    end if
  end function caf_iamn

  pure real(psb_spk_) function caf_samn(a, b)
    implicit none
    real(psb_spk_), intent(in)    :: a
    real(psb_spk_), intent(in) :: b
    real(psb_spk_) :: w, z
    w = abs( a )
    z = abs( b )
    if ( w<z ) then
      caf_samn = a
    else
      caf_samn = b
    end if
  end function caf_samn

  pure real(psb_dpk_) function caf_damn(a, b)
    implicit none
    real(psb_dpk_), intent(in)    :: a
    real(psb_dpk_), intent(in) :: b
    real(psb_dpk_) :: w, z
    w = abs( a )
    z = abs( b )
    if ( w<z ) then
      caf_damn = a
    else
      caf_damn = b
    end if
  end function caf_damn


  pure real(psb_dpk_) function caf_dnrm2(vin, vinout)
    implicit none
    real(psb_dpk_), intent(in)    :: vin
    real(psb_dpk_), intent(in) :: vinout
    integer(psb_ipk_) :: i
    real(psb_dpk_) :: w, z
    w = max( vin, vinout )
    z = min( vin, vinout )
    if ( z == dzero ) then
      caf_dnrm2 = w
    else
      caf_dnrm2 = w*sqrt( done +( z / w )**2 )
    end if
  end function caf_dnrm2

  pure real(psb_spk_) function caf_snrm2(vin, vinout)
    implicit none
    real(psb_spk_), intent(in)    :: vin
    real(psb_spk_), intent(in) :: vinout
    integer(psb_ipk_) :: i
    real(psb_spk_) :: w, z
    w = max( vin, vinout )
    z = min( vin, vinout )
    if ( z == dzero ) then
      caf_snrm2 = w
    else
      caf_snrm2 = w*sqrt( done +( z / w )**2 )
    end if
  end function caf_snrm2

  subroutine caf_camx_reduces(data,result_image)
    implicit none
    complex(psb_spk_) :: data
    complex(psb_spk_) :: data_buf, co_data
    complex(psb_spk_), save :: co_buffer[*]
    integer, optional :: result_image
    integer :: i, np, me
    real(psb_spk_) ::  z, w
    me=this_image()
    np = num_images()
    sync all    
    w = abs(data)
    co_buffer=data
    data_buf = co_buffer
    if (present(result_image)) then
      if (me == result_image) then
        do i=1,np
          co_data=co_buffer[i]
          z = abs(co_data)
          if (z > w) then  
            data_buf = co_data
            w = abs(data_buf)
            !print*,'i, data_buf',i,data_buf
          endif
        enddo
        data=data_buf
      endif
    else
      do i=1,np
        co_data=co_buffer[i]
        z = abs(co_data)
        if (z > w) then  
          data_buf = co_data
          w = abs(data_buf)
        endif
      enddo
      data=data_buf
    endif
    sync all
  end subroutine caf_camx_reduces

  subroutine caf_zamx_reduces(data,result_image)
    implicit none
    complex(psb_dpk_) :: data
    complex(psb_dpk_) :: data_buf, co_data
    integer, optional :: result_image
    complex(psb_dpk_), save :: co_buffer[*]
    integer :: i, np, me
    real(psb_dpk_) ::  z, w
    me=this_image()
    np = num_images()
    sync all    
    w = abs(data)
    co_buffer=data
    data_buf = co_buffer
    if (present(result_image)) then
      if (me == result_image) then
        do i=1,np
          co_data=co_buffer[i]
          z = abs(co_data)
          if (z > w) then  
            data_buf = co_data
            w = abs(data_buf)
            !print*,'i, data_buf',i,data_buf
          endif
        enddo
        data=data_buf
      endif
    else
      do i=1,np
        co_data=co_buffer[i]
        z = abs(co_data)
        if (z > w) then  
          data_buf = co_data
          w = abs(data_buf)
        endif
      enddo
      data=data_buf
    endif
    sync all
  end subroutine caf_zamx_reduces

  subroutine caf_camx_reducev(data,result_image)
    implicit none
    complex(psb_spk_) :: data(:)
    complex(psb_spk_), allocatable :: data_buf(:), co_data(:)
    complex(psb_dpk_), allocatable :: co_buffer(:)[:]
    integer, optional :: result_image
    integer :: i, np, me, size_, j
    real(psb_spk_), allocatable ::  z(:), w(:)
    sync all
    me=this_image()
    np = num_images()
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(data_buf)) deallocate(co_data)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    size_=size(data,1)
    allocate(data_buf(size_),co_data(size_), w(size_), z(size_))
    allocate(co_buffer(size_)[*])
    w = abs(data)
    data_buf = data
    co_buffer = data
    if (present(result_image)) then
      if (me == result_image) then
        do i=1,np
          co_data(1:size_)=co_buffer(1:size_)[i]
          z = abs(co_data)
          do j=1,size(co_data,1)
            if (z(j) > w(j)) then  
              data_buf(j) = co_data(j)
              w(j) = abs(data_buf(j))
            endif
          enddo
        enddo
        data=data_buf
      endif
    else
      do i=1,np
        co_data(1:size_)=co_buffer(1:size_)[i]
        z = abs(co_data)
        do j=1,size(co_data,1)
          if (z(j) > w(j)) then  
            data_buf(j) = co_data(j)
            w(j) = abs(data_buf(j))
          endif
	enddo
      enddo
      data=data_buf
    endif
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(data_buf)) deallocate(co_data)
    if (allocated(co_buffer)) deallocate(co_buffer)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    sync all
  end subroutine caf_camx_reducev

  subroutine caf_zamx_reducev(data,result_image)
    implicit none
    complex(psb_dpk_) :: data(:)
    complex(psb_dpk_), allocatable :: data_buf(:), co_data(:)
    complex(psb_dpk_), allocatable :: co_buffer(:)[:]
    integer, optional :: result_image
    integer :: i, np, me, size_, j
    real(psb_dpk_), allocatable ::  z(:), w(:)

    sync all
    me=this_image()
    np = num_images()
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(data_buf)) deallocate(co_data)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    size_=size(data,1)
    allocate(data_buf(size_),co_data(size_), w(size_), z(size_))
    allocate(co_buffer(size_)[*])
    w = abs(data)
    data_buf = data
    co_buffer = data
    if (present(result_image)) then
      if (me == result_image) then
        do i=1,np
          co_data(1:size_)=co_buffer(1:size_)[i]
          z = abs(co_data)
          do j=1,size(co_data,1)
            if (z(j) > w(j)) then  
              data_buf(j) = co_data(j)
              w(j) = abs(data_buf(j))
            endif
          enddo
        enddo
        data=data_buf
      endif
    else
      do i=1,np
        co_data(1:size_)=co_buffer(1:size_)[i]
        z = abs(co_data)
        do j=1,size(co_data,1)
          if (z(j) > w(j)) then  
            data_buf(j) = co_data(j)
            w(j) = abs(data_buf(j))
          endif
	enddo
      enddo
      data=data_buf
    endif
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(data_buf)) deallocate(co_data)
    if (allocated(co_buffer)) deallocate(co_buffer)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    sync all  
end subroutine caf_zamx_reducev

  subroutine caf_camx_reducem(data,result_image)
    implicit none
    complex(psb_spk_):: data(:,:)
    complex(psb_spk_), allocatable :: data_buf(:,:), co_data(:,:)
    complex(psb_dpk_), allocatable :: co_buffer(:,:)[:]
    integer, optional :: result_image
    integer :: i, np, me, size1, size2, j, k
    real(psb_spk_), allocatable ::  z(:,:), w(:,:)
    sync all
    me=this_image()
    np = num_images()
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(data_buf)) deallocate(co_data)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    if (allocated(co_buffer)) deallocate(co_buffer)
    size1=size(data,1)
    size2=size(data,2)
    allocate(data_buf(size1,size2),co_data(size1,size2), w(size1,size2), z(size1,size2))
    allocate(co_buffer(size1,size2)[*])
    w = abs(data)
    data_buf = data
    co_buffer = data
    if (present(result_image)) then
      if (me == result_image) then
        do i=1,np
          co_data(1:size1,1:size2)=co_buffer(1:size1,1:size2)[i]
          z = abs(co_data)
          do j=1,size1
            do k=1,size2
              if (z(j,k) > w(j,k)) then  
                data_buf(j,k) = co_data(j,k)
                w(j,k) = abs(data_buf(j,k))
              endif
            enddo
          enddo
        enddo
        data=data_buf
      endif
    else
      do i=1,np
        co_data(1:size1,1:size2)=co_buffer(1:size1,1:size2)[i]
        z = abs(co_data)
        do j=1,size1
          do k=1,size2
            if (z(j,k) > w(j,k)) then  
              data_buf(j,k) = co_data(j,k)
              w(j,k) = abs(data_buf(j,k))
            endif
	  enddo
        enddo
      enddo
      data=data_buf
    endif
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(data_buf)) deallocate(co_data)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    if (allocated(co_buffer)) deallocate(co_buffer)
    sync all
  end subroutine caf_camx_reducem

  subroutine caf_zamx_reducem(data,result_image)
    implicit none
    complex(psb_dpk_):: data(:,:)
    complex(psb_dpk_), allocatable :: data_buf(:,:), co_data(:,:)
    complex(psb_dpk_), allocatable:: co_buffer(:,:)[:]
    integer, optional :: result_image
    integer :: i, np, me, size1, size2, j, k
    real(psb_dpk_), allocatable ::  z(:,:), w(:,:)
    sync all
    me=this_image()
    np = num_images()
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(data_buf)) deallocate(co_data)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    if (allocated(co_buffer)) deallocate(co_buffer)
    size1=size(data,1)
    size2=size(data,2)
    allocate(data_buf(size1,size2),co_data(size1,size2), w(size1,size2), z(size1,size2))
    allocate(co_buffer(size1,size2)[*])
    w = abs(data)
    data_buf = data
    co_buffer = data
    if (present(result_image)) then
      if (me == result_image) then
        do i=1,np
          co_data(1:size1,1:size2)=co_buffer(1:size1,1:size2)[i]
          z = abs(co_data)
          do j=1,size1
            do k=1,size2
              if (z(j,k) > w(j,k)) then  
                data_buf(j,k) = co_data(j,k)
                w(j,k) = abs(data_buf(j,k))
              endif
            enddo
          enddo
        enddo
        data=data_buf
      endif
    else
      do i=1,np
        co_data(1:size1,1:size2)=co_buffer(1:size1,1:size2)[i]
        z = abs(co_data)
        do j=1,size1
          do k=1,size2
            if (z(j,k) > w(j,k)) then  
              data_buf(j,k) = co_data(j,k)
              w(j,k) = abs(data_buf(j,k))
            endif
	  enddo
        enddo
      enddo
      data=data_buf
    endif
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(data_buf)) deallocate(co_data)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    if (allocated(co_buffer)) deallocate(co_buffer)
    sync all
  end subroutine caf_zamx_reducem

  subroutine caf_camn_reduces(data,result_image)
    implicit none
    complex(psb_spk_) :: data
    complex(psb_dpk_), save:: co_buffer[*]
    complex(psb_spk_) :: data_buf, co_data
    integer, optional :: result_image
    integer :: i, np, me
    real(psb_spk_) ::  z, w
    me=this_image()
    np = num_images()
    sync all    
    w = abs(data)
    data_buf = data
    co_buffer = data
    if (present(result_image)) then
      if (me == result_image) then
        do i=1,np
          co_data=co_buffer[i]
          z = abs(co_data)
          if (z < w) then  
            data_buf = co_data
            w = abs(data_buf)
            !print*,'i, data_buf',i,data_buf
          endif
        enddo
        data=data_buf
      endif
    else
      do i=1,np
        co_data=co_buffer[i]
        z = abs(co_data)
        if (z < w) then  
          data_buf = co_data
          w = abs(data_buf)
        endif
      enddo
      data=data_buf
    endif
    sync all
  end subroutine caf_camn_reduces

  subroutine caf_zamn_reduces(data,result_image)
    implicit none
    complex(psb_dpk_):: data
    complex(psb_dpk_) :: data_buf, co_data
    complex(psb_dpk_), save :: co_buffer[*]
    integer, optional :: result_image
    integer :: i, np, me
    real(psb_dpk_) ::  z, w
    me=this_image()
    np = num_images()
    sync all    
    w = abs(data)
    data_buf = data
    co_buffer = data
    if (present(result_image)) then
      if (me == result_image) then
        do i=1,np
          co_data=co_buffer[i]
          z = abs(co_data)
          if (z < w) then  
            data_buf = co_data
            w = abs(data_buf)
            !print*,'i, data_buf',i,data_buf
          endif
        enddo
        data=data_buf
      endif
    else
      do i=1,np
        co_data=co_buffer[i]
        z = abs(co_data)
        if (z < w) then  
          data_buf = co_data
          w = abs(data_buf)
        endif
      enddo
      data=data_buf
    endif
    sync all
  end subroutine caf_zamn_reduces

  subroutine caf_camn_reducev(data,result_image)
    implicit none
    complex(psb_spk_) :: data(:)
    complex(psb_spk_), allocatable :: data_buf(:), co_data(:)
    complex(psb_spk_), allocatable :: co_buffer(:)[:]
    integer, optional :: result_image
    integer :: i, np, me, size_, j
    real(psb_spk_), allocatable ::  z(:), w(:)
    sync all
    me=this_image()
    np = num_images()
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(co_data)) deallocate(co_data)
    if (allocated(co_buffer)) deallocate(co_buffer)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    size_=size(data,1)
    allocate(data_buf(size_),co_data(size_), w(size_), z(size_))
    allocate(co_buffer(size_)[*])
    w = abs(data)
    data_buf = data
    co_buffer=data
    if (present(result_image)) then
      if (me == result_image) then
        do i=1,np
          co_data(1:size_)=co_buffer(1:size_)[i]
          z = abs(co_data)
          do j=1,size(co_data,1)
            if (z(j) < w(j)) then  
              data_buf(j) = co_data(j)
              w(j) = abs(data_buf(j))
            endif
          enddo
        enddo
        data=data_buf
      endif
    else
      do i=1,np
        co_data(1:size_)=co_buffer(1:size_)[i]
        z = abs(co_data)
        do j=1,size(co_data,1)
          if (z(j) < w(j)) then  
            data_buf(j) = co_data(j)
            w(j) = abs(data_buf(j))
          endif
	enddo
      enddo
      data=data_buf
    endif
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(data_buf)) deallocate(co_data)
    if (allocated(data_buf)) deallocate(co_buffer)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    sync all
  end subroutine caf_camn_reducev

  subroutine caf_zamn_reducev(data,result_image)
    implicit none
    complex(psb_dpk_) :: data(:)
    complex(psb_dpk_), allocatable :: co_buffer(:)[:]
    complex(psb_dpk_), allocatable :: data_buf(:), co_data(:)
    integer, optional :: result_image
    integer :: i, np, me, size_, j
    real(psb_dpk_), allocatable ::  z(:), w(:)
    sync all
    me=this_image()
    np = num_images()
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(co_data)) deallocate(co_data)
    if (allocated(co_buffer)) deallocate(co_buffer)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    size_=size(data,1)
    allocate(data_buf(size_),co_data(size_), w(size_), z(size_))
    allocate(co_buffer(size_)[*])
    w = abs(data)
    data_buf = data
    co_buffer=data
    if (present(result_image)) then
      if (me == result_image) then
        do i=1,np
          co_data(1:size_)=co_buffer(1:size_)[i]
          z = abs(co_data)
          do j=1,size(co_data,1)
            if (z(j) < w(j)) then  
              data_buf(j) = co_data(j)
              w(j) = abs(data_buf(j))
            endif
          enddo
        enddo
        data=data_buf
      endif
    else
      do i=1,np
        co_data(1:size_)=co_buffer(1:size_)[i]
        z = abs(co_data)
        do j=1,size(co_data,1)
          if (z(j) < w(j)) then  
            data_buf(j) = co_data(j)
            w(j) = abs(data_buf(j))
          endif
	enddo
      enddo
      data=data_buf
    endif
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(data_buf)) deallocate(co_data)
    if (allocated(data_buf)) deallocate(co_buffer)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    sync all
  end subroutine caf_zamn_reducev

  subroutine caf_camn_reducem(data,result_image)
    implicit none
    complex(psb_spk_) :: data(:,:)
    complex(psb_spk_), allocatable :: co_buffer(:,:)[:]
    complex(psb_spk_), allocatable :: data_buf(:,:), co_data(:,:)
    integer, optional :: result_image
    integer :: i, np, me, size1, size2, j, k
    real(psb_spk_), allocatable ::  z(:,:), w(:,:)
    sync all
    me=this_image()
    np = num_images()
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(co_data)) deallocate(co_data)
    if (allocated(co_buffer)) deallocate(co_buffer)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    size1=size(data,1)
    size2=size(data,2)
    allocate(data_buf(size1,size2),co_data(size1,size2), w(size1,size2), z(size1,size2))
    allocate(co_buffer(size1,size2)[*])
    w = abs(data)
    data_buf = data
    co_buffer=data
    if (present(result_image)) then
      if (me == result_image) then
        do i=1,np
          co_data(1:size1,1:size2)=co_buffer(1:size1,1:size2)[i]
          z = abs(co_data)
          do j=1,size1
            do k=1,size2
              if (z(j,k) < w(j,k)) then  
                data_buf(j,k) = co_data(j,k)
                w(j,k) = abs(data_buf(j,k))
              endif
            enddo
          enddo
        enddo
        data=data_buf
      endif
    else
      do i=1,np
        co_data(1:size1,1:size2)=co_buffer(1:size1,1:size2)[i]
        z = abs(co_data)
        do j=1,size1
          do k=1,size2
            if (z(j,k) < w(j,k)) then  
              data_buf(j,k) = co_data(j,k)
              w(j,k) = abs(data_buf(j,k))
            endif
	  enddo
        enddo
      enddo
      data=data_buf
    endif
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(co_data)) deallocate(co_data)
    if (allocated(co_buffer)) deallocate(co_buffer)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    sync all
  end subroutine caf_camn_reducem

  subroutine caf_zamn_reducem(data,result_image)
    implicit none
    complex(psb_dpk_) :: data(:,:)
    complex(psb_dpk_), allocatable :: co_buffer(:,:)[:]
    complex(psb_dpk_), allocatable :: data_buf(:,:), co_data(:,:)
    integer, optional :: result_image
    integer :: i, np, me, size1, size2, j, k
    real(psb_dpk_), allocatable ::  z(:,:), w(:,:)
    sync all
    me=this_image()
    np = num_images()
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(co_data)) deallocate(co_data)
    if (allocated(co_buffer)) deallocate(co_buffer)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    size1=size(data,1)
    size2=size(data,2)
    allocate(data_buf(size1,size2),co_data(size1,size2), w(size1,size2), z(size1,size2))
    allocate(co_buffer(size1,size2)[*])
    w = abs(data)
    data_buf = data
    co_buffer=data
    if (present(result_image)) then
      if (me == result_image) then
        do i=1,np
          co_data(1:size1,1:size2)=co_buffer(1:size1,1:size2)[i]
          z = abs(co_data)
          do j=1,size1
            do k=1,size2
              if (z(j,k) < w(j,k)) then  
                data_buf(j,k) = co_data(j,k)
                w(j,k) = abs(data_buf(j,k))
              endif
            enddo
          enddo
        enddo
        data=data_buf
      endif
    else
      do i=1,np
        co_data(1:size1,1:size2)=co_buffer(1:size1,1:size2)[i]
        z = abs(co_data)
        do j=1,size1
          do k=1,size2
            if (z(j,k) < w(j,k)) then  
              data_buf(j,k) = co_data(j,k)
              w(j,k) = abs(data_buf(j,k))
            endif
	  enddo
        enddo
      enddo
      data=data_buf
    endif
    if (allocated(data_buf)) deallocate(data_buf)
    if (allocated(co_data)) deallocate(co_data)
    if (allocated(co_buffer)) deallocate(co_buffer)
    if (allocated(w)) deallocate(w)
    if (allocated(z)) deallocate(z)
    sync all
  end subroutine caf_zamn_reducem
end module psb_caf_mod

