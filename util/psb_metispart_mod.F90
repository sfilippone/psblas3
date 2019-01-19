!   
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone    
!        Alfredo Buttari      
!   
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!   
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!   
!    
!
! Purpose: 
!  Provide a set of subroutines to define a data distribution based on 
!  a graph partitioning routine from METIS. May serve as the basis
!  for interfacing other graph partitioning tools.
! 
!  Subroutines:
!  
!  BUILD_MTPART(A,NPARTS): This subroutine will be called by the root
!    process to build define the data distribuition mapping. 
!      Input parameters:
!        TYPE(D_SPMAT) :: A   The input matrix. The coefficients are
!                             ignored; only the structure is used.
!        integer(psb_ipk_) :: NPARTS  How many parts we are requiring to the 
!                                 partition utility
! 
!  DISTR_MTPART(ROOT,ICTXT): This subroutine will be called by
!      all processes to distribute the information computed by the root
!      process, to be used subsequently.
!
!
!  PART_GRAPH : The subroutine to be passed to PSBLAS sparse library;
!      uses information prepared by the previous two subroutines.
!
module psb_metispart_mod
  use psb_base_mod, only : psb_ipk_, psb_sspmat_type, psb_cspmat_type,&
       & psb_dspmat_type, psb_zspmat_type, psb_err_unit, psb_mpik_,&
       & psb_s_csr_sparse_mat, psb_d_csr_sparse_mat, &
       & psb_c_csr_sparse_mat, psb_z_csr_sparse_mat
  public part_graph, build_mtpart, distr_mtpart,&
       & getv_mtpart, free_part
  private 
  integer(psb_ipk_), allocatable, save :: graph_vect(:)

  interface build_mtpart
    module procedure build_mtpart,&
         & d_mat_build_mtpart, s_mat_build_mtpart,&
         & z_mat_build_mtpart, c_mat_build_mtpart, &
         & d_csr_build_mtpart, s_csr_build_mtpart,&
         & z_csr_build_mtpart, c_csr_build_mtpart

  end interface

contains
  
  subroutine part_graph(global_indx,n,np,pv,nv)
    implicit none 
    integer(psb_ipk_), intent(in)  :: global_indx, n
    integer(psb_ipk_), intent(in)  :: np
    integer(psb_ipk_), intent(out) :: nv
    integer(psb_ipk_), intent(out) :: pv(*)
    
    IF (.not.allocated(graph_vect)) then
       write(psb_err_unit,*) 'Fatal error in PART_GRAPH: vector GRAPH_VECT ',&
	    & 'not initialized'
       return
    endif
    if ((global_indx<1).or.(global_indx > size(graph_vect))) then       
       write(psb_err_unit,*) 'Fatal error in PART_GRAPH: index GLOBAL_INDX ',&
	    & 'outside GRAPH_VECT bounds',global_indx,size(graph_vect)
       return
    endif
    nv = 1
    pv(nv) = graph_vect(global_indx)
    return
  end subroutine part_graph


  subroutine distr_mtpart(root, ictxt)
    use psb_base_mod
    implicit none 
    integer(psb_ipk_) :: root, ictxt
    integer(psb_ipk_) :: n, me, np, info
    
    call psb_info(ictxt,me,np)

    if (.not.((root>=0).and.(root<np))) then 
      write(psb_err_unit,*) 'Fatal error in DISTR_MTPART: invalid ROOT  ',&
           & 'coordinates '
      call psb_abort(ictxt)
      return
    endif

    if (me == root) then 
      if (.not.allocated(graph_vect)) then
        write(psb_err_unit,*) 'Fatal error in DISTR_MTPART: vector GRAPH_VECT ',&
             & 'not initialized'
        call psb_abort(ictxt)
        return
      endif
      n = size(graph_vect)
      call psb_bcast(ictxt,n,root=root)
    else 
      call psb_bcast(ictxt,n,root=root)

      allocate(graph_vect(n),stat=info)
      if (info /= psb_success_) then
        write(psb_err_unit,*) 'Fatal error in DISTR_MTPART: memory allocation ',&
             & ' failure.'
        return
      endif
    endif
    call psb_bcast(ictxt,graph_vect(1:n),root=root)

    return

  end subroutine distr_mtpart
  
  subroutine  getv_mtpart(ivg)
    implicit none 
    integer(psb_ipk_), allocatable, intent(out)  :: ivg(:)
    if (allocated(graph_vect)) then 
      allocate(ivg(size(graph_vect)))
      ivg(:) = graph_vect(:)
    end if
  end subroutine getv_mtpart
  
  subroutine d_mat_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_dspmat_type), intent(in) :: a
    integer(psb_ipk_) :: nparts
    real(psb_dpk_), optional :: weights(:)
    real(psb_spk_), allocatable :: wgh_(:)
    

    select type (aa=>a%a) 
    type is (psb_d_csr_sparse_mat)
      if (present(weights)) then 
        if (size(weights)==nparts) then 
          wgh_ = weights
        end if
      end if
      if (allocated(wgh_)) then 
        call build_mtpart(aa%get_nrows(),aa%get_fmt(),aa%ja,aa%irp,nparts,wgh_)
      else
        call build_mtpart(aa%get_nrows(),aa%get_fmt(),aa%ja,aa%irp,nparts)
      end if
    class default
      write(psb_err_unit,*) 'Sorry, right now we only take CSR input!'
    end select

  end subroutine d_mat_build_mtpart

  
  subroutine z_mat_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_zspmat_type), intent(in) :: a
    integer(psb_ipk_) :: nparts
    real(psb_dpk_), optional :: weights(:)
    real(psb_spk_), allocatable :: wgh_(:)
    

    select type (aa=>a%a) 
    type is (psb_z_csr_sparse_mat) 
      if (present(weights)) then 
        if (size(weights)==nparts) then 
          wgh_ = weights
        end if
      end if
      if (allocated(wgh_)) then 
        call build_mtpart(aa%get_nrows(),aa%get_fmt(),aa%ja,aa%irp,nparts,wgh_)
      else
        call build_mtpart(aa%get_nrows(),aa%get_fmt(),aa%ja,aa%irp,nparts)
      end if
    class default
      write(psb_err_unit,*) 'Sorry, right now we only take CSR input!'
    end select

  end subroutine z_mat_build_mtpart

  
  subroutine s_mat_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_sspmat_type), intent(in) :: a
    integer(psb_ipk_) :: nparts
    real(psb_spk_), optional :: weights(:)
    

    select type (aa=>a%a) 
    type is (psb_s_csr_sparse_mat)
      call build_mtpart(aa%get_nrows(),aa%get_fmt(),aa%ja,aa%irp,nparts,weights)
    class default
      write(psb_err_unit,*) 'Sorry, right now we only take CSR input!'
    end select

  end subroutine s_mat_build_mtpart

  
  subroutine c_mat_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_cspmat_type), intent(in) :: a
    integer(psb_ipk_) :: nparts
    real(psb_spk_), optional :: weights(:)
    

    select type (aa=>a%a) 
    type is (psb_c_csr_sparse_mat)
      call build_mtpart(aa%get_nrows(),aa%get_fmt(),aa%ja,aa%irp,nparts,weights)
    class default
      write(psb_err_unit,*) 'Sorry, right now we only take CSR input!'
    end select

  end subroutine c_mat_build_mtpart

  
  subroutine d_csr_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_d_csr_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: nparts
    real(psb_dpk_), optional :: weights(:)
    real(psb_spk_), allocatable :: wgh_(:)


    if (present(weights)) then 
      if (size(weights)==nparts) then 
        wgh_ = weights
      end if
    end if
    if (allocated(wgh_)) then 
      call build_mtpart(a%get_nrows(),a%get_fmt(),a%ja,a%irp,nparts,wgh_)
    else
      call build_mtpart(a%get_nrows(),a%get_fmt(),a%ja,a%irp,nparts)
    end if

  end subroutine d_csr_build_mtpart

  subroutine z_csr_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_z_csr_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: nparts
    real(psb_dpk_), optional :: weights(:)
    real(psb_spk_), allocatable :: wgh_(:)
    

    if (present(weights)) then 
      if (size(weights)==nparts) then 
        wgh_ = weights
      end if
    end if
    if (allocated(wgh_)) then 
      call build_mtpart(a%get_nrows(),a%get_fmt(),a%ja,a%irp,nparts,wgh_)
    else
      call build_mtpart(a%get_nrows(),a%get_fmt(),a%ja,a%irp,nparts)
    end if

  end subroutine z_csr_build_mtpart

  
  subroutine s_csr_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_s_csr_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: nparts
    real(psb_spk_), optional :: weights(:)
    

    call build_mtpart(a%get_nrows(),a%get_fmt(),a%ja,a%irp,nparts,weights)

  end subroutine s_csr_build_mtpart

  
  subroutine c_csr_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_c_csr_sparse_mat), intent(in) :: a
    integer(psb_ipk_) :: nparts
    real(psb_spk_), optional :: weights(:)
    
    
    call build_mtpart(a%get_nrows(),a%get_fmt(),a%ja,a%irp,nparts,weights)

  end subroutine c_csr_build_mtpart



  subroutine build_mtpart(n,fida,ja,irp,nparts,weights)
    use psb_base_mod
    implicit none 
    integer(psb_ipk_) :: nparts
    integer(psb_ipk_) :: ja(:), irp(:)
    integer(psb_ipk_) :: n, i,numflag,nedc,wgflag
    character(len=5)     :: fida
    integer(psb_ipk_), parameter :: nb=512
    real(psb_dpk_), parameter :: seed=12345.d0
    integer(psb_ipk_) :: iopt(10),idummy(2),jdummy(2), info
    real(psb_spk_),optional :: weights(:)
    integer(psb_ipk_) :: nl,nptl
    integer(psb_ipk_), allocatable :: irpl(:),jal(:),gvl(:)
    real(psb_spk_),allocatable  :: wgh_(:)

#if defined(HAVE_METIS)
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
      if (psb_toupper(fida) == 'CSR') then 
        iopt(1) = 0
        numflag  = 1
        wgflag   = 0

!!$        write(*,*) 'Before allocation',nparts

        irpl = irp
        jal  = ja
        nl   = n
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
        write(psb_err_unit,*) 'Fatal error in BUILD_MTPART: matrix format ',&
             & ' failure. ', FIDA
        return
      endif
    else
      do i=1, n
        graph_vect(i) = 0
      enddo
    endif
#else
    write(psb_err_unit,*) 'Warning: METIS was not configured at PSBLAS compile time !'
#endif

    return

  end subroutine build_mtpart


  subroutine free_part(info)
    implicit none 
    integer(psb_ipk_) :: info
    
    deallocate(graph_vect,stat=info)
    return
  end subroutine free_part    

end module psb_metispart_mod

