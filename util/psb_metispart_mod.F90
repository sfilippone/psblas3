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
  use psb_base_mod, only : psb_ipk_, psb_lpk_, psb_mpk_, psb_epk_, &
       & psb_err_unit, psb_spk_, psb_dpk_,&
       & psb_lsspmat_type, psb_lcspmat_type,&
       & psb_ldspmat_type, psb_lzspmat_type,  &
       & psb_ls_csr_sparse_mat, psb_ld_csr_sparse_mat, &
       & psb_lc_csr_sparse_mat, psb_lz_csr_sparse_mat
  public part_graph, build_mtpart, distr_mtpart,&
       & getv_mtpart, free_part
  private 
  integer(psb_lpk_), allocatable, save :: graph_vect(:)

  interface build_mtpart
    module procedure ld_mat_build_mtpart, ls_mat_build_mtpart,&
         & lz_mat_build_mtpart, lc_mat_build_mtpart
  end interface

  interface psi_build_mtpart
    subroutine psi_l_build_mtpart(n,ja,irp,nparts,vect, weights)
      import :: psb_lpk_, psb_spk_, psb_dpk_
      implicit none 
      integer(psb_lpk_), intent(in) :: n, nparts
      integer(psb_lpk_), intent(in) :: ja(:), irp(:)
      integer(psb_lpk_), allocatable, intent(inout) :: vect(:)
#if defined(METIS_REAL_32) || !defined(HAVE_METIS)
      real(psb_spk_),optional, intent(in) :: weights(:)
#elif defined(METIS_REAL_64)
      real(psb_dpk_),optional, intent(in) :: weights(:)
#else
      choke on me;
#endif
    end subroutine psi_l_build_mtpart
  end interface
  
contains
  
  subroutine part_graph(global_indx,n,np,pv,nv)
    implicit none
    integer(psb_lpk_), intent(in)  :: global_indx, n
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
    integer(psb_ipk_) :: me, np, info
    integer(psb_lpk_) :: n

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


  subroutine ld_mat_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_ldspmat_type), intent(in) :: a
    integer(psb_lpk_) :: nparts
    real(psb_dpk_), optional :: weights(:)
    

    select type (aa=>a%a) 
    type is (psb_ld_csr_sparse_mat)
      call ld_csr_build_mtpart(aa,nparts,weights)        
    class default
      write(psb_err_unit,*) 'Sorry, right now we only take CSR input!'
    end select

  end subroutine ld_mat_build_mtpart
  
  subroutine ld_csr_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_ld_csr_sparse_mat), intent(in) :: a
    integer(psb_lpk_) :: nparts
    real(psb_dpk_), optional :: weights(:)
#if defined(METIS_REAL_32) || !defined(HAVE_METIS)
    real(psb_spk_), allocatable :: wgh_(:)
#elif defined(METIS_REAL_64)
    real(psb_dpk_), allocatable :: wgh_(:)
#else
      choke on me;
#endif

    if (present(weights)) then 
      if (size(weights)==nparts) then 
        wgh_ = weights
      end if
    end if
    if (allocated(wgh_)) then 
      call psi_build_mtpart(a%get_nrows(),a%ja,a%irp,nparts,graph_vect,wgh_)
    else
      call psi_build_mtpart(a%get_nrows(),a%ja,a%irp,nparts,graph_vect)
    end if

  end subroutine ld_csr_build_mtpart
  
  subroutine lz_mat_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_lzspmat_type), intent(in) :: a
    integer(psb_lpk_) :: nparts
    real(psb_dpk_), optional :: weights(:)
    

    select type (aa=>a%a) 
    type is (psb_lz_csr_sparse_mat) 
      call lz_csr_build_mtpart(aa,nparts,weights)        
    class default
      write(psb_err_unit,*) 'Sorry, right now we only take CSR input!'
    end select

  end subroutine lz_mat_build_mtpart

  subroutine lz_csr_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_lz_csr_sparse_mat), intent(in) :: a
    integer(psb_lpk_) :: nparts
    real(psb_dpk_), optional :: weights(:)
#if defined(METIS_REAL_32) || !defined(HAVE_METIS)
    real(psb_spk_), allocatable :: wgh_(:)
#elif defined(METIS_REAL_64)
    real(psb_dpk_), allocatable :: wgh_(:)
#else
      choke on me;
#endif

    if (present(weights)) then 
      if (size(weights)==nparts) then 
        wgh_ = weights
      end if
    end if
    if (allocated(wgh_)) then 
      call psi_build_mtpart(a%get_nrows(),a%ja,a%irp,nparts,graph_vect,wgh_)
    else
      call psi_build_mtpart(a%get_nrows(),a%ja,a%irp,nparts,graph_vect)
    end if

  end subroutine lz_csr_build_mtpart
  
  subroutine ls_mat_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_lsspmat_type), intent(in) :: a
    integer(psb_lpk_) :: nparts
    real(psb_spk_), optional :: weights(:)
    

    select type (aa=>a%a) 
    type is (psb_ls_csr_sparse_mat)
      call ls_csr_build_mtpart(aa,nparts,weights)        
    class default
      write(psb_err_unit,*) 'Sorry, right now we only take CSR input!'
    end select

  end subroutine ls_mat_build_mtpart

  subroutine lc_mat_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_lcspmat_type), intent(in) :: a
    integer(psb_lpk_) :: nparts
    real(psb_spk_), optional :: weights(:)
    

    select type (aa=>a%a) 
    type is (psb_lc_csr_sparse_mat)
      call lc_csr_build_mtpart(aa,nparts,weights)        
    class default
      write(psb_err_unit,*) 'Sorry, right now we only take CSR input!'
    end select

  end subroutine lc_mat_build_mtpart

  subroutine lc_csr_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_lc_csr_sparse_mat), intent(in) :: a
    integer(psb_lpk_) :: nparts
    real(psb_spk_), optional :: weights(:)
#if defined(METIS_REAL_32) || !defined(HAVE_METIS)
    real(psb_spk_), allocatable :: wgh_(:)
#elif defined(METIS_REAL_64)
    real(psb_dpk_), allocatable :: wgh_(:)
#else
      choke on me;
#endif

  
    if (present(weights)) then 
      if (size(weights)==nparts) then 
        wgh_ = weights
      end if
    end if
    if (allocated(wgh_)) then 
      call psi_build_mtpart(a%get_nrows(),a%ja,a%irp,nparts,graph_vect,wgh_)
    else
      call psi_build_mtpart(a%get_nrows(),a%ja,a%irp,nparts,graph_vect)
    end if

  end subroutine lc_csr_build_mtpart
  
  subroutine ls_csr_build_mtpart(a,nparts,weights)
    use psb_base_mod
    implicit none 
    type(psb_ls_csr_sparse_mat), intent(in) :: a
    integer(psb_lpk_) :: nparts
    real(psb_spk_), optional :: weights(:)
#if defined(METIS_REAL_32) || !defined(HAVE_METIS)
    real(psb_spk_), allocatable :: wgh_(:)
#elif defined(METIS_REAL_64)
    real(psb_dpk_), allocatable :: wgh_(:)
#else
      choke on me;
#endif

  
    if (present(weights)) then 
      if (size(weights)==nparts) then 
        wgh_ = weights
      end if
    end if
    if (allocated(wgh_)) then 
      call psi_build_mtpart(a%get_nrows(),a%ja,a%irp,nparts,graph_vect,wgh_)
    else
      call psi_build_mtpart(a%get_nrows(),a%ja,a%irp,nparts,graph_vect)
    end if
    
  end subroutine ls_csr_build_mtpart
  
  !
  ! WARNING: called IRET otherwise Intel compiler complains,
  ! methinks it's a compiler bug, will need to report. 
  !
  subroutine free_part(iret)
    implicit none 
    integer(psb_ipk_), intent(out) :: iret

    deallocate(graph_vect,stat=iret)
    return
  end subroutine free_part    

end module psb_metispart_mod

