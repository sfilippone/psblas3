  !
  ! Subroutine: psb_cdfree
  !   Frees a descriptor data structure.
  ! 
  ! Arguments: 
  !    desc_a   - type(psb_desc_type).         The communication descriptor to be freed.
subroutine psb_cd_destroy(desc)
  !...free descriptor structure...
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_desc_mod, psb_protect_name => psb_cd_destroy
#ifdef MPI_MOD
  use mpi
#endif
  Implicit None
#ifdef MPI_H
  include 'mpif.h'
#endif
  !....parameters...
  class(psb_desc_type), intent(inout) :: desc
  !...locals....
  integer(psb_ipk_) :: info, i, j


  if (allocated(desc%halo_index)) &
       &  deallocate(desc%halo_index,stat=info)

  if (allocated(desc%bnd_elem)) &
       &    deallocate(desc%bnd_elem,stat=info)

  if (allocated(desc%ovrlap_index)) &
       & deallocate(desc%ovrlap_index,stat=info)

  if (allocated(desc%ovrlap_elem)) &
       & deallocate(desc%ovrlap_elem,stat=info)
  if (allocated(desc%ovr_mst_idx)) &
       & deallocate(desc%ovr_mst_idx,stat=info)

  if (allocated(desc%lprm)) &
       & deallocate(desc%lprm,stat=info)
  if (allocated(desc%idx_space)) &
       & deallocate(desc%idx_space,stat=info)

  if (allocated(desc%sendtypes)) then 
    do j=1, size(desc%sendtypes,2) 
      do i=1, size(desc%sendtypes,1) 
        if (desc%sendtypes(i,j) /= mpi_datatype_null) then 
          call mpi_type_free(desc%sendtypes(i,j),info)
        end if
      end do
    end do
    deallocate(desc%sendtypes,stat=info)
  end if
  

  if (allocated(desc%recvtypes)) then 
    do j=1, size(desc%recvtypes,2) 
      do i=1, size(desc%recvtypes,1) 
        if (desc%recvtypes(i,j) /= mpi_datatype_null) then 
          call mpi_type_free(desc%recvtypes(i,j),info)
        end if
      end do
    end do
    deallocate(desc%recvtypes,stat=info)
  end if
  

  if (allocated(desc%indxmap)) then 
    call desc%indxmap%free()
    deallocate(desc%indxmap, stat=info)
  end if

  call desc%nullify()

  return

end subroutine psb_cd_destroy
