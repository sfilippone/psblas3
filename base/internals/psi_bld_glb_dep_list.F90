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
subroutine psi_i_bld_glb_dep_list(ictxt,loc_dl,length_dl,dep_list,dl_lda,info)
  use psi_mod, psb_protect_name => psi_i_bld_glb_dep_list
#ifdef MPI_MOD
  use mpi
#endif
  use psb_penv_mod
  use psb_const_mod
  use psb_error_mod
  use psb_desc_mod
  use psb_sort_mod
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif
  !     ....scalar parameters...
  integer(psb_ipk_), intent(in)  :: ictxt
  integer(psb_ipk_), intent(out) :: dl_lda
  integer(psb_ipk_), intent(in)  :: loc_dl(:), length_dl(0:)
  integer(psb_ipk_), allocatable, intent(out) :: dep_list(:,:)
  integer(psb_ipk_), intent(out) :: info


  !     .....local arrays....
  integer(psb_ipk_) :: int_err(5)

  !     .....local scalars...
  integer(psb_ipk_) :: i, proc,j,err_act
  integer(psb_ipk_) :: err
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_mpk_) :: iictxt, icomm, me, np, minfo
  logical, parameter :: dist_symm_list=.false., print_dl=.false.
  character  name*20
  name='psi_bld_glb_dep_list'

  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  iictxt = ictxt 
  info = psb_success_

  call psb_info(iictxt,me,np)


  dl_lda = length_dl(me)
  call psb_max(iictxt, dl_lda)

  if (debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),': Dep_list length ',length_dl(me),dl_lda
  dl_lda = max(dl_lda,1) 
  allocate(dep_list(dl_lda,0:np),stat=info)
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if
  icomm = psb_get_mpi_comm(iictxt)
  call mpi_allgather(loc_dl,dl_lda,psb_mpi_ipk_,&
       & dep_list,dl_lda,psb_mpi_ipk_,icomm,minfo)

  info = minfo  
  if (info /= psb_success_) then 
    info=psb_err_internal_error_
    goto 9999
  endif
  if (print_dl) then
    if (me == 0) then
      write(0,*) ' Dep_list '
      do i=0,np-1
        j = length_dl(i) 
        write(0,*) 'Proc ',i,':',dep_list(1:j,i)
      end do
      flush(0)
    end if
    call psb_barrier(ictxt)
  end if

  call psb_erractionrestore(err_act)
  return


9999 continue

  call psb_errpush(info,name,i_err=int_err)
  call psb_error_handler(err_act)

  return

end subroutine psi_i_bld_glb_dep_list

subroutine psi_i_bld_glb_csr_dep_list(ictxt,loc_dl,length_dl,c_dep_list,dl_ptr,info)
  use psi_mod, psb_protect_name => psi_i_bld_glb_csr_dep_list
#ifdef MPI_MOD
  use mpi
#endif
  use psb_penv_mod
  use psb_const_mod
  use psb_error_mod
  use psb_desc_mod
  use psb_sort_mod
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif
  !     ....scalar parameters...
  integer(psb_ipk_), intent(in)  :: ictxt
  integer(psb_ipk_), intent(in)  :: loc_dl(:), length_dl(0:)
  integer(psb_ipk_), allocatable, intent(out) :: c_dep_list(:), dl_ptr(:) 
  integer(psb_ipk_), intent(out) :: info


  !     .....local arrays....
  integer(psb_ipk_) :: int_err(5)

  !     .....local scalars...
  integer(psb_ipk_) :: i, proc,j,err_act, length, myld
  integer(psb_ipk_) :: err
  integer(psb_ipk_) :: debug_level, debug_unit
  integer(psb_mpk_) :: iictxt, icomm, me, np, minfo
  logical, parameter :: dist_symm_list=.false., print_dl=.true. 
  character  name*20
  name='psi_bld_glb_csr_dep_list'

  call psb_erractionsave(err_act)
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  iictxt = ictxt 
  info = psb_success_

  call psb_info(iictxt,me,np)

  myld = length_dl(me)
  length = sum(length_dl(0:np-1)) 
  allocate(dl_ptr(0:np),stat=info)
  dl_ptr(0) = 0
  do i=1, np
    dl_ptr(i) = dl_ptr(i-1) + length_dl(i-1)
  end do

  if (length /= dl_ptr(np)) then
    write(0,*) me,trim(name),' Inconsistency: ',length,dl_ptr(np)
  end if
  if (debug_level >= psb_debug_inner_) &
       & write(debug_unit,*) me,' ',trim(name),': Dep_list length ',length_dl(me)

  allocate(c_dep_list(length),stat=info)
  
  if (info /= psb_success_) then 
    call psb_errpush(psb_err_from_subroutine_,name,a_err='Allocate')
    goto 9999      
  end if
  icomm = psb_get_mpi_comm(iictxt)
  call mpi_allgatherv(loc_dl,myld,psb_mpi_ipk_,&
       & c_dep_list,length_dl,dl_ptr,psb_mpi_ipk_,icomm,minfo)
  
  info = minfo  
  if (info /= psb_success_) then 
    info=psb_err_internal_error_
    goto 9999
  endif
  dl_ptr = dl_ptr + 1
  
  if (print_dl) then
    if (me == 0) then
      block
        character(len=80) :: fname, frmt 
        integer :: ni, nl,ldl, lname, iout = 87
        ldl = dl_ptr(np) 
        ni  = floor(log10(1.0*np)) + 1
        nl  = floor(log10(1.0*ldl)) + 1
        write(frmt,'(a,i3.3,a,i3.3,a)') '(a,i',ni,'.',ni,')'
        write(fname,frmt) 'dep_list_p_',np
        lname = len_trim(fname) 
        write(frmt,'(a,i3.3,a,i3.3,a)') '(a,i',nl,'.',nl,')'
        write(fname(lname+1:lname+nl+3),frmt) '_l_',ldl
        fname = trim(fname)//'.mtx'
        open(iout,file=fname)
        if (.true.) then
          write(iout,*)   np, np, ldl-1
          do i=0,np-1
            do j=dl_ptr(i),dl_ptr(i+1)-1
              write(iout,*) i+1,c_dep_list(j)+1
            end do
          end do
        else
          write(iout,*) ' Dep_list '
          do i=0,np-1
            write(iout,*) 'Proc ',i,':',c_dep_list(dl_ptr(i):dl_ptr(i+1)-1) 
          end do
        end if
        close(iout)
      end block
    end if
    call psb_barrier(ictxt)
  end if

  call psb_erractionrestore(err_act)
  return


9999 continue

  call psb_errpush(info,name,i_err=int_err)
  call psb_error_handler(err_act)

  return

end subroutine psi_i_bld_glb_csr_dep_list
