!!$ 
!!$              Parallel Sparse BLAS  v2.0
!!$    (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        University of Rome Tor Vergata
!!$ 
!!$  Redistribution and use in source and binary forms, with or without
!!$  modification, are permitted provided that the following conditions
!!$  are met:
!!$    1. Redistributions of source code must retain the above copyright
!!$       notice, this list of conditions and the following disclaimer.
!!$    2. Redistributions in binary form must reproduce the above copyright
!!$       notice, this list of conditions, and the following disclaimer in the
!!$       documentation and/or other materials provided with the distribution.
!!$    3. The name of the PSBLAS group or the names of its contributors may
!!$       not be used to endorse or promote products derived from this
!!$       software without specific written permission.
!!$ 
!!$  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!!$  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!!$  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!!$  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!!$  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!!$  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!!$  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!!$  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!!$  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!!$  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!!$  POSSIBILITY OF SUCH DAMAGE.
!!$ 
!!$  
! File: psb_cdren.f90
!
! Subroutine: psb_cdren
!    Updates a communication descriptor according to a renumbering scheme.
! 
! Parameters: 
!    trans    - character.                     Whether iperm or its transpose should be applied.
!    iperm    - integer,dimension(:).          The renumbering scheme.
!    desc_a   - type(<psb_desc_type>).         The communication descriptor to be updated.
!    info     - integer.                       Eventually returns an error code.
!
subroutine psb_cdren(trans,iperm,desc_a,info)
  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  use psb_penv_mod
  use psb_string_mod
  implicit none

  interface isaperm
    logical function isaperm(n,ip)
      integer, intent(in)    :: n   
      integer, intent(inout) :: ip(*)
    end function isaperm
  end interface

  !...parameters....
  type(psb_desc_type), intent(inout)    ::  desc_a
  integer, intent(inout)                ::  iperm(:)
  character, intent(in)                 :: trans
  integer, intent(out)                  :: info
  !....locals....
  integer                       :: i,j,np,me, n_col, kh, nh
  integer                       :: dectype
  integer                       :: ictxt,n_row, int_err(5), err_act
  real(kind(1.d0))              :: time(10), mpi_wtime, real_err(6)
  external mpi_wtime
  logical, parameter            :: debug=.false.
  character(len=20)             :: name

  if(psb_get_errstatus() /= 0) return 
  info=0
  call psb_erractionsave(err_act)
  name = 'psb_dcren'

  time(1) = mpi_wtime()

  ictxt   = psb_cd_get_context(desc_a)
  dectype = psb_cd_get_dectype(desc_a)
  n_row   = psb_cd_get_local_rows(desc_a)
  n_col   = psb_cd_get_local_cols(desc_a)

  ! check on blacs grid 
  call psb_info(ictxt, me, np)
  if (np == -1) then
    info = 2010
    call psb_errpush(info,name)
    goto 9999
  endif

  if (.not.psb_is_asb_desc(desc_a)) then 
    info = 600
    int_err(1) = dectype
    call psb_errpush(info,name,int_err)
    goto 9999
  endif

  if (iperm(1) /= 0) then 
    if (.not.isaperm(n_row,iperm)) then
      info = 610
      int_err(1) = iperm(1)
      call psb_errpush(info,name,int_err)
      goto 9999
    endif
  endif

  if (debug) write (*, *) '   begin matrix assembly...'

  !check on errors encountered in psdspins

  if ((iperm(1) /= 0))   then 

    if (debug) write(0,*) 'spasb: here we go with ',iperm(1) 
    deallocate(desc_a%lprm)
    allocate(desc_a%lprm(n_col))
    if (toupper(trans) == 'N') then 
      do i=1, n_row
        desc_a%lprm(iperm(i)) = i
      enddo
      do i=n_row+1,n_col
        desc_a%lprm(i) = i
      enddo
    else if (toupper(trans) == 'T') then 
      do i=1, n_row
        desc_a%lprm(i) = iperm(i)
      enddo
      do i=n_row+1,n_col
        desc_a%lprm(i) = i
      enddo
    endif
    ! crossed fingers.....
    ! fix glob_to_loc/loc_to_glob  mappings, then indices lists
    ! hmm, maybe we should just moe all of this onto a different level,
    ! have a specialized subroutine, and do it in the solver context???? 
    if (debug) write(0,*) 'spasb: renumbering glob_to_loc'
    do i=1, n_col
      desc_a%glob_to_loc(desc_a%loc_to_glob(desc_a%lprm(i))) = i  
    enddo
    if (debug) write(0,*) 'spasb: renumbering loc_to_glob'
    do i=1,psb_cd_get_global_rows(desc_a) 
      j = desc_a%glob_to_loc(i)
      if (j>0) then 
        desc_a%loc_to_glob(j) = i
      endif
    enddo
    if (debug) write(0,*) 'spasb: renumbering halo_index'
    i=1
    kh=desc_a%halo_index(i)
    do while (kh /= -1) 
      i = i+1
      nh = desc_a%halo_index(i)
      do j = i+1, i+nh
        desc_a%halo_index(j) = &
             &desc_a%lprm(desc_a%halo_index(j))
      enddo
      i = i + nh + 1
      nh = desc_a%halo_index(i)
      do j= i+1, i+nh
        desc_a%halo_index(j) = &
             &desc_a%lprm(desc_a%halo_index(j))
      enddo
      i = i + nh + 1
      kh=desc_a%halo_index(i)
    enddo
    if (debug) write(0,*) 'spasb: renumbering ovrlap_index'
    i=1
    kh=desc_a%ovrlap_index(i)
    do while (kh /= -1) 
      i = i + 1
      nh = desc_a%ovrlap_index(i)
      do j= i+1, i+nh
        desc_a%ovrlap_index(j) = &
             &desc_a%lprm(desc_a%ovrlap_index(j))
      enddo
      i = i + nh + 1
      kh=desc_a%ovrlap_index(i)
    enddo
    if (debug) write(0,*) 'spasb: renumbering ovrlap_elem'
    i = 1
    kh=desc_a%ovrlap_elem(i)
    do while (kh /= -1)          
      desc_a%ovrlap_elem(i) = &
           &desc_a%lprm(desc_a%ovrlap_elem(i))
      i = i+2
      kh=desc_a%ovrlap_elem(i)         
    enddo
    if (debug) write(0,*) 'spasb: done renumbering'
    if (debug) then
      write(60+me,*) 'n_row ',n_row,' n_col',n_col, ' trans: ',trans
      do i=1,n_col
        write(60+me,*)i, ' lprm ', desc_a%lprm(i), ' iperm',iperm(i)
      enddo
      i=1
      kh = desc_a%halo_index(i)
      do while (kh /= -1) 
        write(60+me,*) i, kh 
        i = i+1
        kh = desc_a%halo_index(i)
      enddo
      close(60+me)
    end if

!!$    iperm(1) = 0
  else 
!!$    allocate(desc_a%lprm(1))
!!$    desc_a%lprm(1) = 0       
  endif


  time(4) = mpi_wtime()
  time(4) = time(4) - time(3)
  if (debug) then 
    call psb_amx(ictxt, time(4))

    write (*, *) '         comm structs assembly: ', time(4)*1.d-3
  end if

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act == psb_act_ret_) then
    return
  else
    call psb_error(ictxt)
  end if
  return

end subroutine psb_cdren
