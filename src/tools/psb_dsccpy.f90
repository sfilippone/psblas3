! File: psb_dsccpy.f90
!
! Subroutine: psb_dsccpy
!   Produces a clone of a descriptor.
! 
! Parameters: 
!    desc_out - type(<psb_desc_type>).         The output communication descriptor.
!    desc_a   - type(<psb_desc_type>).         The communication descriptor to be cloned.
!    info     - integer.                       Eventually returns an error code.
subroutine psb_dsccpy(desc_out, desc_a, info)

  use psb_descriptor_type
  use psb_serial_mod
  use psb_realloc_mod
  use psb_const_mod
  use psb_error_mod

  implicit none
  !....parameters...

  type(psb_desc_type), intent(out) :: desc_out
  type(psb_desc_type), intent(in)  :: desc_a
  integer, intent(out)             :: info

  !locals
  integer             :: nprow,npcol,me,mypcol,&
       & icontxt, isz, dectype, err_act, err
  integer             :: int_err(5),temp(1)
  real(kind(1.d0))    :: real_err(5)
  integer, parameter  :: ione=1, itwo=2,root=0
  logical, parameter  :: debug=.false.
  character(len=20)   :: name, char_err

  info=0
  call psb_erractionsave(err_act)
  name = 'psb_dsccpy'

  icontxt=desc_a%matrix_data(psb_ctxt_)
  
  ! check on blacs grid 
  call blacs_gridinfo(icontxt, nprow, npcol, me, mypcol)
  if (nprow.eq.-1) then
     info = 2010
     call psb_errpush(info,name)
     goto 9999
  else if (npcol.ne.1) then
     info = 2030
     int_err(1) = npcol
     call psb_errpush(info,name,int_err)
     goto 9999
  endif

  call psb_nullify_desc(desc_out)

  if (associated(desc_a%matrix_data)) then 
     isz = size(desc_a%matrix_data)
     !    allocate(desc_out%matrix_data(isz),stat=info)
     call psb_realloc(isz,desc_out%matrix_data,info)
     if(debug) write(0,*) 'dsccpy: m_data',isz,':',desc_a%matrix_data(:)
     if (info.ne.0) then     
        info=4010
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
     else
        desc_out%matrix_data(:) = desc_a%matrix_data(:)
     endif
  endif
  
  if (associated(desc_a%halo_index)) then 
     isz = size(desc_a%halo_index)
     !    allocate(desc_out%matrix_data(isz),stat=info)
     call psb_realloc(isz,desc_out%halo_index,info)
     if(debug) write(0,*) 'dsccpy: h_idx',isz,':',desc_a%halo_index(:)
     if (info.ne.0) then     
        info=4010
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
     else
        desc_out%halo_index(:) = desc_a%halo_index(:)
     endif
  endif

  
  if (associated(desc_a%bnd_elem)) then 
     isz = size(desc_a%bnd_elem)
     !    allocate(desc_out%matrix_data(isz),stat=info)
     call psb_realloc(isz,desc_out%bnd_elem,info)
     if(debug) write(0,*) 'dsccpy: bnd_elem',isz,':',desc_a%bnd_elem(:)
     if (info.ne.0) then     
        info=4010
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
     else
        desc_out%bnd_elem(:) = desc_a%bnd_elem(:)
     endif
  endif

  
  if (associated(desc_a%ovrlap_elem)) then 
     isz = size(desc_a%ovrlap_elem)
     !    allocate(desc_out%matrix_data(isz),stat=info)
     call psb_realloc(isz,desc_out%ovrlap_elem,info)
     if(debug) write(0,*) 'dsccpy: ovrlap_elem',isz,':',desc_a%ovrlap_elem(:)
     if (info.ne.0) then     
        info=4010
        char_err='psrealloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
     else
        desc_out%ovrlap_elem(:) = desc_a%ovrlap_elem(:)
     endif
  endif

  if (associated(desc_a%ovrlap_index)) then 
     isz = size(desc_a%ovrlap_index)
     !    allocate(desc_out%matrix_data(isz),stat=info)
     call psb_realloc(isz,desc_out%ovrlap_index,info)
     if(debug) write(0,*) 'dsccpy: ovrlap_index',isz,':',desc_a%ovrlap_index(:)
     if (info.ne.0) then     
        info=4010
        char_err='psrealloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
     else
        desc_out%ovrlap_index(:) = desc_a%ovrlap_index(:)
     endif
  endif


  if (associated(desc_a%loc_to_glob)) then 
     isz = size(desc_a%loc_to_glob)
     !    allocate(desc_out%matrix_data(isz),stat=info)
     call psb_realloc(isz,desc_out%loc_to_glob,info)
     if(debug) write(0,*) 'dsccpy: loc_to_glob',isz,':',desc_a%loc_to_glob(:)
     if (info.ne.0) then     
        info=4010
        char_err='psrealloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
     else
        desc_out%loc_to_glob(:) = desc_a%loc_to_glob(:)
     endif
  endif

  if (associated(desc_a%glob_to_loc)) then 
     isz = size(desc_a%glob_to_loc)
     !    allocate(desc_out%matrix_data(isz),stat=info)
     call psb_realloc(isz,desc_out%glob_to_loc,info)
     if(debug) write(0,*) 'dsccpy: glob_to_loc',isz,':',desc_a%glob_to_loc(:)
     if (info.ne.0) then     
        info=4010
        char_err='psrealloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
     else
        desc_out%glob_to_loc(:) = desc_a%glob_to_loc(:)
     endif
  endif

  if (associated(desc_a%lprm)) then 
     isz = size(desc_a%lprm)
     !    allocate(desc_out%matrix_data(isz),stat=info)
     call psb_realloc(isz,desc_out%lprm,info)
     if(debug) write(0,*) 'dsccpy: lprm',isz,':',desc_a%lprm(:)
     if (info.ne.0) then     
        info=4010
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
     else
        desc_out%lprm(:) = desc_a%lprm(:)
     endif
  endif

  if (associated(desc_a%idx_space)) then 
     isz = size(desc_a%idx_space)
     !    allocate(desc_out%matrix_data(isz),stat=info)
     call psb_realloc(isz,desc_out%idx_space,info)
     if(debug) write(0,*) 'dsccpy: idx_space',isz,':',desc_a%idx_space(:)
     if (info.ne.0) then     
        info=4010
        char_err='psb_realloc'
        call psb_errpush(info,name,a_err=char_err)
        goto 9999
     else
        desc_out%idx_space(:) = desc_a%idx_space(:)
     endif
  endif

  call psb_erractionrestore(err_act)
  return
  
9999 continue
  call psb_erractionrestore(err_act)
  
  if (err_act.eq.act_ret) then
     return
  else
     call psb_error(icontxt)
  end if
  return

end subroutine psb_dsccpy
