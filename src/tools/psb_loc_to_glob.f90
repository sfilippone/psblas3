! File: psb_loc_to_glob.f90
!
! Subroutine: psb_loc_to_glob2
!    Performs local to global indexes translation
! 
! Parameters: 
!    x        - integer, dimension(:).    Array containing the indices to be translated.
!    y        - integer, dimension(:).    Array containing the indices to be translated.
!    desc_a   - type(<psb_desc_type>).    The communication descriptor.        
!    info     - integer.                  Eventually returns an error code.
!    iact     - integer(optional).        A character defining the behaviour of this subroutine when is found an index not belonging to the calling process
!
subroutine psb_loc_to_glob2(x,y,desc_a,info,iact)

  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  implicit none

  !...parameters....
  type(psb_desc_type), intent(in)    ::  desc_a
  integer, intent(in)                ::  x(:)  
  integer, intent(out)               ::  y(:)  
  integer, intent(out)               ::  info
  character, intent(in), optional    ::  iact

  !....locals....
  integer                            ::  err, n, i, tmp, icontxt
  character                          ::  strings,act
  integer                            ::  int_err(5), err_act
  real(kind(1.d0))                   ::  real_val
  integer, parameter                 ::  zero=0
  character(len=20)   :: name, char_err

  info=0
  name='psb_loc_to_glob2'
  call psb_erractionsave(err_act)

  if (present(iact)) then
     act=iact
  else
     act='A'
  endif   
  
  real_val = 0.d0

  n=size(x)
  do i=1,n
     if ((x(i).gt.desc_a%matrix_data(psb_n_col_)).or.&
          &  (x(i).le.zero)) then
        info=140
        int_err(1)=tmp
        int_err(2)=desc_a%matrix_data(m_)
        exit
     else
        tmp=desc_a%loc_to_glob(x(i))
        if((tmp.gt.zero).or.(tmp.le.desc_a%matrix_data(m_))) then
           y(i)=tmp
        else
           info = 140
           int_err(1)=tmp
           int_err(2)=desc_a%matrix_data(psb_n_col_)
           exit
        end if
     end if
  enddo
  
  if (info.ne.0) then
     select case(act)
     case('E')
        call psb_erractionrestore(err_act)
        return
     case('W')
        write(0,'("Error ",i5," in subroutine glob_to_loc")') info
     case('A')
        call psb_errpush(info,name)
        goto 9999
     end select
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_ret) then
     return
  else
     call psb_error()
  end if
  return

end subroutine psb_loc_to_glob2


! Subroutine: psb_loc_to_glob
!    Performs local to global indexes translation
! 
! Parameters: 
!    x        - integer, dimension(:).    Array containing the indices to be translated.
!    desc_a   - type(<psb_desc_type>).    The communication descriptor.        
!    info     - integer.                  Eventually returns an error code.
!    iact     - integer(optional).        A character defining the behaviour of this subroutine when is found an index not belonging to the calling process
!
subroutine psb_loc_to_glob(x,desc_a,info,iact)

  use psb_descriptor_type
  use psb_const_mod
  use psb_error_mod
  implicit none

  !...parameters....
  type(psb_desc_type), intent(in)    ::  desc_a
  integer, intent(inout)             ::  x(:)  
  integer, intent(out)               ::  info
  character, intent(in), optional    ::  iact

  !....locals....
  integer                            ::  err, n ,i, tmp,icontxt, err_act
  character                          ::  act
  integer                            ::  int_err(5)
  real(kind(1.d0))                   ::  real_val
  integer, parameter                 ::  zero=0
  character(len=20)   :: name, char_err

  info=0
  name='psb_loc_to_glob'
  call psb_erractionsave(err_act)

  if (present(iact)) then
     act=iact
  else
     act='A'
  endif   

  real_val = 0.d0

  n=size(x)
  do i=1,n
    if ((x(i).gt.desc_a%matrix_data(psb_n_col_)).or.&
         &  (x(i).le.zero)) then
      info=140
      int_err(1)=x(i)
      int_err(2)=desc_a%matrix_data(psb_n_col_)  
      exit
    else
      tmp=desc_a%loc_to_glob(x(i))
      if((tmp.gt.zero).or.(tmp.le.desc_a%matrix_data(m_))) then
        x(i)=tmp
      else
        info = 140
        exit
      end if
    end if
  enddo

  if (info.ne.0) then
     select case(act)
     case('E')
        call psb_erractionrestore(err_act)
        return
     case('W')
        write(0,'("Error ",i5," in subroutine glob_to_loc")') info
     case('A')
        call psb_errpush(info,name)
        goto 9999
     end select
  endif

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_ret) then
     return
  else
     call psb_error()
  end if
  return

end subroutine psb_loc_to_glob

