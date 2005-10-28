! File: psb_chkmat.f90
! File: psb_chkmat.f90
!
! Subroutine: psb_chkmat
!    pbmatvect checks the validity of a descriptor vector DESCDEC, the
!    related global indexes IA, JA. It also computes the starting local
!    indexes (IIA,JJA) corresponding to the submatrix starting globally at
!    the entry pointed by (IA,JA). Finally, if an inconsitency is found among 
!    its parameters ia, ja and desc_A, the routine returns an error code in
!    info.
!
! Parameters:
!  m        - integer.               The number of rows of the matrix being operated on.    
!  n        - integer.               The number of columns of the matrix being operated on.    
!  ia       - integer.               a's global row index, which points to the beginning 
!                                    of the submatrix which is to be operated on.      
!  ja       - integer.               a's global column index, which points to the beginning 
!                                    of the submatrix which is to be operated on.      
!  desc_dec - integer,dimension(:).  Is the matrix_data array.
!  info     - integer.               Eventually returns an error code.
!  iia      - integer(optional).     The local rows starting index of the submatrix.
!  jja      - integer(optional).     The local columns starting index of the submatrix.
!
subroutine psb_chkmat( m, n, ia, ja, desc_dec, info, iia, jja)

  use psb_error_mod
  implicit none

  integer, intent(in)    ::  m,n,ia,ja
  integer, intent(in)    ::  desc_dec(:)
  integer, intent(out)   ::  info
  integer, optional      ::  iia, jja

  !  locals
  integer           :: err_act, int_err(5)
  character(len=20) :: name, ch_err

  if(psb_get_errstatus().ne.0) return 
  info=0
  name='psb_chkmat'
  call psb_erractionsave(err_act)

  if (m.lt.0) then
     info=10
     int_err(1) = 1
     int_err(2) = m
  else if (n.lt.0) then
     info=10
     int_err(1) = 3
     int_err(2) = n
  else if ((ix.lt.1) .and. (m.ne.0)) then
     info=20
     int_err(1) = 4
     int_err(2) = ix
  else if ((jx.lt.1) .and. (n.ne.0)) then
     info=20
     int_err(1) = 5
     int_err(2) = jx
  else if (desc_dec(psb_n_col_).lt.0) then
     info=40
     int_err(1) = 6
     int_err(2) = psb_n_col_ 
     int_err(3) = desc_dec(psb_n_col_)
  else if (desc_dec(psb_n_row_).lt.0) then
     info=40
     int_err(1) = 6
     int_err(2) = psb_n_row_ 
     int_err(3) = desc_dec(psb_n_row_)
  else if (desc_dec(m_).lt.m) then
     info=60
     int_err(1) = 1
     int_err(2) = m
     int_err(3) = 5
     int_err(4) = m_
     int_err(5) = desc_dec(m_)
  else if (desc_dec(m_).lt.m) then
     info=60
     int_err(1) = 2
     int_err(2) = n
     int_err(3) = 5
     int_err(4) = n_
     int_err(5) = desc_dec(n_)
  else if (desc_dec(m_).lt.ia) then
     info=60
     int_err(1) = 3 
     int_err(2) = ia
     int_err(3) = 5
     int_err(4) = m_
     int_err(5) = desc_dec(m_)
  else if (desc_dec(n_).lt.ja) then
     info=60
     int_err(1) = 4 
     int_err(2) = ja
     int_err(3) = 5
     int_err(4) = n_
     int_err(5) = desc_dec(n_)
  else if (desc_dec(m_).lt.(ia+m-1)) then
     info=80
     int_err(1) = 1
     int_err(2) = m
     int_err(3) = 3
     int_err(4) = ia
  else if (desc_dec(n_).lt.(ja+n-1)) then
     info=80
     int_err(1) = 2
     int_err(2) = n
     int_err(3) = 4
     int_err(4) = ja
  end if

  if (info.ne.0) then
     call psb_errpush(info,name,i_err=int_err)
     goto 9999
  end if

  ! Compute local indices for submatrix starting
  ! at global indices ix and jx
  if(present(iia).and.present(jja)) then
     if (desc_dec(psb_n_row_).gt.0) then
        iia=1
        jja=1
     else
        iia=desc_dec(psb_n_row_)+1
        jja=desc_dec(psb_n_col_)+1
     end if
  end if

  call psb_erractionrestore(err_act)
  return  

9999 continue
  call psb_erractionrestore(err_act)

  if (err_act.eq.act_abort) then
     call psb_error()
     return
  end if
  return
end subroutine psb_chkmat
