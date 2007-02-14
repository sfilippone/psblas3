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
! File: psb_check_mod.f90

module psb_check_mod

!   interface
!      module procedure psb_chkvect
!   end interface

!   interface
!      module procedure psb_chkglobvect
!   end interface

!   interface
!      module procedure psb_chkmat
!   end interface

contains
  ! Subroutine: psb_chkvect
  !    psb_chkvect checks the validity of a descriptor vector desc_dec, the
  !    related global indexes ix, jx and the leading dimension lldx. It also
  !    eventually computes the starting local indexes (iix,jjx) corresponding
  !    to the submatrix starting globally at the entry pointed by (ix,jx).
  !    Finally, if an inconsistency is found among its parameters ix, jx,
  !    descdec and lldx, the routine returns an error code in info.
  !
  ! Parameters:
  !  m        - integer.               The number of rows of the dense matrix X being operated on.    
  !  n        - integer.               The number of columns of the dense matrix X being operated on.    
  !  lldx     - integer.               The leading dimension of the local dense matrix X.
  !  ix       - integer.               X's global row index, which points to the beginning 
  !                                    of the dense submatrix which is to be operated on.      
  !  jx       - integer.               X's global column index, which points to the beginning 
  !                                    of the dense submatrix which is to be operated on.      
  !  desc_dec - integer,dimension(:).  Is the matrix_data array.
  !  info     - integer.               Eventually returns an error code.
  !  iix      - integer(optional).     The local rows starting index of the submatrix.
  !  jjx      - integer(optional).     The local columns starting index of the submatrix.
  subroutine psb_chkvect( m, n, lldx, ix, jx, desc_dec, info, iix, jjx)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    integer, intent(in)    ::  m,n,ix,jx,lldx
    type(psb_desc_type), intent(in)    ::  desc_dec
    integer, intent(out)   ::  info
    integer, optional      ::  iix, jjx

    !  locals
    integer           :: err_act, int_err(5)
    character(len=20) :: name, ch_err

    if(psb_get_errstatus() /= 0) return 
    info=0
    name='psb_chkvect'
    call psb_erractionsave(err_act)


    if (m < 0) then
       info=10
       int_err(1) = 1
       int_err(2) = m
    else if (n < 0) then
       info=10
       int_err(1) = 3
       int_err(2) = n
    else if ((ix < 1) .and. (m /= 0)) then
       info=20
       int_err(1) = 4
       int_err(2) = ix
    else if ((jx < 1) .and. (n /= 0)) then
       info=20
       int_err(1) = 5
       int_err(2) = jx
    else if (psb_cd_get_local_cols(desc_dec) < 0) then
       info=40
       int_err(1) = 6
       int_err(2) = psb_n_col_ 
       int_err(3) = psb_cd_get_local_cols(desc_dec)
    else if (psb_cd_get_local_rows(desc_dec) < 0) then
       info=40
       int_err(1) = 6
       int_err(2) = psb_n_row_ 
       int_err(3) = psb_cd_get_local_cols(desc_dec)
    else if (lldx < psb_cd_get_local_cols(desc_dec)) then
       info=50
       int_err(1) = 3
       int_err(2) = lldx
       int_err(3) = 6
       int_err(4) = psb_n_col_
       int_err(5) = psb_cd_get_local_cols(desc_dec)
    else if (psb_cd_get_global_cols(desc_dec) < m) then
       info=60
       int_err(1) = 1
       int_err(2) = m
       int_err(3) = 6
       int_err(4) = psb_n_
       int_err(5) = psb_cd_get_global_cols(desc_dec)
    else if (psb_cd_get_global_cols(desc_dec) < ix) then
       info=60
       int_err(1) = 4
       int_err(2) = ix
       int_err(3) = 6
       int_err(4) = psb_n_
       int_err(5) = psb_cd_get_global_cols(desc_dec)
    else if (psb_cd_get_global_rows(desc_dec) < jx) then
       info=60
       int_err(1) = 5
       int_err(2) = jx
       int_err(3) = 6
       int_err(4) = psb_m_
       int_err(5) = psb_cd_get_global_rows(desc_dec)
    else if (psb_cd_get_global_cols(desc_dec) < (ix+m-1)) then
       info=80
       int_err(1) = 1
       int_err(2) = m
       int_err(3) = 4
       int_err(4) = ix
    end if

    if (info /= 0) then
       call psb_errpush(info,name,i_err=int_err)
       goto 9999
    end if

    ! Compute local indices for submatrix starting
    ! at global indices ix and jx
    if(present(iix)) iix=ix  ! (for our applications iix=ix))
    if(present(jjx)) jjx=ix  ! (for our applications jjx=jx))

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
       call psb_error()
       return
    end if
    return

  end subroutine psb_chkvect

  !
  ! Subroutine: psb_chkglobvect
  !    psb_chkglobvect checks the validity of a descriptor vector desc_dec, the
  !    related global indexes ix, jx and the leading dimension lldx.
  !    If an inconsistency is found among its parameters ix, jx,
  !    descdec and lldx, the routine returns an error code in info.
  !
  ! Parameters:
  !  m        - integer.               The number of rows of the dense matrix X being operated on.    
  !  n        - integer.               The number of columns of the dense matrix X being operated on.    
  !  lldx     - integer.               The leading dimension of the local dense matrix X.
  !  ix       - integer.               X's global row index, which points to the beginning 
  !                                    of the dense submatrix which is to be operated on.      
  !  jx       - integer.               X's global column index, which points to the beginning 
  !                                    of the dense submatrix which is to be operated on.      
  !  desc_dec - integer,dimension(:).  Is the matrix_data array.
  !  info     - integer.               Eventually returns an error code.
  !
  subroutine psb_chkglobvect( m, n, lldx, ix, jx, desc_dec, info)

    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    integer, intent(in)    ::  m,n,ix,jx,lldx
    type(psb_desc_type), intent(in)    ::  desc_dec
    integer, intent(out)   ::  info

    !  locals
    integer           :: err_act, int_err(5)
    character(len=20) :: name, ch_err

    if(psb_get_errstatus() /= 0) return 
    info=0
    name='psb_chkglobvect'
    call psb_erractionsave(err_act)


    if (m < 0) then
       info=10
       int_err(1) = 1
       int_err(2) = m
    else if (n < 0) then
       info=10
       int_err(1) = 3
       int_err(2) = n
    else if ((ix < 1) .and. (m /= 0)) then
       info=20
       int_err(1) = 4
       int_err(2) = ix
    else if ((jx < 1) .and. (n /= 0)) then
       info=20
       int_err(1) = 5
       int_err(2) = jx
    else if (psb_cd_get_local_cols(desc_dec) < 0) then
       info=40
       int_err(1) = 6
       int_err(2) = psb_n_col_ 
       int_err(3) = psb_cd_get_local_cols(desc_dec)
    else if (psb_cd_get_local_rows(desc_dec) < 0) then
       info=40
       int_err(1) = 6
       int_err(2) = psb_n_row_ 
       int_err(3) = psb_cd_get_local_rows(desc_dec)
    else if (lldx < psb_cd_get_global_rows(desc_dec)) then
       info=50
       int_err(1) = 3
       int_err(2) = lldx
       int_err(3) = 6
       int_err(4) = psb_n_col_
       int_err(5) = psb_cd_get_global_rows(desc_dec)
    else if (psb_cd_get_global_cols(desc_dec) < m) then
       info=60
       int_err(1) = 1
       int_err(2) = m
       int_err(3) = 6
       int_err(4) = psb_n_
       int_err(5) = psb_cd_get_global_cols(desc_dec)
    else if (psb_cd_get_global_cols(desc_dec) < ix) then
       info=60
       int_err(1) = 4
       int_err(2) = ix
       int_err(3) = 6
       int_err(4) = psb_n_
       int_err(5) = psb_cd_get_global_cols(desc_dec)
    else if (psb_cd_get_global_rows(desc_dec) < jx) then
       info=60
       int_err(1) = 5
       int_err(2) = jx
       int_err(3) = 6
       int_err(4) = psb_m_
       int_err(5) = psb_cd_get_global_rows(desc_dec)
    else if (psb_cd_get_global_cols(desc_dec) < (ix+m-1)) then
       info=80
       int_err(1) = 1
       int_err(2) = m
       int_err(3) = 4
       int_err(4) = ix
    end if

    if (info /= 0) then
       call psb_errpush(info,name,i_err=int_err)
       goto 9999
    end if

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
       call psb_error()
       return
    end if
    return

  end subroutine psb_chkglobvect

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

    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    integer, intent(in)    ::  m,n,ia,ja
    type(psb_desc_type), intent(in)    ::  desc_dec
    integer, intent(out)   ::  info
    integer, optional      ::  iia, jja

    !  locals
    integer           :: err_act, int_err(5)
    character(len=20) :: name, ch_err

    if(psb_get_errstatus() /= 0) return 
    info=0
    name='psb_chkmat'
    call psb_erractionsave(err_act)

    if (m < 0) then
      info=10
      int_err(1) = 1
      int_err(2) = m
    else if (n < 0) then
      info=10
      int_err(1) = 3
      int_err(2) = n
    else if ((ia < 1) .and. (m /= 0)) then
      info=20
      int_err(1) = 4
      int_err(2) = ia
    else if ((ja < 1) .and. (n /= 0)) then
      info=20
      int_err(1) = 5
      int_err(2) = ja
    else if (psb_cd_get_local_cols(desc_dec) < 0) then
      info=40
      int_err(1) = 6
      int_err(2) = psb_n_col_ 
      int_err(3) = psb_cd_get_local_cols(desc_dec)
    else if (psb_cd_get_local_rows(desc_dec) < 0) then
      info=40
      int_err(1) = 6
      int_err(2) = psb_n_row_ 
      int_err(3) = psb_cd_get_local_rows(desc_dec)
    else if (psb_cd_get_global_rows(desc_dec) < m) then
      info=60
      int_err(1) = 1
      int_err(2) = m
      int_err(3) = 5
      int_err(4) = psb_m_
      int_err(5) = psb_cd_get_global_rows(desc_dec)
    else if (psb_cd_get_global_rows(desc_dec) < m) then
      info=60
      int_err(1) = 2
      int_err(2) = n
      int_err(3) = 5
      int_err(4) = psb_m_
      int_err(5) = psb_cd_get_global_rows(desc_dec)
    else if (psb_cd_get_global_rows(desc_dec) < ia) then
      info=60
      int_err(1) = 3 
      int_err(2) = ia
      int_err(3) = 5
      int_err(4) = psb_m_
      int_err(5) = psb_cd_get_global_rows(desc_dec)
    else if (psb_cd_get_global_cols(desc_dec) < ja) then
      info=60
      int_err(1) = 4 
      int_err(2) = ja
      int_err(3) = 5
      int_err(4) = psb_n_
      int_err(5) = psb_cd_get_global_cols(desc_dec)
    else if (psb_cd_get_global_rows(desc_dec) < (ia+m-1)) then
      info=80
      int_err(1) = 1
      int_err(2) = m
      int_err(3) = 3
      int_err(4) = ia
    else if (psb_cd_get_global_cols(desc_dec) < (ja+n-1)) then
      info=80
      int_err(1) = 2
      int_err(2) = n
      int_err(3) = 4
      int_err(4) = ja
    end if

    if (info /= 0) then
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

    ! Compute local indices for submatrix starting
    ! at global indices ix and jx
    if(present(iia).and.present(jja)) then
      if (psb_cd_get_local_rows(desc_dec) > 0) then
        iia=1
        jja=1
      else
        iia=psb_cd_get_local_rows(desc_dec)+1
        jja=psb_cd_get_local_cols(desc_dec)+1
      end if
    end if

    call psb_erractionrestore(err_act)
    return  

9999 continue
    call psb_erractionrestore(err_act)

    if (err_act == psb_act_abort_) then
       call psb_error()
       return
    end if
    return
  end subroutine psb_chkmat

end module psb_check_mod
