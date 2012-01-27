!!$ 
!!$              Parallel Sparse BLAS  version 3.0
!!$    (C) Copyright 2006, 2007, 2008, 2009, 2010
!!$                       Salvatore Filippone    University of Rome Tor Vergata
!!$                       Alfredo Buttari        CNRS-IRIT, Toulouse
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
  !  desc_dec - integer(psb_ipk_),dimension(:).  Is the matrix_data array.
  !  info     - integer.               Return code
  !  iix      - integer(optional).     The local rows starting index of the submatrix.
  !  jjx      - integer(optional).     The local columns starting index of the submatrix.
  subroutine psb_chkvect( m, n, lldx, ix, jx, desc_dec, info, iix, jjx)
    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    integer(psb_ipk_), intent(in)    ::  m,n,ix,jx,lldx
    type(psb_desc_type), intent(in)    ::  desc_dec
    integer(psb_ipk_), intent(out)   ::  info
    integer(psb_ipk_), optional      ::  iix, jjx

    !  locals
    integer(psb_ipk_) :: err_act, int_err(5)
    character(len=20) :: name

    if(psb_get_errstatus() /= 0) return 
    info=psb_success_
    name='psb_chkvect'
    call psb_erractionsave(err_act)


    if (m < 0) then
       info=psb_err_iarg_neg_
       int_err(1) = 1
       int_err(2) = m
    else if (n < 0) then
       info=psb_err_iarg_neg_
       int_err(1) = 3
       int_err(2) = n
    else if ((ix < 1) .and. (m /= 0)) then
       info=psb_err_iarg_pos_
       int_err(1) = 4
       int_err(2) = ix
    else if ((jx < 1) .and. (n /= 0)) then
       info=psb_err_iarg_pos_
       int_err(1) = 5
       int_err(2) = jx
    else if (desc_dec%get_local_cols() < 0) then
       info=psb_err_iarg_invalid_i_
       int_err(1) = 6
       int_err(2) = psb_n_col_ 
       int_err(3) = desc_dec%get_local_cols()
    else if (desc_dec%get_local_rows() < 0) then
       info=psb_err_iarg_invalid_i_
       int_err(1) = 6
       int_err(2) = psb_n_row_ 
       int_err(3) = desc_dec%get_local_cols()
    else if (lldx < desc_dec%get_local_cols()) then
       info=psb_err_iarg_not_gtia_ii_
       int_err(1) = 3
       int_err(2) = lldx
       int_err(3) = 6
       int_err(4) = psb_n_col_
       int_err(5) = desc_dec%get_local_cols()
    else if (desc_dec%get_global_cols() < m) then
       info=psb_err_iarg_not_gteia_ii_
       int_err(1) = 1
       int_err(2) = m
       int_err(3) = 6
       int_err(4) = psb_n_
       int_err(5) = desc_dec%get_global_cols()
    else if (desc_dec%get_global_cols() < ix) then
       info=psb_err_iarg_not_gteia_ii_
       int_err(1) = 4
       int_err(2) = ix
       int_err(3) = 6
       int_err(4) = psb_n_
       int_err(5) = desc_dec%get_global_cols()
    else if (desc_dec%get_global_rows() < jx) then
       info=psb_err_iarg_not_gteia_ii_
       int_err(1) = 5
       int_err(2) = jx
       int_err(3) = 6
       int_err(4) = psb_m_
       int_err(5) = desc_dec%get_global_rows()
    else if (desc_dec%get_global_cols() < (ix+m-1)) then
       info=psb_err_iarg2_neg_
       int_err(1) = 1
       int_err(2) = m
       int_err(3) = 4
       int_err(4) = ix
    end if

    if (info /= psb_success_) then
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
  !  desc_dec - integer(psb_ipk_),dimension(:).  Is the matrix_data array.
  !  info     - integer.               Return code
  !
  subroutine psb_chkglobvect( m, n, lldx, ix, jx, desc_dec, info)

    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    integer(psb_ipk_), intent(in)    ::  m,n,ix,jx,lldx
    type(psb_desc_type), intent(in)    ::  desc_dec
    integer(psb_ipk_), intent(out)   ::  info

    !  locals
    integer(psb_ipk_) :: err_act, int_err(5)
    character(len=20) :: name

    if(psb_get_errstatus() /= 0) return 
    info=psb_success_
    name='psb_chkglobvect'
    call psb_erractionsave(err_act)


    if (m < 0) then
       info=psb_err_iarg_neg_
       int_err(1) = 1
       int_err(2) = m
    else if (n < 0) then
       info=psb_err_iarg_neg_
       int_err(1) = 3
       int_err(2) = n
    else if ((ix < 1) .and. (m /= 0)) then
       info=psb_err_iarg_pos_
       int_err(1) = 4
       int_err(2) = ix
    else if ((jx < 1) .and. (n /= 0)) then
       info=psb_err_iarg_pos_
       int_err(1) = 5
       int_err(2) = jx
    else if (desc_dec%get_local_cols() < 0) then
       info=psb_err_iarg_invalid_i_
       int_err(1) = 6
       int_err(2) = psb_n_col_ 
       int_err(3) = desc_dec%get_local_cols()
    else if (desc_dec%get_local_rows() < 0) then
       info=psb_err_iarg_invalid_i_
       int_err(1) = 6
       int_err(2) = psb_n_row_ 
       int_err(3) = desc_dec%get_local_rows()
    else if (lldx < desc_dec%get_global_rows()) then
       info=psb_err_iarg_not_gtia_ii_
       int_err(1) = 3
       int_err(2) = lldx
       int_err(3) = 6
       int_err(4) = psb_n_col_
       int_err(5) = desc_dec%get_global_rows()
    else if (desc_dec%get_global_cols() < m) then
       info=psb_err_iarg_not_gteia_ii_
       int_err(1) = 1
       int_err(2) = m
       int_err(3) = 6
       int_err(4) = psb_n_
       int_err(5) = desc_dec%get_global_cols()
    else if (desc_dec%get_global_cols() < ix) then
       info=psb_err_iarg_not_gteia_ii_
       int_err(1) = 4
       int_err(2) = ix
       int_err(3) = 6
       int_err(4) = psb_n_
       int_err(5) = desc_dec%get_global_cols()
    else if (desc_dec%get_global_rows() < jx) then
       info=psb_err_iarg_not_gteia_ii_
       int_err(1) = 5
       int_err(2) = jx
       int_err(3) = 6
       int_err(4) = psb_m_
       int_err(5) = desc_dec%get_global_rows()
    else if (desc_dec%get_global_cols() < (ix+m-1)) then
       info=psb_err_iarg2_neg_
       int_err(1) = 1
       int_err(2) = m
       int_err(3) = 4
       int_err(4) = ix
    end if

    if (info /= psb_success_) then
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
  !  desc_dec - integer(psb_ipk_),dimension(:).  Is the matrix_data array.
  !  info     - integer.               Return code
  !  iia      - integer(optional).     The local rows starting index of the submatrix.
  !  jja      - integer(optional).     The local columns starting index of the submatrix.
  !
  subroutine psb_chkmat( m, n, ia, ja, desc_dec, info, iia, jja)

    use psb_descriptor_type
    use psb_const_mod
    use psb_error_mod
    implicit none

    integer(psb_ipk_), intent(in)    ::  m,n,ia,ja
    type(psb_desc_type), intent(in)    ::  desc_dec
    integer(psb_ipk_), intent(out)   ::  info
    integer(psb_ipk_), optional      ::  iia, jja

    !  locals
    integer(psb_ipk_) :: err_act, int_err(5)
    character(len=20) :: name

    if(psb_get_errstatus() /= 0) return 
    info=psb_success_
    name='psb_chkmat'
    call psb_erractionsave(err_act)

    if (m < 0) then
      info=psb_err_iarg_neg_
      int_err(1) = 1
      int_err(2) = m
    else if (n < 0) then
      info=psb_err_iarg_neg_
      int_err(1) = 3
      int_err(2) = n
    else if ((ia < 1) .and. (m /= 0)) then
      info=psb_err_iarg_pos_
      int_err(1) = 4
      int_err(2) = ia
    else if ((ja < 1) .and. (n /= 0)) then
      info=psb_err_iarg_pos_
      int_err(1) = 5
      int_err(2) = ja
    else if (desc_dec%get_local_cols() < 0) then
      info=psb_err_iarg_invalid_i_
      int_err(1) = 6
      int_err(2) = psb_n_col_ 
      int_err(3) = desc_dec%get_local_cols()
    else if (desc_dec%get_local_rows() < 0) then
      info=psb_err_iarg_invalid_i_
      int_err(1) = 6
      int_err(2) = psb_n_row_ 
      int_err(3) = desc_dec%get_local_rows()
    else if (desc_dec%get_global_rows() < m) then
      info=psb_err_iarg_not_gteia_ii_
      int_err(1) = 1
      int_err(2) = m
      int_err(3) = 5
      int_err(4) = psb_m_
      int_err(5) = desc_dec%get_global_rows()
    else if (desc_dec%get_global_rows() < m) then
      info=psb_err_iarg_not_gteia_ii_
      int_err(1) = 2
      int_err(2) = n
      int_err(3) = 5
      int_err(4) = psb_m_
      int_err(5) = desc_dec%get_global_rows()
    else if (desc_dec%get_global_rows() < ia) then
      info=psb_err_iarg_not_gteia_ii_
      int_err(1) = 3 
      int_err(2) = ia
      int_err(3) = 5
      int_err(4) = psb_m_
      int_err(5) = desc_dec%get_global_rows()
    else if (desc_dec%get_global_cols() < ja) then
      info=psb_err_iarg_not_gteia_ii_
      int_err(1) = 4 
      int_err(2) = ja
      int_err(3) = 5
      int_err(4) = psb_n_
      int_err(5) = desc_dec%get_global_cols()
    else if (desc_dec%get_global_rows() < (ia+m-1)) then
      info=psb_err_iarg2_neg_
      int_err(1) = 1
      int_err(2) = m
      int_err(3) = 3
      int_err(4) = ia
    else if (desc_dec%get_global_cols() < (ja+n-1)) then
      info=psb_err_iarg2_neg_
      int_err(1) = 2
      int_err(2) = n
      int_err(3) = 4
      int_err(4) = ja
    end if

    if (info /= psb_success_) then
      call psb_errpush(info,name,i_err=int_err)
      goto 9999
    end if

    ! Compute local indices for submatrix starting
    ! at global indices ix and jx
    if(present(iia).and.present(jja)) then
      if (desc_dec%get_local_rows() > 0) then
        iia=1
        jja=1
      else
        iia=desc_dec%get_local_rows()+1
        jja=desc_dec%get_local_cols()+1
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
