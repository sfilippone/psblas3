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
  !  Purpose
  !  == = ====
  !  
  !  Allocate special descriptor for replicated index space.
  ! 
  !
  !
  ! INPUT
  ! == ====
  ! M                 :(Global Input) Integer 
  !                    Total number of  equations
  !                    required.
  !
  ! ictxt      : (Global Input)Integer BLACS context for an NPx1 grid 
  !                required.
  !
  ! OUTPUT
  ! == =======
  ! desc   : TYPEDESC
  ! desc OUTPUT FIELDS:
  !
  ! MATRIX_DATA   : Pointer to integer Array 
  !                contains some
  !               local and global information about matrix:
  !
  !  NOTATION        STORED IN		     EXPLANATION
  !  ------------ ---------------------- -------------------------------------
  !  DEC_TYPE        MATRIX_DATA[DEC_TYPE_]   Decomposition type, temporarly is
  !                      setted to 1( matrix not yet assembled)
  !  M 	             MATRIX_DATA[M_]          Total number of equations
  !  N 	             MATRIX_DATA[N_]          Total number of variables
  !  N_ROW           MATRIX_DATA[N_ROW_]      Number of local equations
  !  N_COL           MATRIX_DATA[N_COL_]      Number of local columns (see below)
  !  CTXT_A          MATRIX_DATA[CTXT_]     The BLACS context handle, 
  !                                         indicating
  !	  			            the global context of the operation
  !					    on the matrix.
  !					    The context itself is global.
  !
  !  GLOB_TO_LOC     Array of dimension equal to number of global 
  !                  rows/cols (MATRIX_DATA[M_]). On exit,
  !                  for all global indices either:
  !                  1. The index belongs to the current process; the entry
  !                     is set to the next free local row index.
  !                  2. The index belongs to process P (0<=P<=NP-1); the entry 
  !                     is set to 
  !                     -(NP+P+1)
  !
  !  LOC_TO_GLOB     An array of dimension equal to number of local cols N_COL
  !                  i.e. all columns of the matrix such that there is at least
  !                  one nonzero entry within the local row range. At the time 
  !                  this routine is called N_COL cannot be know, so we set 
  !                  N_COL=N_ROW, and dimension this vector on N_ROW plus an 
  !                  estimate. On exit the vector elements are set
  !                  to the index of the corresponding entry in GLOB_TO_LOC, or  
  !                  to -1 for indices I>N_ROW.
  !
  !
  !  HALO_INDEX      Not touched here, as it depends on the matrix pattern
  !
  !  OVRLAP_INDEX    On exit from this routine, the overlap indices are stored in
  !                  triples (Proc, 1, Index), similar to the assembled format 
  !                  but neither optimized, nor deadlock free. 
  !                  List is terminated with -1
  !
  !  OVRLAP_ELEM     On exit from this routine, just a list of pairs (index,#p).
  !                  List is terminated with -1.
  !                  
  !
  ! END OF desc OUTPUT FIELDS
  !
  !
subroutine psb_cdrep(m, ictxt, desc, info)
  use psb_base_mod
  use psi_mod
  use psb_repl_map_mod
  implicit None
  !....Parameters...
  Integer, intent(in)               :: m,ictxt
  integer, intent(out)              :: info
  Type(psb_desc_type), intent(out)  :: desc

  !locals
  Integer             :: i,np,me,err,n,err_act
  integer             :: int_err(5),exch(2), thalo(1), tovr(1), text(1)
  integer              :: debug_level, debug_unit
  character(len=20)   :: name

  if(psb_get_errstatus() /= 0) return 
  info=psb_success_
  err=0
  name = 'psb_cdrep'
  debug_unit  = psb_get_debug_unit()
  debug_level = psb_get_debug_level()

  call psb_info(ictxt, me, np)
  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': ',np

  n = m
  !... check m and n parameters....
  if (m < 1) then
    info = psb_err_iarg_neg_
    int_err(1) = 1
    int_err(2) = m
  else if (n < 1) then
    info = psb_err_iarg_neg_
    int_err(1) = 2
    int_err(2) = n
  endif

  if (info /= psb_success_) then 
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if

  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),':  doing global checks'  
  !global check on m and n parameters
  if (me == psb_root_) then
    exch(1)=m
    exch(2)=n
    call psb_bcast(ictxt,exch(1:2),root=psb_root_)
  else
    call psb_bcast(ictxt,exch(1:2),root=psb_root_)
    if (exch(1) /= m) then
      info=psb_err_parm_differs_among_procs_
      int_err(1)=1
    else if (exch(2) /= n) then
      info=psb_err_parm_differs_among_procs_
      int_err(1)=2
    endif
  endif

  if (info /= psb_success_) then 
    call psb_errpush(info,name,i_err=int_err)
    goto 9999
  end if


  call psb_nullify_desc(desc)


  !count local rows number
  ! allocate work vector
!!$  allocate(desc%matrix_data(psb_mdata_size_),&
!!$       &   desc%ovrlap_elem(0,3),stat=info)
!!$  if (info /= psb_success_) then     
!!$    info=psb_err_alloc_request_
!!$    int_err(1)=2*m+psb_mdata_size_+1
!!$    call psb_errpush(info,name,i_err=int_err,a_err='integer')
!!$    goto 9999
!!$  endif
!!$  ! If the index space is replicated there's no point in not having 
!!$  ! the full map on the current process. 
!!$
!!$  desc%matrix_data(psb_m_)        = m
!!$  desc%matrix_data(psb_n_)        = n
!!$  desc%matrix_data(psb_n_row_)    = m
!!$  desc%matrix_data(psb_n_col_)    = n
!!$  desc%matrix_data(psb_ctxt_)     = ictxt
!!$  call psb_get_mpicomm(ictxt,desc%matrix_data(psb_mpi_c_))
!!$  desc%matrix_data(psb_dec_type_) = psb_desc_bld_


  allocate(psb_repl_map :: desc%indxmap, stat=info)
  select type(aa => desc%indxmap) 
  type is (psb_repl_map) 
    call aa%repl_map_init(ictxt,m,info)
  class default 
    ! This cannot happen 
    info = psb_err_internal_error_
    call psb_errpush(info,name)
    Goto 9999
  end select
   

  tovr  = -1 
  call psi_bld_tmpovrl(tovr,desc,info)
!!$  desc%matrix_data(psb_dec_type_) = psb_desc_bld_
  

  if (debug_level >= psb_debug_ext_) &
       & write(debug_unit,*) me,' ',trim(name),': end'

  call psb_erractionrestore(err_act)
  return

9999 continue
  call psb_erractionrestore(err_act)
  if (err_act == psb_act_abort_) then
    call psb_error(ictxt)
    return
  end if
  return

end subroutine psb_cdrep
