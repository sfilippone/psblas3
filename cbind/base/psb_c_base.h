#ifndef PSB_C_BASE__
#define PSB_C_BASE__
#ifdef __cplusplus
extern "C" {
  /*typedef char _Bool;*/
#endif

#include <float.h>
#ifdef __cplusplus
#include <complex>
#else
#include <complex.h>
#endif
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "psb_config.h"
#include "psb_types.h"


  typedef struct PSB_C_DESCRIPTOR {
    void *descriptor;
  } psb_c_descriptor;


  typedef struct PSB_C_CTXT {
    psb_i_t *ctxt;
  } psb_c_ctxt;



  void  psb_c_check_error(psb_c_ctxt cctxt);
  psb_i_t psb_c_error();
  psb_i_t psb_c_clean_errstack();
  void psb_c_print_errmsg();
  char *psb_c_pop_errmsg();
  psb_i_t  psb_c_f2c_errmsg(char *, psb_i_t);
  void psb_c_seterraction_ret();
  void psb_c_seterraction_print();
  void psb_c_seterraction_abort();

  /* Environment routines */
  void    psb_c_init(psb_c_ctxt *cctxt);
  void    psb_c_init_from_fint(psb_c_ctxt *cctxt, psb_i_t f_comm);
  void    psb_c_exit(psb_c_ctxt cctxt);
  void    psb_c_exit_ctxt(psb_c_ctxt cctxt);
  void    psb_c_abort(psb_c_ctxt cctxt);
  void    psb_c_barrier(psb_c_ctxt cctxt);
  void    psb_c_info(psb_c_ctxt cctxt, psb_i_t *iam, psb_i_t *np);
  void    psb_c_get_i_ctxt(psb_c_ctxt cctxt, psb_i_t *ictxt, psb_i_t *info);
  bool    psb_c_cmp_ctxt(psb_c_ctxt cctxt1, psb_c_ctxt cctxt2);
  psb_d_t  psb_c_wtime();
  psb_i_t psb_c_get_errstatus();

  psb_i_t psb_c_get_index_base();
  void psb_c_set_index_base(psb_i_t base);

  void   psb_c_mbcast(psb_c_ctxt cctxt, psb_i_t n, psb_m_t *v, psb_i_t root);
  void   psb_c_ibcast(psb_c_ctxt cctxt, psb_i_t n, psb_i_t *v, psb_i_t root);
  void   psb_c_lbcast(psb_c_ctxt cctxt, psb_i_t n, psb_l_t *v, psb_i_t root);
  void   psb_c_ebcast(psb_c_ctxt cctxt, psb_i_t n, psb_e_t *v, psb_i_t root);
  void   psb_c_sbcast(psb_c_ctxt cctxt, psb_i_t n, psb_s_t *v, psb_i_t root);
  void   psb_c_dbcast(psb_c_ctxt cctxt, psb_i_t n, psb_d_t *v, psb_i_t root);
  void   psb_c_cbcast(psb_c_ctxt cctxt, psb_i_t n, psb_c_t *v, psb_i_t root);
  void   psb_c_zbcast(psb_c_ctxt cctxt, psb_i_t n, psb_z_t *v, psb_i_t root);
  void   psb_c_hbcast(psb_c_ctxt cctxt, const char *v, psb_i_t root);

  /* Descriptor/integer routines */
  psb_c_descriptor* psb_c_new_descriptor();
  void psb_c_delete_descriptor(psb_c_descriptor *);
  psb_c_ctxt* psb_c_new_ctxt();
  void psb_c_delete_ctxt(psb_c_ctxt *);
  psb_i_t    psb_c_cdall_vg(psb_l_t ng, psb_i_t *vg, psb_c_ctxt cctxt, psb_c_descriptor *cd);
  psb_i_t    psb_c_cdall_vl(psb_i_t nl, psb_l_t *vl, psb_c_ctxt cctxt, psb_c_descriptor *cd);
  psb_i_t    psb_c_cdall_nl(psb_i_t nl, psb_c_ctxt cctxt, psb_c_descriptor *cd);
  psb_i_t    psb_c_cdall_repl(psb_l_t n, psb_c_ctxt cctxt, psb_c_descriptor *cd);
  psb_i_t    psb_c_cdasb(psb_c_descriptor *cd);
  psb_i_t    psb_c_cdfree(psb_c_descriptor *cd);
  psb_i_t    psb_c_cdins(psb_i_t nz, const psb_l_t *ia, const psb_l_t *ja, psb_c_descriptor *cd);
  bool       psb_c_is_owned(psb_l_t gindex, psb_c_descriptor *cd);

  psb_i_t    psb_c_cd_get_local_rows(psb_c_descriptor *cd);
  psb_i_t    psb_c_cd_get_local_cols(psb_c_descriptor *cd);
  psb_l_t    psb_c_cd_get_global_rows(psb_c_descriptor *cd);
  psb_l_t    psb_c_cd_get_global_cols(psb_c_descriptor *cd);
  psb_i_t    psb_c_cd_get_global_indices(psb_l_t idx[], psb_i_t nidx, bool owned, psb_c_descriptor *cd);
  psb_i_t    psb_c_g2l(psb_c_descriptor *cdh,psb_l_t gindex,bool cowned);


  /*  legal values for upd argument */
#define psb_upd_srch_   98764
#define psb_upd_perm_   98765
#define psb_upd_def_   psb_upd_srch_
  /*  legal values for dupl argument */
#define psb_dupl_ovwrt_  0
#define psb_dupl_add_    1
#define psb_dupl_err_    2
#define psb_dupl_def_    psb_dupl_ovwrt_

  /* legal values for afmt */
#define PSB_AFMT_CSR     "CSR"
#define PSB_AFMT_CSC     "CSC"
#define PSB_AFMT_COO     "COO"
#define PSB_AFMT_RSB     "RSB"

  /* Transpose argument */
#define psb_NoTrans_    "N"
#define psb_Trans_      "T"
#define psb_ConjTrans_  "C"

  /*  legal values for halo swap modes argument */
#define  psb_swap_send_  1
#define  psb_swap_recv_  2
#define  psb_swap_sync_  4
#define  psb_swap_mpi_   8

  /*  legal values for ovrl update argument */
#define psb_none_        0
#define psb_sum_         1
#define psb_avg_         2
#define psb_square_root_ 3
#define psb_setzero_     4


#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
