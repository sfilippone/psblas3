
#ifndef PSB_C_BASE__
#define PSB_C_BASE__
#ifdef __cplusplus
extern "C" {
  /*typedef char _Bool;*/
#endif

#include <float.h>
#include <complex.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>


  typedef int32_t psb_m_t;

#if defined(IPK4) &&  defined(LPK4)
  typedef int32_t psb_i_t;
  typedef int32_t psb_l_t;
#elif defined(IPK4) &&  defined(LPK8)
  typedef int32_t psb_i_t;
  typedef int64_t psb_l_t;
#elif defined(IPK8) &&  defined(LPK8)
  typedef int64_t psb_i_t;
  typedef int64_t psb_l_t;
#else
#endif
  typedef int64_t psb_e_t;

  typedef float  psb_s_t;
  typedef double psb_d_t;
  typedef float  complex psb_c_t;
  typedef double complex psb_z_t;
#define PSB_ERR_ERROR  -1
#define PSB_ERR_SUCCESS 0


  typedef struct PSB_C_DESCRIPTOR {
    void *descriptor;
  } psb_c_descriptor;



  psb_i_t psb_c_error();
  psb_i_t psb_c_clean_errstack();
  void psb_c_print_errmsg();
  char *psb_c_pop_errmsg();
  psb_i_t  psb_c_f2c_errmsg(char *, psb_i_t);
  void psb_c_seterraction_ret();
  void psb_c_seterraction_print();
  void psb_c_seterraction_abort();

  /* Environment routines */
  psb_i_t psb_c_init();
  void    psb_c_exit_ctxt(psb_i_t ictxt);
  void    psb_c_exit(psb_i_t ictxt);
  void    psb_c_abort(psb_i_t ictxt);
  void    psb_c_barrier(psb_i_t ictxt);
  void    psb_c_info(psb_i_t ictxt, psb_i_t *iam, psb_i_t *np);
  psb_d_t  psb_c_wtime();
  psb_i_t psb_c_get_errstatus();

  psb_i_t psb_c_get_index_base();
  void psb_c_set_index_base(psb_i_t base);

  void   psb_c_mbcast(psb_i_t ictxt, psb_i_t n, psb_m_t *v, psb_i_t root);
  void   psb_c_ibcast(psb_i_t ictxt, psb_i_t n, psb_i_t *v, psb_i_t root);
  void   psb_c_lbcast(psb_i_t ictxt, psb_i_t n, psb_l_t *v, psb_i_t root);
  void   psb_c_ebcast(psb_i_t ictxt, psb_i_t n, psb_e_t *v, psb_i_t root);
  void   psb_c_sbcast(psb_i_t ictxt, psb_i_t n, psb_s_t *v, psb_i_t root);
  void   psb_c_dbcast(psb_i_t ictxt, psb_i_t n, psb_d_t *v, psb_i_t root);
  void   psb_c_cbcast(psb_i_t ictxt, psb_i_t n, psb_c_t *v, psb_i_t root);
  void   psb_c_zbcast(psb_i_t ictxt, psb_i_t n, psb_z_t *v, psb_i_t root);
  void   psb_c_hbcast(psb_i_t ictxt, const char *v, psb_i_t root);

  /* Descriptor/integer routines */
  psb_c_descriptor* psb_c_new_descriptor();
  psb_i_t    psb_c_cdall_vg(psb_l_t ng, psb_i_t *vg, psb_i_t ictxt, psb_c_descriptor *cd);
  psb_i_t    psb_c_cdall_vl(psb_i_t nl, psb_l_t *vl, psb_i_t ictxt, psb_c_descriptor *cd);
  psb_i_t    psb_c_cdall_nl(psb_i_t nl, psb_i_t ictxt, psb_c_descriptor *cd);
  psb_i_t    psb_c_cdall_repl(psb_l_t n, psb_i_t ictxt, psb_c_descriptor *cd);
  psb_i_t    psb_c_cdasb(psb_c_descriptor *cd);
  psb_i_t    psb_c_cdfree(psb_c_descriptor *cd);
  psb_i_t    psb_c_cdins(psb_i_t nz, const psb_l_t *ia, const psb_l_t *ja, psb_c_descriptor *cd);


  psb_i_t    psb_c_cd_get_local_rows(psb_c_descriptor *cd);
  psb_i_t    psb_c_cd_get_local_cols(psb_c_descriptor *cd);
  psb_l_t    psb_c_cd_get_global_rows(psb_c_descriptor *cd);
  psb_i_t    psb_c_cd_get_global_indices(psb_l_t idx[], psb_i_t nidx, bool owned, psb_c_descriptor *cd);

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
