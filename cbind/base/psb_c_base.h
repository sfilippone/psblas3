#ifndef PSB_C_BASE__
#define PSB_C_BASE__
#ifdef __cplusplus
extern "C" {
  typedef char _Bool;
#endif
  
  typedef int psb_err_t;
  typedef int psb_ctx_t;
#define PSB_ERR_ERROR  -1
#define PSB_ERR_SUCCESS 0
  
  
  typedef struct PSB_C_DESCRIPTOR {
    void *descriptor;
  } psb_c_descriptor; 
  
  

  int psb_c_error();
  int psb_c_clean_errstack();
  void psb_c_print_errmsg();
  char *psb_c_pop_errmsg();
  int  psb_c_f2c_errmsg(char *, int);
  void psb_c_seterraction_ret();
  void psb_c_seterraction_print();
  void psb_c_seterraction_abort();

  /* Environment routines */ 
  int    psb_c_init();
  void   psb_c_exit_ctxt(int ictxt);
  void   psb_c_exit(int ictxt);
  void   psb_c_abort(int ictxt);
  void   psb_c_barrier(int ictxt);
  void   psb_c_info(int ictxt, int *iam, int *np);
  double psb_c_wtime();
  int    psb_c_get_errstatus();

  void   psb_c_ibcast(int ictxt, int n, int *v, int root);
  void   psb_c_dbcast(int ictxt, int n, double *v, int root);
  void   psb_c_hbcast(int ictxt, const char *v, int root);
  
  /* Descriptor/integer routines */ 
  psb_c_descriptor* psb_c_new_descriptor();
  int    psb_c_cdall_vg(int ng, int *vg, int ictxt, psb_c_descriptor *cd);
  int    psb_c_cdall_vl(int nl, int *vl, int ictxt, psb_c_descriptor *cd);
  int    psb_c_cdall_nl(int nl, int ictxt, psb_c_descriptor *cd);
  int    psb_c_cdall_repl(int n, int ictxt, psb_c_descriptor *cd);
  int    psb_c_cdasb(psb_c_descriptor *cd);
  int    psb_c_cdfree(psb_c_descriptor *cd);
  int    psb_c_cdins(int nz, const int *ia, const int *ja, psb_c_descriptor *cd);
  

  int    psb_c_cd_get_local_rows(psb_c_descriptor *cd);
  int    psb_c_cd_get_local_cols(psb_c_descriptor *cd);
  int    psb_c_cd_get_global_rows(psb_c_descriptor *cd);
  int    psb_c_cd_get_global_rows(psb_c_descriptor *cd);


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

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif
