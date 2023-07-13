#include "../include/Sp3MM_CSR_OMP_Multi.h"
#include "../include/utils.h"
#include "../include/ompChunksDivide.h"
#include <stdio.h>

enum impl_types
{
    ROW_BY_ROW_UB
};

CHUNKS_DISTR chunksFair;

void psb_f_spmm_build_spacc(idx_t a_m, idx_t a_n, idx_t a_nz,
                             double **a_as_ptr, idx_t **a_ja_ptr,
                             idx_t **a_irp_ptr, idx_t a_max_row_nz,
                             idx_t b_m, idx_t b_n, idx_t b_nz,
                             double **b_as_ptr, idx_t **b_ja_ptr,
                             idx_t **b_irp_ptr, idx_t b_max_row_nz,
                             enum impl_types impl_choice,
                             void **accumul,
                             void **rows_sizes,
                             void **tmp_matrix,
                             idx_t *nnz)
{
    int rc;
    spmat a, b;
    CONFIG cfg;

    double *a_as = *a_as_ptr;
    idx_t *a_ja = *a_ja_ptr;
    idx_t *a_irp = *a_irp_ptr;

    double *b_as = *b_as_ptr;
    idx_t *b_ja = *b_ja_ptr;
    idx_t *b_irp = *b_irp_ptr;
    
#ifdef ROWLENS
    a->RL = a_rl;
    b->RL = b_rl;
#endif // ROWLENS

    // setting up cfg
    // TODO : CHECK THAT THIS IS COMPATIBLE WITH PSB
    rc = getConfig(&cfg);

    // TODO : change chunk distribution with a choice ?
    cfg.chunkDistrbFunc = &chunksFair;
    cfg.threadNum = 1;

    printf("irp %d %d %d %d %d\n", a_irp[0], a_irp[1], a_irp[2], a_irp[3], a_irp[4]);
    // setting up spmat type matrices
    a.M = a_m;
    a.N = a_n;
    a.NZ = a_nz;
    a.AS = a_as;
    a.JA = a_ja;
    a.IRP = a_irp;
    a.MAX_ROW_NZ = a_max_row_nz;

    b.M = b_m;
    b.N = b_n;
    b.NZ = b_nz;
    b.AS = b_as;
    b.JA = b_ja;
    b.IRP = b_irp;
    b.MAX_ROW_NZ = b_max_row_nz;

    // computing the size
    switch (impl_choice)
    {
    case ROW_BY_ROW_UB:
        *nnz = spmmRowByRowCalculateSize_1(&a, &b, &cfg, accumul, rows_sizes, tmp_matrix);
    default:
        break;
    }
}

void psb_f_spmm_merge_spacc(void **accumul,
                            void **rows_sizes,
                            void **tmp_matrix,
                            enum impl_types impl_choice,
                            double **c_as,
                            idx_t **c_ja,
                            idx_t **c_irp,
                            int *info)
{
    // merging the rows into the correct arrays
    switch (impl_choice)
    {
    case ROW_BY_ROW_UB:
        spmmRowByRowPopulate_1(accumul, rows_sizes, tmp_matrix, c_as, c_ja, c_irp);
        break;
    default:
        break;
    }
}