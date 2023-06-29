#include "../include/Sp3MM_CSR_OMP_Multi.h"
#include "../include/utils.h"

/**
 * @brief performs multiplication of two sparse matrices 
 * A and B stored in the CRS format. The resulting matrix is C
 * 
 * @param[in] a_m           number of columns of A
 * @param[in] a_n           number of rows of A
 * @param[in] a_nz          number of non zero elements of A
 * @param[in] a_as          array with the non zero coefficients of A
 * @param[in] a_ja          array with the column indices of the non 
 *                      zero coefficients of A
 * @param[in] a_irp         array with the indices of the beginning of 
 *                      each row in a_as and a_ja
 * @param[in] a_rl          array with the lengths of each rows of A
 * @param[in] a_max_row_nz  maximum number of coefficients in a row in A
 * @param[in] b_m           number of columns of B
 * @param[in] b_n           number of rows of B
 * @param[in] b_nz          number of non zero elements of B
 * @param[in] b_as          array with the non zero coefficients of B
 * @param[in] b_ja          array with the column indices of the non
 *                      zero coefficients of B
 * @param[in] b_irp         array with the indices of the beginning of
 *                      each row in b_as and b_ja
 * @param[in] b_rl          array with the lengths of each rows of B
 * @param[in] b_max_row_nz  maximum number of coefficients in a row in B
 * @param[out] c_m           number of colums of C
 * @param[out] c_n           number of rows of C
 * @param[out] c_nz          number of non zero elements in C
 * @param[out] c_as          array with the non zero coefficients of C
 * @param[out] c_ja          array with the column indices of the non
 *                      zero coefficients of C
 * @param[out] c_irp         array with the indices of the beginning of
 *                      each row in c_as and c_ja
 * @param[out] c_rl          array with the lengths of each rows of C
 * @param[out] c_max_row_nz  maximum number of coefficients in a row in C
 * @param[out] info          return value to check if the operation was successful
 */
#ifdef ROWLENS
void psb_f_spmm_row_by_row_ub_0(idx_t a_m, idx_t a_n, idx_t a_nz, 
                                double *a_as, idx_t *a_ja, 
                                idx_t *a_irp, idx_t *a_rl, idx_t a_max_row_nz, 
                                idx_t b_m, idx_t b_n, idx_t b_nz, 
                                double *b_as, idx_t *b_ja, 
                                idx_t *b_irp, idx_t *b_rl, idx_t b_max_row_nz,
                                idx_t *c_m, idx_t *c_n, idx_t *c_nz, 
                                double **c_as, idx_t **c_ja, 
                                idx_t **c_irp, idx_t **c_rl, idx_t *c_max_row_nz,
                                int info)
#else
void psb_f_spmm_row_by_row_ub_0(idx_t a_m, idx_t a_n, idx_t a_nz, 
                                double *a_as, idx_t *a_ja, 
                                idx_t *a_irp, idx_t a_max_row_nz, 
                                idx_t b_m, idx_t b_n, idx_t b_nz, 
                                double *b_as, idx_t *b_ja, 
                                idx_t *b_irp, idx_t b_max_row_nz,
                                idx_t *c_m, idx_t *c_n, idx_t *c_nz, 
                                double **c_as, idx_t **c_ja, 
                                idx_t **c_irp,  idx_t *c_max_row_nz,
                                int info)
#endif                            
{
    int rc;
    spmat *a, *b, *c;
    CONFIG *cfg;

    #ifdef ROWLENS
    a->RL = a_rl;
    b->RL = b_rl;
    #endif//ROWLENS

    // setting up cfg
    // TODO : CHECK THAT THIS IS COMPATIBLE WITH PSB
    rc = getConfig(cfg);

    // setting up spmat type matrices
    a->M = a_m;
    a->N = a_n;
    a->NZ = a_nz;
    a->AS = a_as;
    a->JA = a_ja;
    a->IRP = a_irp;
    a->MAX_ROW_NZ = a_max_row_nz;

    b->M = b_m;
    b->N = b_n;
    b->NZ = b_nz;
    b->AS = b_as;
    b->JA = b_ja;
    b->IRP = b_irp;
    b->MAX_ROW_NZ = b_max_row_nz;

    // performing spmm
    c = spmmRowByRow_0(a, b, cfg);

    // output result
    *(c_m) = c->M;
    *(c_n) = c->N;
    *(c_nz)= c->NZ;
    *(c_as)= c->AS;
    *(c_ja)= c->JA;
    *(c_irp)=c->IRP;
    #ifdef ROWLENS
    *(c_rl)= c->RL;
    #endif
}