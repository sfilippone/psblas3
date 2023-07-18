#include "../include/Sp3MM_CSR_OMP_Multi.h"
#include "../include/utils.h"
#include "../include/ompChunksDivide.h"
#include <stdio.h>

enum impl_types
{
    ROW_BY_ROW_UB
};

CHUNKS_DISTR chunksFair;
CHUNKS_DISTR chunksFairFolded;
CHUNKS_DISTR chunksNOOP;

CHUNKS_DISTR_INTERF chunkDistrbFunc=&chunksFairFolded;

static CONFIG Conf = {
	.gridRows  = 20,
	.gridCols  = 2,
	.symbMMRowImplID = RBTREE,
};

void psb_f_spmm_build_spacc(idx_t a_m, idx_t a_n, idx_t a_nz,
                             double **a_as_ptr, idx_t **a_ja_ptr,
                             idx_t **a_irp_ptr,
                             idx_t b_m, idx_t b_n, idx_t b_nz,
                             double **b_as_ptr, idx_t **b_ja_ptr,
                             idx_t **b_irp_ptr,
                             enum impl_types impl_choice,
                             void **accumul,
                             void **rows_sizes,
                             void **tmp_matrix,
                             idx_t *nnz)
{
    int rc;
    spmat a, b;

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

	if (!getConfig(&Conf)){
		VERBOSE printf("configuration changed from env");
	}
	//save exported parallelization grid configured ... see Config.md in the top folder
	//const ushort GRIDROWS = Conf.gridRows;
	//const ushort GRIDCOLS = Conf.gridCols;
	const ushort ORIG_GRID_CONF[2] = {Conf.gridRows,Conf.gridCols};

	int maxThreads = omp_get_max_threads();
	Conf.threadNum = (uint) maxThreads;
	DEBUG   printf("omp_get_max_threads:\t%d\n",maxThreads); 

	/*
	 * get exported schedule configuration, 
	 * if schedule != static -> dyn like -> set a chunk division function before omp for
	 */
	int schedKind_chunk_monotonic[3];
	ompGetRuntimeSchedule(schedKind_chunk_monotonic);
	Conf.chunkDistrbFunc = &chunksNOOP; 
	if (schedKind_chunk_monotonic[0] != omp_sched_static)
		Conf.chunkDistrbFunc = chunkDistrbFunc;
	VERBOSE 
	  printf("%s",Conf.chunkDistrbFunc == &chunksNOOP?"static schedule =>chunkDistrbFunc NOOP\n":"");

    int i, row_nz;

    // setting up spmat type matrices
    a.M = a_m;
    a.N = a_n;
    a.NZ = a_nz;
    a.AS = a_as;
    a.JA = a_ja;
    a.IRP = a_irp;
    a.MAX_ROW_NZ = a.IRP[1] - a.IRP[0];
    for (i = 1; i < a.M; i ++){
        row_nz = a.IRP[i+1] - a.IRP[i];
        if (row_nz > a.MAX_ROW_NZ){
            a.MAX_ROW_NZ = row_nz;
        }
    }

    b.M = b_m;
    b.N = b_n;
    b.NZ = b_nz;
    b.AS = b_as;
    b.JA = b_ja;
    b.IRP = b_irp;
    b.MAX_ROW_NZ = b.IRP[1] - b.IRP[0];
    for (i = 1; i < b.M; i ++){
        row_nz = b.IRP[i+1] - b.IRP[i];
        if (row_nz > b.MAX_ROW_NZ){
            b.MAX_ROW_NZ = row_nz;
        }
    }


    // computing the size
    switch (impl_choice)
    {
    case ROW_BY_ROW_UB:
        *nnz = spmmRowByRowCalculateSize_1(&a, &b, &Conf, accumul, rows_sizes, tmp_matrix);
    default:
        break;
    }
}

void psb_f_spmm_merge_spacc(void **accumul,
                            void **rows_sizes,
                            void **tmp_matrix,
                            enum impl_types impl_choice,
                            double **c_as_ptr,
                            idx_t **c_ja_ptr,
                            idx_t **c_irp_ptr,
                            int *info)
{
    double *c_as = *c_as_ptr;
    idx_t *c_ja = *c_ja_ptr;
    idx_t *c_irp = *c_irp_ptr;
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