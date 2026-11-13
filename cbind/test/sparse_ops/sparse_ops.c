/* Distributed rectangular sparse-operation C binding regression.
 * Compile the same test for each scalar family; no internal Fortran objects
 * are accessed.  The reference products are computed as small dense matrices.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "psb_base_cbind.h"

#if defined(TEST_SINGLE_REAL)
#define SCALAR_TAG s
#define SCALAR_NAME "s"
typedef psb_c_s_t scalar_t;
#define SCALAR_ABS(x) fabs((double)(x))
#define SCALAR_CONJ(x) (x)
#define SCALAR_VALUE(re, im) ((scalar_t)(re))
#define TEST_TOL 5.e-5
#elif defined(TEST_SINGLE_COMPLEX)
#define SCALAR_TAG c
#define SCALAR_NAME "c"
typedef psb_c_c_t scalar_t;
#define SCALAR_ABS(x) ((double)cabsf(x))
#define SCALAR_CONJ(x) conjf(x)
#define SCALAR_VALUE(re, im) ((scalar_t)((re) + I * (im)))
#define TEST_TOL 5.e-5
#elif defined(TEST_DOUBLE_COMPLEX)
#define SCALAR_TAG z
#define SCALAR_NAME "z"
typedef psb_c_z_t scalar_t;
#define SCALAR_ABS(x) cabs(x)
#define SCALAR_CONJ(x) conj(x)
#define SCALAR_VALUE(re, im) ((scalar_t)((re) + I * (im)))
#define TEST_TOL 2.e-12
#else
#define SCALAR_TAG d
#define SCALAR_NAME "d"
typedef psb_c_d_t scalar_t;
#define SCALAR_ABS(x) fabs(x)
#define SCALAR_CONJ(x) (x)
#define SCALAR_VALUE(re, im) ((scalar_t)(re))
#define TEST_TOL 2.e-12
#endif

#define PSB_NAME_(tag, name) psb_c_##tag##name
#define PSB_NAME_EXPAND(tag, name) PSB_NAME_(tag, name)
#define PSB_NAME(name) PSB_NAME_EXPAND(SCALAR_TAG, name)
typedef PSB_NAME(spmat) matrix_t;
typedef PSB_NAME(vector) vector_t;

enum { PRESSURE_SIZE = 3, VELOCITY_SIZE = 5, MAX_ENTRIES = 64 };

static int rank_id;
static int process_count;
static int index_base;
static psb_c_ctxt context = {NULL};

static void require(int condition, const char *what)
{
  int failed = !condition;
  int any_failed = 0;
  if (failed)
    fprintf(stderr, "[%d] %s (%s, index base %d)\n",
            rank_id, what, SCALAR_NAME, index_base);
  MPI_Allreduce(&failed, &any_failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (any_failed) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

static void check(psb_c_i_t status, const char *what)
{
  if (status != 0)
    fprintf(stderr, "[%d] %s returned %lld\n", rank_id, what,
            (long long)status);
  require(status == 0, what);
}

#define CHECK(call) check((call), #call)

static void expect_error(psb_c_i_t status, const char *what)
{
  require(status != 0, what);
  psb_c_clean_errstack();
}

static int local_count(int global_size)
{
  return global_size / process_count +
         (rank_id < global_size % process_count);
}

static int first_owned(int global_size)
{
  int remainder = global_size % process_count;
  return rank_id * (global_size / process_count) +
         (rank_id < remainder ? rank_id : remainder);
}

/* Exercise both contiguous and explicit-index allocation, including nl=0. */
static void make_descriptor(int n, int explicit_indices,
                            psb_c_descriptor *descriptor)
{
  psb_c_l_t owned[VELOCITY_SIZE] = {0};
  const int nl = local_count(n);
  const int first = first_owned(n);
  for (int i = 0; i < nl; ++i) owned[i] = first + i + index_base;
  if (explicit_indices)
    CHECK(psb_c_cdall_vl(nl, owned, context, descriptor));
  else
    CHECK(psb_c_cdall_nl(nl, context, descriptor));
  require(psb_c_cd_get_local_rows(descriptor) == nl,
          "wrong owned row count");
  require(psb_c_cd_get_global_rows(descriptor) == n,
          "wrong global row count");
}

static void assemble(const scalar_t *dense, int nr, int nc,
                     matrix_t *a, psb_c_descriptor *rows,
                     psb_c_descriptor *columns, int explicit_rows)
{
  psb_c_l_t ir[MAX_ENTRIES] = {0}, jc[MAX_ENTRIES] = {0};
  scalar_t value[MAX_ENTRIES] = {0};
  psb_c_i_t nz = 0;
  make_descriptor(nr, explicit_rows, rows);
  make_descriptor(nc, !explicit_rows, columns);
  CHECK(PSB_NAME(spall_2desc)(a, rows, columns));
  for (int i = first_owned(nr); i < first_owned(nr) + local_count(nr); ++i)
    for (int j = 0; j < nc; ++j)
      if (SCALAR_ABS(dense[i * nc + j]) != 0.) {
        ir[nz] = i + index_base;
        jc[nz] = j + index_base;
        value[nz++] = dense[i * nc + j];
      }
  CHECK(PSB_NAME(spins_2desc)(nz, ir, jc, value, a, rows, columns));
  CHECK(psb_c_cdasb(rows));
  CHECK(psb_c_cdasb(columns));
  CHECK(PSB_NAME(spasb_2desc)(a, rows, columns));
}

static void transpose_dense(const scalar_t *a, int nr, int nc,
                            scalar_t *at, int conjugate)
{
  for (int i = 0; i < nr; ++i)
    for (int j = 0; j < nc; ++j)
      at[j * nr + i] = conjugate ? SCALAR_CONJ(a[i * nc + j])
                                  : a[i * nc + j];
}

static void multiply_dense(const scalar_t *a, const scalar_t *b,
                           int nr, int inner, int nc, scalar_t *c)
{
  for (int i = 0; i < nr; ++i)
    for (int j = 0; j < nc; ++j) {
      c[i * nc + j] = 0;
      for (int k = 0; k < inner; ++k)
        c[i * nc + j] += a[i * inner + k] * b[k * nc + j];
    }
}

/* Check exported global indices and every owned dense entry, including zeros. */
static void check_matrix(matrix_t *a, psb_c_descriptor *rows,
                         psb_c_descriptor *columns, int nr, int nc,
                         const scalar_t *expected, const char *what)
{
  psb_c_i_t nz = -1, returned_nz = -1;
  psb_c_l_t ir[MAX_ENTRIES] = {0}, jc[MAX_ENTRIES] = {0};
  scalar_t value[MAX_ENTRIES] = {0};
  scalar_t actual[MAX_ENTRIES] = {0};
  int valid = 1;
  CHECK(PSB_NAME(spget_coo)(a, rows, columns, 0, &nz, ir, jc, value));
  require(nz >= 0 && nz <= MAX_ENTRIES, "invalid exported size");
  CHECK(PSB_NAME(spget_coo)(a, rows, columns, MAX_ENTRIES,
                          &returned_nz, ir, jc, value));
  require(returned_nz == nz, "export size changed");
  for (psb_c_i_t k = 0; k < nz; ++k) {
    const psb_c_l_t i = ir[k] - index_base;
    const psb_c_l_t j = jc[k] - index_base;
    if (i < 0 || i >= nr || j < 0 || j >= nc ||
        !psb_c_is_owned(ir[k], rows)) {
      valid = 0;
    } else {
      actual[i * nc + j] += value[k];
    }
  }
  for (int i = 0; i < nr; ++i)
    if (psb_c_is_owned(i + index_base, rows))
      for (int j = 0; j < nc; ++j) {
        double difference = SCALAR_ABS(actual[i * nc + j] - expected[i * nc + j]);
        if (!(difference <= TEST_TOL * (1. + SCALAR_ABS(expected[i * nc + j]))))
          valid = 0;
      }
  require(valid, what);
  require(PSB_NAME(mat_get_nrows)(a) == psb_c_cd_get_local_rows(rows),
          "matrix has unexpected local row size");
  require(PSB_NAME(mat_get_ncols)(a) == psb_c_cd_get_local_cols(columns),
          "matrix and column descriptor disagree");
}

struct descriptor_snapshot {
  psb_c_i_t rows, columns;
  psb_c_l_t global_rows;
  psb_c_l_t indices[MAX_ENTRIES];
};

static struct descriptor_snapshot snapshot(psb_c_descriptor *descriptor)
{
  struct descriptor_snapshot s;
  memset(&s, 0, sizeof(s));
  s.rows = psb_c_cd_get_local_rows(descriptor);
  s.columns = psb_c_cd_get_local_cols(descriptor);
  s.global_rows = psb_c_cd_get_global_rows(descriptor);
  require(s.columns >= 0 && s.columns <= MAX_ENTRIES,
          "unexpected descriptor size");
  CHECK(psb_c_cd_get_global_indices(s.indices, s.columns, false, descriptor));
  return s;
}

static void check_descriptor(psb_c_descriptor *descriptor,
                             const struct descriptor_snapshot *before)
{
  const struct descriptor_snapshot after = snapshot(descriptor);
  require(before->rows == after.rows && before->columns == after.columns &&
          before->global_rows == after.global_rows &&
          memcmp(before->indices, after.indices,
                 (size_t)after.columns * sizeof(after.indices[0])) == 0,
          "input descriptor was modified");
}

static void copy_and_scale(matrix_t *source, psb_c_descriptor *rows,
                           matrix_t *destination)
{
  psb_c_l_t global[VELOCITY_SIZE] = {0};
  scalar_t weights[VELOCITY_SIZE] = {0};
  const psb_c_i_t nl = psb_c_cd_get_local_rows(rows);
  CHECK(PSB_NAME(spall)(destination, rows));
  CHECK(PSB_NAME(copy_mat)(source, destination, rows));
  CHECK(psb_c_cd_get_global_indices(global, nl, true, rows));
  for (psb_c_i_t i = 0; i < nl; ++i)
    weights[i] = (scalar_t)(global[i] - index_base + 1);
  CHECK(PSB_NAME(sprowscl)(destination, nl, weights));
}

static void check_export_capacity(matrix_t *a, psb_c_descriptor *rows,
                                  psb_c_descriptor *columns)
{
  psb_c_i_t nz = 0, required = -1;
  psb_c_l_t ir[MAX_ENTRIES], jc[MAX_ENTRIES];
  scalar_t values[MAX_ENTRIES];
  int valid = 1;
  CHECK(PSB_NAME(spget_coo)(a, rows, columns, 0, &nz, ir, jc, values));
  for (int k = 0; k < MAX_ENTRIES; ++k) {
    ir[k] = -17;
    jc[k] = -29;
    values[k] = (scalar_t)-41;
  }
  /* Export is local; ranks with zero entries do not need a failure case. */
  if (nz > 1) {
    if (PSB_NAME(spget_coo)(a, rows, columns, nz - 1,
                           &required, ir, jc, values) == 0) valid = 0;
    if (required != nz) valid = 0;
    for (int k = 0; k < MAX_ENTRIES; ++k)
      if (ir[k] != -17 || jc[k] != -29 ||
          SCALAR_ABS(values[k] - (scalar_t)-41) != 0.) valid = 0;
  }
  require(valid, "short export buffer was accepted or overwritten");
  psb_c_clean_errstack();
}

static void set_vector(vector_t *x, psb_c_descriptor *descriptor,
                       const scalar_t *values, scalar_t halo_value)
{
  scalar_t local[MAX_ENTRIES];
  psb_c_l_t owned[MAX_ENTRIES] = {0};
  const psb_c_i_t nr = psb_c_cd_get_local_rows(descriptor);
  const psb_c_i_t nc = psb_c_cd_get_local_cols(descriptor);
  require(nc >= 0 && nc <= MAX_ENTRIES, "unexpected vector size");
  CHECK(psb_c_cd_get_global_indices(owned, nr, true, descriptor));
  for (psb_c_i_t i = 0; i < nc; ++i) local[i] = halo_value;
  for (psb_c_i_t i = 0; i < nr; ++i) local[i] = values[owned[i] - index_base];
  CHECK(PSB_NAME(vect_set_vect)(x, local, nc));
}

static void check_vector(vector_t *x, psb_c_descriptor *descriptor,
                         const scalar_t *expected, const char *what)
{
  scalar_t local[MAX_ENTRIES] = {0};
  psb_c_l_t owned[MAX_ENTRIES] = {0};
  const psb_c_i_t nr = psb_c_cd_get_local_rows(descriptor);
  const psb_c_i_t n = PSB_NAME(vect_get_nrows)(x);
  int valid = 1;
  require(n >= nr && n <= MAX_ENTRIES, "wrong vector storage size");
  CHECK(psb_c_cd_get_global_indices(owned, nr, true, descriptor));
  CHECK(PSB_NAME(vect_f_get_cpy)(local, x));
  for (psb_c_i_t i = 0; i < nr; ++i) {
    const scalar_t wanted = expected[owned[i] - index_base];
    const double difference = SCALAR_ABS(local[i] - wanted);
    if (!(difference <= TEST_TOL * (1. + SCALAR_ABS(wanted)))) valid = 0;
  }
  require(valid, what);
}

static void check_matvec(matrix_t *a, psb_c_descriptor *rows,
                         psb_c_descriptor *columns, const scalar_t *dense,
                         int nr, int nc)
{
  const char *transposes[] = {"N", "T", "C"};
  const scalar_t alpha = SCALAR_VALUE(1.25, .5);
  const scalar_t beta = SCALAR_VALUE(-.5, .25);
  for (int operation = 0; operation < 3; ++operation) {
    vector_t x = {NULL}, y = {NULL};
    psb_c_descriptor *xd = operation == 0 ? columns : rows;
    psb_c_descriptor *yd = operation == 0 ? rows : columns;
    const int nx = operation == 0 ? nc : nr;
    const int ny = operation == 0 ? nr : nc;
    scalar_t global_x[VELOCITY_SIZE] = {0};
    scalar_t expected_y[VELOCITY_SIZE] = {0};
    scalar_t next_y[VELOCITY_SIZE] = {0};
    const scalar_t zeros[VELOCITY_SIZE] = {0};
    CHECK(PSB_NAME(geall)(&x, xd));
    CHECK(PSB_NAME(geall)(&y, yd));
    CHECK(PSB_NAME(geasb)(&x, xd));
    CHECK(PSB_NAME(geasb)(&y, yd));
    for (int i = 0; i < nx; ++i) global_x[i] = SCALAR_VALUE(i + 1., .2 * i);
    for (int i = 0; i < ny; ++i) expected_y[i] = SCALAR_VALUE(2. - .5 * i, -.1 * i);
    set_vector(&y, yd, expected_y, SCALAR_VALUE(73., -19.));
    for (int repeat = 0; repeat < 3; ++repeat) {
      /* Nonzero halo garbage checks that every application exchanges/reset
       * the relevant halo, including repeated transposed operations.
       */
      set_vector(&x, xd, global_x, SCALAR_VALUE(101. + repeat, 31.));
      for (int i = 0; i < ny; ++i) {
        scalar_t sum = 0;
        for (int j = 0; j < nx; ++j) {
          scalar_t entry = operation == 0 ? dense[i * nc + j]
                                          : dense[j * nc + i];
          if (operation == 2) entry = SCALAR_CONJ(entry);
          sum += entry * global_x[j];
        }
        next_y[i] = alpha * sum + beta * expected_y[i];
      }
      CHECK(PSB_NAME(spmm_2desc)(alpha, a, &x, beta, &y,
                                rows, columns, transposes[operation]));
      check_vector(&y, yd, next_y, "rectangular matrix-vector product");
      check_vector(&x, xd, global_x, "matrix-vector product changed owned input");
      memcpy(expected_y, next_y, (size_t)ny * sizeof(expected_y[0]));
    }
    CHECK(PSB_NAME(vect_zero)(&x));
    check_vector(&x, xd, zeros, "vector zero retained nonzero entries");
    set_vector(&x, xd, zeros, SCALAR_VALUE(-67., 23.));
    CHECK(PSB_NAME(spmm_2desc)((scalar_t)1, a, &x, (scalar_t)0, &y,
                              rows, columns, transposes[operation]));
    check_vector(&y, yd, zeros, "zero input retained stale matrix-vector result");
    expect_error(PSB_NAME(spmm_2desc)(alpha, a, &x, beta, &x,
                                     rows, columns, transposes[operation]),
                 "matrix-vector input/output alias accepted");
    expect_error(PSB_NAME(spmm_2desc)(alpha, a, &x, beta, &y,
                                     rows, columns, "X"),
                 "invalid matrix-vector transpose accepted");
    expect_error(PSB_NAME(spmm_2desc)(alpha, a, &x, beta, &y,
                                     rows, columns,
                                     rank_id == process_count - 1 ? "X" : transposes[operation]),
                 "matrix-vector error was not agreed across ranks");
    CHECK(PSB_NAME(gefree)(&y, yd));
    CHECK(PSB_NAME(gefree)(&x, xd));
  }
}

static void run_case(int base)
{
  /* G is deliberately different from transpose(B).  In the complex builds,
   * nonzero imaginary parts distinguish transpose from conjugate transpose.
   */
  const scalar_t b[PRESSURE_SIZE * VELOCITY_SIZE] = {
    SCALAR_VALUE(1., .25), 0, 2, 0, -1,
    0, SCALAR_VALUE(3., -.5), 0, 4, 0,
    5, 0, 0, -2, SCALAR_VALUE(1., .75)
  };
  const scalar_t g[VELOCITY_SIZE * PRESSURE_SIZE] = {
    2, 0, SCALAR_VALUE(5., -.5),
    0, 1, 0,
    SCALAR_VALUE(1., .25), 0, 1,
    0, 2, -1,
    -1, 0, SCALAR_VALUE(3., .5)
  };
  scalar_t bt[VELOCITY_SIZE * PRESSURE_SIZE];
  scalar_t bh[VELOCITY_SIZE * PRESSURE_SIZE];
  scalar_t gt[PRESSURE_SIZE * VELOCITY_SIZE];
  scalar_t dbt[VELOCITY_SIZE * PRESSURE_SIZE];
  scalar_t dg[VELOCITY_SIZE * PRESSURE_SIZE];
  scalar_t lb[PRESSURE_SIZE * PRESSURE_SIZE];
  scalar_t lg[PRESSURE_SIZE * PRESSURE_SIZE];
  scalar_t lbg[PRESSURE_SIZE * PRESSURE_SIZE];
  matrix_t B = {NULL}, G = {NULL};
  matrix_t BT = {NULL}, BH = {NULL}, GT = {NULL};
  matrix_t DBT = {NULL}, DG = {NULL};
  matrix_t LB = {NULL}, LG = {NULL}, LBG = {NULL}, extra = {NULL};
  psb_c_descriptor br = {NULL}, bc = {NULL}, gr = {NULL}, gc = {NULL};
  psb_c_descriptor btc = {NULL}, bhc = {NULL}, gtc = {NULL};
  psb_c_descriptor lbc = {NULL}, lgc = {NULL}, lbgc = {NULL}, extra_c = {NULL};
  struct descriptor_snapshot before[4];
  scalar_t dummy = 1;

  index_base = base;
  psb_c_set_index_base(base);
  transpose_dense(b, PRESSURE_SIZE, VELOCITY_SIZE, bt, 0);
  transpose_dense(b, PRESSURE_SIZE, VELOCITY_SIZE, bh, 1);
  transpose_dense(g, VELOCITY_SIZE, PRESSURE_SIZE, gt, 0);
  for (int i = 0; i < VELOCITY_SIZE; ++i)
    for (int j = 0; j < PRESSURE_SIZE; ++j) {
      dbt[i * PRESSURE_SIZE + j] = (scalar_t)(i + 1) * bt[i * PRESSURE_SIZE + j];
      dg[i * PRESSURE_SIZE + j] = (scalar_t)(i + 1) * g[i * PRESSURE_SIZE + j];
    }
  multiply_dense(b, dbt, PRESSURE_SIZE, VELOCITY_SIZE, PRESSURE_SIZE, lb);
  multiply_dense(gt, dg, PRESSURE_SIZE, VELOCITY_SIZE, PRESSURE_SIZE, lg);
  multiply_dense(b, dg, PRESSURE_SIZE, VELOCITY_SIZE, PRESSURE_SIZE, lbg);

  assemble(b, PRESSURE_SIZE, VELOCITY_SIZE, &B, &br, &bc, 0);
  assemble(g, VELOCITY_SIZE, PRESSURE_SIZE, &G, &gr, &gc, 1);
  before[0] = snapshot(&br);
  before[1] = snapshot(&bc);
  before[2] = snapshot(&gr);
  before[3] = snapshot(&gc);
  check_matrix(&B, &br, &bc, PRESSURE_SIZE, VELOCITY_SIZE, b, "B assembly");
  check_matrix(&G, &gr, &gc, VELOCITY_SIZE, PRESSURE_SIZE, g, "G assembly");
  check_export_capacity(&B, &br, &bc);
  check_matvec(&B, &br, &bc, b, PRESSURE_SIZE, VELOCITY_SIZE);
  check_matvec(&G, &gr, &gc, g, VELOCITY_SIZE, PRESSURE_SIZE);

  CHECK(PSB_NAME(glob_transpose)(&B, &br, &bc, &BT, &btc, "T"));
  CHECK(PSB_NAME(glob_transpose)(&B, &br, &bc, &BH, &bhc, "C"));
  CHECK(PSB_NAME(glob_transpose)(&G, &gr, &gc, &GT, &gtc, "T"));
  check_matrix(&BT, &bc, &btc, VELOCITY_SIZE, PRESSURE_SIZE, bt, "B transpose");
  check_matrix(&BH, &bc, &bhc, VELOCITY_SIZE, PRESSURE_SIZE, bh,
               "B conjugate transpose");
  check_matrix(&GT, &gc, &gtc, PRESSURE_SIZE, VELOCITY_SIZE, gt, "G transpose");

  copy_and_scale(&BT, &bc, &DBT);
  copy_and_scale(&G, &gr, &DG);
  check_matrix(&DBT, &bc, &btc, VELOCITY_SIZE, PRESSURE_SIZE, dbt, "D B transpose");
  check_matrix(&DG, &gr, &gc, VELOCITY_SIZE, PRESSURE_SIZE, dg, "D G");

  CHECK(PSB_NAME(par_spspmm)(&B, &br, &bc, &DBT, &bc, &btc, &LB, &lbc));
  CHECK(PSB_NAME(par_spspmm)(&GT, &gc, &gtc, &DG, &gr, &gc, &LG, &lgc));
  CHECK(PSB_NAME(par_spspmm)(&B, &br, &bc, &DG, &gr, &gc, &LBG, &lbgc));
  check_matrix(&LB, &br, &lbc, PRESSURE_SIZE, PRESSURE_SIZE, lb, "B D B transpose");
  check_matrix(&LG, &gc, &lgc, PRESSURE_SIZE, PRESSURE_SIZE, lg, "G transpose D G");
  check_matrix(&LBG, &br, &lbgc, PRESSURE_SIZE, PRESSURE_SIZE, lbg, "B D G");
  check_matrix(&B, &br, &bc, PRESSURE_SIZE, VELOCITY_SIZE, b, "B was modified");
  check_matrix(&G, &gr, &gc, VELOCITY_SIZE, PRESSURE_SIZE, g, "G was modified");
  check_matrix(&BT, &bc, &btc, VELOCITY_SIZE, PRESSURE_SIZE, bt, "B transpose was modified");
  check_matrix(&DG, &gr, &gc, VELOCITY_SIZE, PRESSURE_SIZE, dg, "product changed right operand");
  check_descriptor(&br, &before[0]);
  check_descriptor(&bc, &before[1]);
  check_descriptor(&gr, &before[2]);
  check_descriptor(&gc, &before[3]);

  /* Repeated construction returns independently owned output objects. */
  for (int repeat = 0; repeat < 3; ++repeat) {
    CHECK(PSB_NAME(par_spspmm)(&B, &br, &bc, &DG, &gr, &gc, &extra, &extra_c));
    check_matrix(&extra, &br, &extra_c, PRESSURE_SIZE, PRESSURE_SIZE, lbg,
                 "repeated product");
    CHECK(PSB_NAME(spfree)(&extra, &br));
    CHECK(psb_c_cdfree(&extra_c));
  }

  /* Validation errors, including errors present on only one rank, must be
   * agreed before entering the native distributed numerical operation.
   */
  expect_error(PSB_NAME(glob_transpose)(&B, &br, &bc, &extra, &extra_c, "X"),
               "invalid transpose flag accepted");
  expect_error(PSB_NAME(glob_transpose)(&B, &br, &bc, &extra, &extra_c,
                                       rank_id == process_count - 1 ? "X" : "T"),
               "transpose error was not agreed across ranks");
  expect_error(PSB_NAME(glob_transpose)(&B, &br, &bc, &B, &extra_c, "T"),
               "matrix output alias accepted");
  expect_error(PSB_NAME(glob_transpose)(&B, &br, &bc, &extra, &br, "T"),
               "descriptor output alias accepted");
  expect_error(PSB_NAME(par_spspmm)(&B, &br, &bc, &G, &gr, &gc, &G, &extra_c),
               "product output alias accepted");
  expect_error(PSB_NAME(par_spspmm)(&B, &br, &bc, &G, &gr, &gc,
                                   rank_id == process_count - 1 ? &G : &extra, &extra_c),
               "product alias error was not agreed across ranks");
  expect_error(PSB_NAME(par_spspmm)(&B, &br, &bc, &B, &br, &bc, &extra, &extra_c),
               "incompatible inner dimension accepted");
  expect_error(PSB_NAME(sprowscl)(&B, -1, &dummy),
               "negative row scaling size accepted");
  check_matrix(&B, &br, &bc, PRESSURE_SIZE, VELOCITY_SIZE, b,
               "invalid call changed input matrix");

  CHECK(PSB_NAME(spfree)(&LBG, &br));
  CHECK(PSB_NAME(spfree)(&LG, &gc));
  CHECK(PSB_NAME(spfree)(&LB, &br));
  CHECK(PSB_NAME(spfree)(&DG, &gr));
  CHECK(PSB_NAME(spfree)(&DBT, &bc));
  CHECK(PSB_NAME(spfree)(&GT, &gc));
  CHECK(PSB_NAME(spfree)(&BH, &bc));
  CHECK(PSB_NAME(spfree)(&BT, &bc));
  CHECK(PSB_NAME(spfree)(&G, &gr));
  CHECK(PSB_NAME(spfree)(&B, &br));
  CHECK(psb_c_cdfree(&lbgc));
  CHECK(psb_c_cdfree(&lgc));
  CHECK(psb_c_cdfree(&lbc));
  CHECK(psb_c_cdfree(&gtc));
  CHECK(psb_c_cdfree(&bhc));
  CHECK(psb_c_cdfree(&btc));
  CHECK(psb_c_cdfree(&gc));
  CHECK(psb_c_cdfree(&gr));
  CHECK(psb_c_cdfree(&bc));
  CHECK(psb_c_cdfree(&br));
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
  MPI_Comm_size(MPI_COMM_WORLD, &process_count);
  psb_c_init(&context);
  psb_c_seterraction_ret();
  run_case(0);
  run_case(1);
  if (rank_id == 0)
    printf("PASS sparse_ops[%s]: %d rank(s), index bases 0 and 1\n",
           SCALAR_NAME, process_count);
  psb_c_exit_ctxt(context);
  MPI_Finalize();
  return EXIT_SUCCESS;
}
