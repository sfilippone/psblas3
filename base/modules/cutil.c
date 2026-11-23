#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

void psi_c_diffadd(void *p1, void *p2, int64_t *ret)
{
  *ret = (int)((char *)p2-(char *)p1);
  return;
}

/*
 * Point-to-point MPI on raw buffers given as C pointers, used by the halo
 * exchange when the communication buffer lives in device memory and MPI is
 * CUDA-aware (the Fortran MPI bindings cannot take a device address).
 * FIRST is the 1-based index of the first element of the message in BUF;
 * DTYPE, COMM and REQ are Fortran MPI handles.
 */
#include "psb_config.h"
#include <stdlib.h>
#include <strings.h>
#if !defined(PSB_SERIAL_MPI)
#include <mpi.h>
#if defined(__has_include)
#if __has_include(<mpi-ext.h>)
#include <mpi-ext.h>
#endif
#endif

static char *psi_c_elem_ptr(void *buf, int first, MPI_Datatype dt)
{
  int sz = 0;
  MPI_Type_size(dt, &sz);
  return ((char *) buf) + ((size_t) (first - 1)) * ((size_t) sz);
}

int psi_c_dev_isend(void *buf, int first, int count, MPI_Fint dtype, int dest,
                    int tag, MPI_Fint comm, MPI_Fint *req)
{
  MPI_Datatype dt = MPI_Type_f2c(dtype);
  MPI_Request  r;
  int ret = MPI_Isend(psi_c_elem_ptr(buf, first, dt), count, dt, dest, tag,
                      MPI_Comm_f2c(comm), &r);
  *req = MPI_Request_c2f(r);
  return ret;
}

int psi_c_dev_irecv(void *buf, int first, int count, MPI_Fint dtype, int src,
                    int tag, MPI_Fint comm, MPI_Fint *req)
{
  MPI_Datatype dt = MPI_Type_f2c(dtype);
  MPI_Request  r;
  int ret = MPI_Irecv(psi_c_elem_ptr(buf, first, dt), count, dt, src, tag,
                      MPI_Comm_f2c(comm), &r);
  *req = MPI_Request_c2f(r);
  return ret;
}

/*
 * Does the MPI library report CUDA support at run time?
 *   1 yes, 0 no, -1 unknown (no query available for this MPI)
 */
int psi_c_mpi_cuda_query(void)
{
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  /* Open MPI: compile-time macro + run-time query */
  return (MPIX_Query_cuda_support() ? 1 : 0);
#elif defined(MPIX_GPU_SUPPORT_CUDA)
  /* MPICH >= 4 */
  int flag = 0;
  if (MPIX_GPU_query_support(MPIX_GPU_SUPPORT_CUDA, &flag) != MPI_SUCCESS) return -1;
  return (flag ? 1 : 0);
#elif defined(OPEN_MPI) && OPEN_MPI
  /* Open MPI built without CUDA support */
  return 0;
#else
  return -1;
#endif
}

/*
 * Requested mode for the halo exchange on device buffers, from
 * PSB_CUDA_AWARE and the MPI query (the self test is done afterwards):
 *   1  try it
 *   0  off (PSB_CUDA_AWARE=0, or auto without reported support)
 *   2  off: PSB_CUDA_AWARE=1 but MPI reports no CUDA support
 *   3  try it: PSB_CUDA_AWARE=1 and this MPI cannot be queried
 */
int psi_c_cuda_aware_mpi(void)
{
  const char *e = getenv("PSB_CUDA_AWARE");
  int q = psi_c_mpi_cuda_query();
  if (e != NULL) {
    if (!strcasecmp(e, "0") || !strcasecmp(e, "off") || !strcasecmp(e, "no"))
      return 0;
    if (!strcasecmp(e, "1") || !strcasecmp(e, "on") || !strcasecmp(e, "yes")) {
      if (q == 0)  return 2;
      if (q == -1) return 3;
      return 1;
    }
  }
  return (q == 1 ? 1 : 0);
}

/* Rank of the calling process among the processes of its node */
int psi_c_node_local_rank(void)
{
  MPI_Comm node;
  int rank = 0, init = 0;
  MPI_Initialized(&init);
  if (!init) return 0;
  if (MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                          MPI_INFO_NULL, &node) != MPI_SUCCESS) return 0;
  MPI_Comm_rank(node, &rank);
  MPI_Comm_free(&node);
  return rank;
}
#else
int psi_c_dev_isend(void *buf, int first, int count, int dtype, int dest,
                    int tag, int comm, int *req) { return -1; }
int psi_c_dev_irecv(void *buf, int first, int count, int dtype, int src,
                    int tag, int comm, int *req) { return -1; }
int psi_c_mpi_cuda_query(void) { return 0; }
int psi_c_cuda_aware_mpi(void) { return 0; }
int psi_c_node_local_rank(void) { return 0; }
#endif
