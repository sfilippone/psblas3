#ifndef PSB_FAKEMPI_H
#define PSB_FAKEMPI_H
#include "psb_config.h"

#define MPI_INTEGER        1
#define MPI_INTEGER8       2
#define MPI_REAL           3
#define MPI_DOUBLE         4
#define MPI_COMPLEX        5
#define MPI_DOUBLE_COMPLEX 6
#define MPI_CHARACTER      7
#define MPI_LOGICAL        8
#define MPI_INTEGER2       9
#define MPI_INTEGER4       10
#define MPI_COMM_NULL      -1
#define MPI_COMM_WORLD     1


double mpi_wtime();
void mpi_wait(int *request, int* status, int *ierr);
void mpi_send(void* buf, int* count, int* datatype,
	      int *dest, int *tag, int *comm, int *ierr);
void mpi_isend(void* buf, int* count, int* datatype,
	       int *dest, int *tag, int *comm, int *request,
	       int *ierr);
void mpi_irecv(void* buf, int* count, int* datatype,
	       int *src, int *tag, int *comm, int *request,
	       int *ierr);
void mpi_alltoall(void* sdb, int* sdc, int* sdt,
		  void* rvb, int* rvc, int* rvt, int* comm, int* ierr);
void mpi_alltoallv(void* sdb, int* sdc, int* sdspl, int* sdt,
		   void* rvb, int* rvc, int* rdspl, int* rvt, int* comm, int* ierr);
void mpi_gather(void* sdb, int* sdc, int* sdt,
		void* rvb, int* rvc, int* rvt, int *root, int* comm, int* ierr);
void mpi_gatherv(void* sdb, int* sdc, int* sdt,
		 void* rvb, int* rvc, int* rdspl,
		 int* rvt, int* comm, int *root, int* ierr);
void mpi_scatter(void* sdb, int* sdc, int* sdt,
		 void* rvb, int* rvc, int* rvt, int *root, int* comm, int* ierr);
void mpi_scatterv(void* sdb, int* sdc, int* sdspl, int* sdt,
		  void* rvb, int* rvc, 
		  int* rvt, int* comm, int *root, int* ierr);
void mpi_allgather(void* sdb, int* sdc, int* sdt,
		   void* rvb, int* rvc, int* rvt, int* comm, int* ierr);
void mpi_allgatherv(void* sdb, int* sdc, int* sdt,
		    void* rvb, int* rvc, int* rdspl,
		    int* rvt, int* comm, int* ierr);
#endif
