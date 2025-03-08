#include <sys/time.h>
#include <stdio.h>
#include <string.h>
#include "psb_internals.h"

double mpi_wtime() 
{
  struct timeval tt;
  struct timezone tz;
  double temp;
  if (gettimeofday(&tt,&tz) != 0) {
    fprintf(stderr,"Fatal error for gettimeofday ??? \n");
    temp=0.0;
  } else {
    temp = ((double)tt.tv_sec) + ((double)tt.tv_usec)*1.0e-6;
  }
  return(temp);
}


void mpi_wait(int *request, int* status, int *ierr)

{
  *ierr = 0;  
  return;
}
void mpi_send(void* buf, int* count, int* datatype,
	      int *dest, int *tag, int *comm, int *ierr)
{
  *ierr = 0;
  return;
}
void mpi_isend(void* buf, int* count, int* datatype,
	       int *dest, int *tag, int *comm, int *request,
	       int *ierr)
{
  *ierr = 0;
  return;
}
void mpi_irecv(void* buf, int* count, int* datatype,
	       int *src, int *tag, int *comm, int *request,
	       int *ierr)
{
  *ierr = 0;
  return;
}


void mpi_alltoall(void* sdb, int* sdc, int* sdt,
		  void* rvb, int* rvc, int* rvt, int* comm, int* ierr)
{
  int i,j,k; 
  
  if ((*sdt == MPI_INTEGER)||(*sdt == MPI_INTEGER4)||(*sdt == MPI_LOGICAL)) {
    memcpy(rvb,sdb, (*sdc)*sizeof(int32_t));
  }
  if (*sdt == MPI_CHARACTER) {
    memcpy(rvb,sdb, (*sdc)*sizeof(char));
  }
  if (*sdt == MPI_INTEGER8) {
    memcpy(rvb,sdb, (*sdc)*sizeof(int64_t));
  }
  if (*sdt == MPI_INTEGER2) {
    memcpy(rvb,sdb, (*sdc)*sizeof(int16_t));
  }
  if (*sdt == MPI_REAL) {
    memcpy(rvb,sdb, (*sdc)*sizeof(float));
  } 
  if (*sdt == MPI_DOUBLE) {
    memcpy(rvb,sdb, (*sdc)*sizeof(double));
  } 
  if (*sdt == MPI_COMPLEX) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(float));
  }
  if (*sdt == MPI_DOUBLE_COMPLEX) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(double));
  }
  *ierr = 0;
}

void mpi_alltoallv(void* sdb, int* sdc, int* sdspl, int* sdt,
		   void* rvb, int* rvc, int* rdspl, int* rvt, int* comm, int* ierr)
{
  int i,j,k; 

  
  if ((*sdt == MPI_INTEGER)||(*sdt == MPI_INTEGER4)||(*sdt == MPI_LOGICAL)) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(int32_t)),
	   (void *)((char *)sdb+sdspl[0]*sizeof(int32_t)),(*sdc)*sizeof(int32_t));
  }
  if (*sdt == MPI_CHARACTER) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(char)),
	   (void *)((char *)sdb+sdspl[0]*sizeof(char)),(*sdc)*sizeof(char));
  }
  if (*sdt == MPI_INTEGER8) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(int64_t)),
	   (void *)((char *)sdb+sdspl[0]*sizeof(int64_t)),(*sdc)*sizeof(int64_t));
  }
  if (*sdt == MPI_REAL) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(float)),
	   (void *)((char *)sdb+sdspl[0]*sizeof(float)),(*sdc)*sizeof(float));
  } 
  if (*sdt == MPI_DOUBLE) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(double)),
	   (void *)((char *)sdb+sdspl[0]*sizeof(double)),(*sdc)*sizeof(double));
  } 
  if (*sdt == MPI_COMPLEX) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(float)),
	   (void *)((char *)sdb+sdspl[0]*2*sizeof(float)),(*sdc)*2*sizeof(float));
  }
  if (*sdt == MPI_DOUBLE_COMPLEX) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(double)),
	   (void *)((char *)sdb+sdspl[0]*2*sizeof(double)),(*sdc)*2*sizeof(double));
  }
  *ierr = 0;
}


void mpi_gather(void* sdb, int* sdc, int* sdt,
		void* rvb, int* rvc, int* rvt, int *root, int* comm, int* ierr)
{
  int i,j,k; 
  
  if ((*sdt == MPI_INTEGER)||(*sdt == MPI_INTEGER4)||(*sdt == MPI_LOGICAL)) {
    memcpy(rvb,sdb, (*sdc)*sizeof(int32_t));
  }
  if (*sdt == MPI_INTEGER8) {
    memcpy(rvb,sdb, (*sdc)*sizeof(int64_t));
  }
  if (*sdt == MPI_CHARACTER) {
    memcpy(rvb,sdb, (*sdc)*sizeof(char));
  }
  if (*sdt == MPI_REAL) {
    memcpy(rvb,sdb, (*sdc)*sizeof(float));
  } 
  if (*sdt == MPI_DOUBLE) {
    memcpy(rvb,sdb, (*sdc)*sizeof(double));
  } 
  if (*sdt == MPI_COMPLEX) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(float));
  }
  if (*sdt == MPI_DOUBLE_COMPLEX) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(double));
  }
  *ierr = 0;
}


void mpi_gatherv(void* sdb, int* sdc, int* sdt,
		 void* rvb, int* rvc, int* rdspl,
		 int* rvt, int* comm, int *root, int* ierr)
{
  int i,j,k; 
  
  if ((*sdt == MPI_INTEGER)||(*sdt == MPI_INTEGER4)||(*sdt == MPI_LOGICAL)) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(int32_t)),
	   (void *)((char *)sdb),(*sdc)*sizeof(int32_t));
  }
  if (*sdt == MPI_INTEGER8) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(int64_t)),
	   (void *)((char *)sdb),(*sdc)*sizeof(int64_t));
  }
  if (*sdt == MPI_CHARACTER) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(char)),
	   (void *)((char *)sdb),(*sdc)*sizeof(char));
  }
  if (*sdt == MPI_REAL) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(float)),
	   (void *)((char *)sdb),(*sdc)*sizeof(float));
  } 
  if (*sdt == MPI_DOUBLE) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(double)),
	   (void *)((char *)sdb),(*sdc)*sizeof(double));
  } 
  if (*sdt == MPI_COMPLEX) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(float)),
	   (void *)((char *)sdb),(*sdc)*2*sizeof(float));
  }
  if (*sdt == MPI_DOUBLE_COMPLEX) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(double)),
	   (void *)((char *)sdb),(*sdc)*2*sizeof(double));
  }


  *ierr = 0;
}


void mpi_scatter(void* sdb, int* sdc, int* sdt,
		 void* rvb, int* rvc, int* rvt, int *root, int* comm, int* ierr)
{
  int i,j,k; 
  
  if ((*sdt == MPI_INTEGER)||(*sdt == MPI_INTEGER4)||(*sdt == MPI_LOGICAL)) {
    memcpy(rvb,sdb, (*sdc)*sizeof(int32_t));
  }
  if (*sdt == MPI_CHARACTER) {
    memcpy(rvb,sdb, (*sdc)*sizeof(char));
  }
  if (*sdt == MPI_INTEGER8) {
    memcpy(rvb,sdb, (*sdc)*sizeof(int64_t));
  }
  if (*sdt == MPI_REAL) {
    memcpy(rvb,sdb, (*sdc)*sizeof(float));
  } 
  if (*sdt == MPI_DOUBLE) {
    memcpy(rvb,sdb, (*sdc)*sizeof(double));
  } 
  if (*sdt == MPI_COMPLEX) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(float));
  }
  if (*sdt == MPI_DOUBLE_COMPLEX) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(double));
  }
  *ierr = 0;
}


void mpi_scatterv(void* sdb, int* sdc, int* sdspl, int* sdt,
		  void* rvb, int* rvc, 
		  int* rvt, int* comm, int *root, int* ierr)
{
  int i,j,k; 
  
  if ((*sdt == MPI_INTEGER)||(*sdt == MPI_INTEGER4)||(*sdt == MPI_LOGICAL)) {
    memcpy((void *)((char *)rvb+sdspl[0]*sizeof(int32_t)),
	   (void *)((char *)sdb),(*sdc)*sizeof(int32_t));
  }
  if (*sdt == MPI_CHARACTER) {
    memcpy((void *)((char *)rvb+sdspl[0]*sizeof(char)),
	   (void *)((char *)sdb),(*sdc)*sizeof(char));
  }
  if (*sdt == MPI_INTEGER8) {
    memcpy((void *)((char *)rvb+sdspl[0]*sizeof(int64_t)),
	   (void *)((char *)sdb),(*sdc)*sizeof(int64_t));
  }
  if (*sdt == MPI_REAL) {
    memcpy((void *)((char *)rvb+sdspl[0]*sizeof(float)),
	   (void *)((char *)sdb),(*sdc)*sizeof(float));
  } 
  if (*sdt == MPI_DOUBLE) {
    memcpy((void *)((char *)rvb+sdspl[0]*sizeof(double)),
	   (void *)((char *)sdb),(*sdc)*sizeof(double));
  } 
  if (*sdt == MPI_COMPLEX) {
    memcpy((void *)((char *)rvb+sdspl[0]*2*sizeof(float)),
	   (void *)((char *)sdb),(*sdc)*2*sizeof(float));
  }
  if (*sdt == MPI_DOUBLE_COMPLEX) {
    memcpy((void *)((char *)rvb+sdspl[0]*2*sizeof(double)),
	   (void *)((char *)sdb),(*sdc)*2*sizeof(double));
  }


  *ierr = 0;
}


void mpi_allgather(void* sdb, int* sdc, int* sdt,
		   void* rvb, int* rvc, int* rvt, int* comm, int* ierr)
{
  int i,j,k; 
  
  if ((*sdt == MPI_INTEGER)||(*sdt == MPI_INTEGER4)||(*sdt == MPI_LOGICAL)) {
    memcpy(rvb,sdb, (*sdc)*sizeof(int32_t));
  }
  if (*sdt == MPI_CHARACTER) {
    memcpy(rvb,sdb, (*sdc)*sizeof(char));
  }
  if (*sdt == MPI_INTEGER8) {
    memcpy(rvb,sdb, (*sdc)*sizeof(int64_t));
  }
  if (*sdt == MPI_REAL) {
    memcpy(rvb,sdb, (*sdc)*sizeof(float));
  } 
  if (*sdt == MPI_DOUBLE) {
    memcpy(rvb,sdb, (*sdc)*sizeof(double));
  } 
  if (*sdt == MPI_COMPLEX) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(float));
  }
  if (*sdt == MPI_DOUBLE_COMPLEX) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(double));
  }
  *ierr = 0;
}

void mpi_allgatherv(void* sdb, int* sdc, int* sdt,
		    void* rvb, int* rvc, int* rdspl,
		    int* rvt, int* comm, int* ierr)
{
  int i,j,k; 
  
  if ((*sdt == MPI_INTEGER)||(*sdt == MPI_INTEGER4)||(*sdt == MPI_LOGICAL)) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(int32_t)),
	   (void *)((char *)sdb),(*sdc)*sizeof(int32_t));
  }
  if (*sdt == MPI_CHARACTER) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(char)),
	   (void *)((char *)sdb),(*sdc)*sizeof(char));
  }
  if (*sdt == MPI_INTEGER8) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(int64_t)),
	   (void *)((char *)sdb),(*sdc)*sizeof(int64_t));
  }
  if (*sdt == MPI_REAL) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(float)),
	   (void *)((char *)sdb),(*sdc)*sizeof(float));
  } 
  if (*sdt == MPI_DOUBLE) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(double)),
	   (void *)((char *)sdb),(*sdc)*sizeof(double));
  } 
  if (*sdt == MPI_COMPLEX) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(float)),
	   (void *)((char *)sdb),(*sdc)*2*sizeof(float));
  }
  if (*sdt == MPI_DOUBLE_COMPLEX) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(double)),
	   (void *)((char *)sdb),(*sdc)*2*sizeof(double));
  }

  *ierr = 0;
}
