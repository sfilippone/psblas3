#if (defined(_WIN32) || defined(WIN32))
#include <Windows.h>
#else
#include <sys/time.h>
#endif
#include <stdio.h>
#include <string.h>
#include "psb_internals.h"


#ifdef LowerUnderscore
#define mpi_wtime       mpi_wtime_
#define mpi_send        mpi_send_
#define mpi_isend       mpi_isend_
#define mpi_irecv       mpi_irecv_
#define mpi_wait        mpi_wait_
#define mpi_alltoall    mpi_alltoall_
#define mpi_alltoallv   mpi_alltoallv_
#define mpi_gather      mpi_gather_
#define mpi_gatherv     mpi_gatherv_
#define mpi_allgather   mpi_allgather_
#define mpi_allgatherv  mpi_allgatherv_
#define mpi_scatterv    mpi_scatterv_
#define mpi_scatter     mpi_scatter_
#endif
#ifdef LowerDoubleUnderscore
#define mpi_wtime       mpi_wtime__
#define mpi_send        mpi_send__
#define mpi_isend       mpi_isend__
#define mpi_irecv       mpi_irecv__
#define mpi_wait        mpi_wait__
#define mpi_alltoall    mpi_alltoall__
#define mpi_alltoallv   mpi_alltoallv__
#define mpi_gather      mpi_gather__
#define mpi_gatherv     mpi_gatherv__
#define mpi_allgather   mpi_allgather__
#define mpi_allgatherv  mpi_allgatherv__
#define mpi_scatterv    mpi_scatterv__
#define mpi_scatter     mpi_scatter__
#endif
#ifdef LowerCase
#define mpi_wtime       mpi_wtime
#define mpi_send        mpi_send
#define mpi_isend       mpi_isend
#define mpi_irecv       mpi_irecv
#define mpi_wait        mpi_wait
#define mpi_alltoall    mpi_alltoall
#define mpi_alltoallv   mpi_alltoallv
#define mpi_gather      mpi_gather
#define mpi_gatherv     mpi_gatherv
#define mpi_allgather   mpi_allgather
#define mpi_allgatherv  mpi_allgatherv
#define mpi_scatterv    mpi_scatterv
#define mpi_scatter     mpi_scatter
#endif
#ifdef UpperUnderscore
#define mpi_wtime       MPI_WTIME_
#define mpi_send        MPI_SEND_
#define mpi_isend       MPI_ISEND_
#define mpi_irecv       MPI_IRECV_
#define mpi_wait        MPI_WAIT_
#define mpi_alltoall    MPI_ALLTOALL_
#define mpi_alltoallv   MPI_ALLTOALLV_
#define mpi_gather      MPI_GATHER_
#define mpi_gatherv     MPI_GATHERV_
#define mpi_allgather   MPI_ALLGATHER_
#define mpi_allgatherv  MPI_ALLGATHERV_
#define mpi_scatterv    MPI_SCATTERV_
#define mpi_scatter     MPI_SCATTER_
#endif
#ifdef UpperDoubleUnderscore 
#define mpi_wtime       MPI_WTIME__
#define mpi_send        MPI_SEND__
#define mpi_isend       MPI_ISEND__
#define mpi_irecv       MPI_IRECV__
#define mpi_wait        MPI_WAIT__
#define mpi_alltoall    MPI_ALLTOALL__
#define mpi_alltoallv   MPI_ALLTOALLV__
#define mpi_gather      MPI_GATHER__
#define mpi_gatherv     MPI_GATHERV__
#define mpi_allgather   MPI_ALLGATHER__
#define mpi_allgatherv  MPI_ALLGATHERV__
#define mpi_scatterv    MPI_SCATTERV__
#define mpi_scatter     MPI_SCATTER__
#endif
#ifdef UpperCase
#define mpi_wtime       MPI_WTIME
#define mpi_send        MPI_SEND
#define mpi_isend       MPI_ISEND
#define mpi_irecv       MPI_IRECV
#define mpi_wait        MPI_WAIT
#define mpi_alltoall    MPI_ALLTOALL
#define mpi_alltoallv   MPI_ALLTOALLV
#define mpi_gather      MPI_GATHER
#define mpi_gatherv     MPI_GATHERV
#define mpi_allgather   MPI_ALLGATHER
#define mpi_allgatherv  MPI_ALLGATHERV
#define mpi_scatterv    MPI_SCATTERV
#define mpi_scatter     MPI_SCATTER
#endif

#define mpi_integer        1
#define mpi_integer8       2
#define mpi_real           3
#define mpi_double         4
#define mpi_complex        5
#define mpi_double_complex 6

double mpi_wtime()
{
#if defined(WIN32) || defined(_WIN32)
  LARGE_INTEGER tim, freq;
  double seconds;

  QueryPerformanceCounter(&tim);
  QeryPerformanceFrequency(&freq);
  seconds = (double)tim / (double)freq;
  return(seconds);
#else
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
#endif
}


void mpi_wait()
{
  return;
}
void mpi_send()
{
  return;
}
void mpi_isend()
{
  return;
}
void mpi_irecv()
{
  return;
}


void mpi_alltoall(void* sdb, int* sdc, int* sdt,
		  void* rvb, int* rvc, int* rvt, int* comm, int* ierr)
{
  int i,j,k; 
  
  if (*sdt == mpi_integer) {
    memcpy(rvb,sdb, (*sdc)*sizeof(int));
  }
  if (*sdt == mpi_integer8) {
    memcpy(rvb,sdb, (*sdc)*sizeof(long long));
  }
  if (*sdt == mpi_real) {
    memcpy(rvb,sdb, (*sdc)*sizeof(float));
  } 
  if (*sdt == mpi_double) {
    memcpy(rvb,sdb, (*sdc)*sizeof(double));
  } 
  if (*sdt == mpi_complex) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(float));
  }
  if (*sdt == mpi_double_complex) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(double));
  }
  *ierr = 0;
}

void mpi_alltoallv(void* sdb, int* sdc, int* sdspl, int* sdt,
		   void* rvb, int* rvc, int* rdspl, int* rvt, int* comm, int* ierr)
{
  int i,j,k; 

  
  if (*sdt == mpi_integer) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(int)),
	   (void *)((char *)sdb+sdspl[0]*sizeof(int)),(*sdc)*sizeof(int));
  }
  if (*sdt == mpi_integer8) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(long long)),
	   (void *)((char *)sdb+sdspl[0]*sizeof(long long)),(*sdc)*sizeof(long long));
  }
  if (*sdt == mpi_real) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(float)),
	   (void *)((char *)sdb+sdspl[0]*sizeof(float)),(*sdc)*sizeof(float));
  } 
  if (*sdt == mpi_double) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(double)),
	   (void *)((char *)sdb+sdspl[0]*sizeof(double)),(*sdc)*sizeof(double));
  } 
  if (*sdt == mpi_complex) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(float)),
	   (void *)((char *)sdb+sdspl[0]*2*sizeof(float)),(*sdc)*2*sizeof(float));
  }
  if (*sdt == mpi_double_complex) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(double)),
	   (void *)((char *)sdb+sdspl[0]*2*sizeof(double)),(*sdc)*2*sizeof(double));
  }
  *ierr = 0;
}


void mpi_gather(void* sdb, int* sdc, int* sdt,
		void* rvb, int* rvc, int* rvt, int *root, int* comm, int* ierr)
{
  int i,j,k; 
  
  if (*sdt == mpi_integer) {
    memcpy(rvb,sdb, (*sdc)*sizeof(int));
  }
  if (*sdt == mpi_integer8) {
    memcpy(rvb,sdb, (*sdc)*sizeof(long long));
  }
  if (*sdt == mpi_real) {
    memcpy(rvb,sdb, (*sdc)*sizeof(float));
  } 
  if (*sdt == mpi_double) {
    memcpy(rvb,sdb, (*sdc)*sizeof(double));
  } 
  if (*sdt == mpi_complex) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(float));
  }
  if (*sdt == mpi_double_complex) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(double));
  }
  *ierr = 0;
}


void mpi_gatherv(void* sdb, int* sdc, int* sdt,
		 void* rvb, int* rvc, int* rdspl,
		 int* rvt, int* comm, int *root, int* ierr)
{
  int i,j,k; 
  
  if (*sdt == mpi_integer) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(int)),
	   (void *)((char *)sdb),(*sdc)*sizeof(int));
  }
  if (*sdt == mpi_integer8) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(long long)),
	   (void *)((char *)sdb),(*sdc)*sizeof(long long));
  }
  if (*sdt == mpi_real) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(float)),
	   (void *)((char *)sdb),(*sdc)*sizeof(float));
  } 
  if (*sdt == mpi_double) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(double)),
	   (void *)((char *)sdb),(*sdc)*sizeof(double));
  } 
  if (*sdt == mpi_complex) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(float)),
	   (void *)((char *)sdb),(*sdc)*2*sizeof(float));
  }
  if (*sdt == mpi_double_complex) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(double)),
	   (void *)((char *)sdb),(*sdc)*2*sizeof(double));
  }


  *ierr = 0;
}


void mpi_scatter(void* sdb, int* sdc, int* sdt,
		 void* rvb, int* rvc, int* rvt, int *root, int* comm, int* ierr)
{
  int i,j,k; 
  
  if (*sdt == mpi_integer) {
    memcpy(rvb,sdb, (*sdc)*sizeof(int));
  }
  if (*sdt == mpi_integer8) {
    memcpy(rvb,sdb, (*sdc)*sizeof(long long));
  }
  if (*sdt == mpi_real) {
    memcpy(rvb,sdb, (*sdc)*sizeof(float));
  } 
  if (*sdt == mpi_double) {
    memcpy(rvb,sdb, (*sdc)*sizeof(double));
  } 
  if (*sdt == mpi_complex) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(float));
  }
  if (*sdt == mpi_double_complex) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(double));
  }
  *ierr = 0;
}


void mpi_scatterv(void* sdb, int* sdc, int* sdt,
		  void* rvb, int* rvc, int* rdspl,
		  int* rvt, int* comm, int *root, int* ierr)
{
  int i,j,k; 
  
  if (*sdt == mpi_integer) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(int)),
	   (void *)((char *)sdb),(*sdc)*sizeof(int));
  }
  if (*sdt == mpi_integer8) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(long long)),
	   (void *)((char *)sdb),(*sdc)*sizeof(long long));
  }
  if (*sdt == mpi_real) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(float)),
	   (void *)((char *)sdb),(*sdc)*sizeof(float));
  } 
  if (*sdt == mpi_double) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(double)),
	   (void *)((char *)sdb),(*sdc)*sizeof(double));
  } 
  if (*sdt == mpi_complex) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(float)),
	   (void *)((char *)sdb),(*sdc)*2*sizeof(float));
  }
  if (*sdt == mpi_double_complex) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(double)),
	   (void *)((char *)sdb),(*sdc)*2*sizeof(double));
  }


  *ierr = 0;
}


void mpi_allgather(void* sdb, int* sdc, int* sdt,
		   void* rvb, int* rvc, int* rvt, int* comm, int* ierr)
{
  int i,j,k; 
  
  if (*sdt == mpi_integer) {
    memcpy(rvb,sdb, (*sdc)*sizeof(int));
  }
  if (*sdt == mpi_integer8) {
    memcpy(rvb,sdb, (*sdc)*sizeof(long long));
  }
  if (*sdt == mpi_real) {
    memcpy(rvb,sdb, (*sdc)*sizeof(float));
  } 
  if (*sdt == mpi_double) {
    memcpy(rvb,sdb, (*sdc)*sizeof(double));
  } 
  if (*sdt == mpi_complex) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(float));
  }
  if (*sdt == mpi_double_complex) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(double));
  }
  *ierr = 0;
}

void mpi_allgatherv(void* sdb, int* sdc, int* sdt,
		    void* rvb, int* rvc, int* rdspl,
		    int* rvt, int* comm, int* ierr)
{
  int i,j,k; 
  
  if (*sdt == mpi_integer) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(int)),
	   (void *)((char *)sdb),(*sdc)*sizeof(int));
  }
  if (*sdt == mpi_integer8) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(long long)),
	   (void *)((char *)sdb),(*sdc)*sizeof(long long));
  }
  if (*sdt == mpi_real) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(float)),
	   (void *)((char *)sdb),(*sdc)*sizeof(float));
  } 
  if (*sdt == mpi_double) {
    memcpy((void *)((char *)rvb+rdspl[0]*sizeof(double)),
	   (void *)((char *)sdb),(*sdc)*sizeof(double));
  } 
  if (*sdt == mpi_complex) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(float)),
	   (void *)((char *)sdb),(*sdc)*2*sizeof(float));
  }
  if (*sdt == mpi_double_complex) {
    memcpy((void *)((char *)rvb+rdspl[0]*2*sizeof(double)),
	   (void *)((char *)sdb),(*sdc)*2*sizeof(double));
  }


  *ierr = 0;
}
