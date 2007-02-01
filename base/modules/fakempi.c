#include <sys/time.h>
#include <stdio.h>
#include <string.h>

#ifdef Add_
#define mpi_wtime       mpi_wtime_
#define mpi_send        mpi_send_
#define mpi_irecv       mpi_irecv_
#define mpi_wait        mpi_wait_
#define mpi_alltoall    mpi_alltoall_
#define mpi_alltoallv   mpi_alltoallv_
#define mpi_allgather   mpi_allgather_
#define mpi_allgatherv  mpi_allgatherv_
#endif
#define mpi_integer        1
#define mpi_double         3
#define mpi_double_complex 5

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


void mpi_wait()
{
  return;
}
void mpi_send()
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
  if (*sdt == mpi_double) {
    memcpy(rvb,sdb, (*sdc)*sizeof(double));
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
    memcpy((rvb+rdspl[0]*sizeof(int)),
	   (sdb+sdspl[0]*sizeof(int)),(*sdc)*sizeof(int));
  }
  if (*sdt == mpi_double) {
    memcpy((rvb+rdspl[0]*sizeof(double)),
	   (sdb+sdspl[0]*sizeof(double)),(*sdc)*sizeof(double));
  } 
  if (*sdt == mpi_double_complex) {
    memcpy((rvb+rdspl[0]*2*sizeof(double)),
	   (sdb+sdspl[0]*2*sizeof(double)),(*sdc)*2*sizeof(double));
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
  if (*sdt == mpi_double) {
    memcpy(rvb,sdb, (*sdc)*sizeof(double));
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
    memcpy(rvb,sdb, (*sdc)*sizeof(int));
  }
  if (*sdt == mpi_double) {
    memcpy(rvb,sdb, (*sdc)*sizeof(double));
  } 
  if (*sdt == mpi_double_complex) {
    memcpy(rvb,sdb, (*sdc)*2*sizeof(double));
  }
  if (*sdt == mpi_integer) {
    memcpy((rvb+rdspl[0]*sizeof(int)),
	   (sdb),(*sdc)*sizeof(int));
  }
  if (*sdt == mpi_double) {
    memcpy((rvb+rdspl[0]*sizeof(double)),
	   (sdb),(*sdc)*sizeof(double));
  } 
  if (*sdt == mpi_double_complex) {
    memcpy((rvb+rdspl[0]*2*sizeof(double)),
	   (sdb),(*sdc)*2*sizeof(double));
  }


  *ierr = 0;
}
