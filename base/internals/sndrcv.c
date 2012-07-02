#include <stdio.h>
#include <mpi.h>

int receiveRoutine(double * y, int recvtype, int procSender,
		   int tag, int comm, int *handle){

	MPI_Comm co = MPI_Comm_f2c(comm);
	MPI_Datatype dt = MPI_Type_f2c(recvtype);
	MPI_Request req;// = MPI_Request_f2c(*handle);
	MPI_Irecv(y, 1, dt, procSender,tag, co, &req);
	*handle = MPI_Request_c2f(req);
	return 0;

}

int sendRoutine(double * y, int sendtype, int procToSend,int tag, int comm){

	MPI_Comm co = MPI_Comm_f2c(comm);
	MPI_Datatype dt = MPI_Type_f2c(sendtype);
	MPI_Send(y, 1, dt, procToSend,tag,co);
	return 0;
}

