#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv)
{
	int n = 0; 
	int myrank=0, numproc = 0;
	MPI_Status estado;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   	MPI_Comm_size(MPI_COMM_WORLD, &numproc);
	if (myrank == 0){
		printf("Introduce un entero para transimitir: \n");
		scanf ("%u",&n);
		printf("El entero introducido es %d\n", n);
		MPI_Send(&n, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(&n, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &estado);
		printf("Soy el proceso 0. El entero que he recibido es: %d\n", n);
	} 
	else{
		MPI_Recv(&n, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &estado);
		printf("Soy el proceso %d. El entero que he recibido es: %d\n", myrank, n);
		++n;
		MPI_Send(&n, 1, MPI_INT, myrank+1<numproc ? myrank+1 : 0, 0, MPI_COMM_WORLD);
	}
	MPI_Finalize();
}