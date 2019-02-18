#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv)
{
	int n = 1000; 
	int myrank=0, numproc = 0;
    char c = 'a';
	double beta, start_time, end_time, times;
	MPI_Status estado;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   	MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    if (myrank == 0){
		for (int i =0;i<n;++i){
			start_time = MPI_Wtime() * 1000000;
			MPI_Send(&c, 1, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
			MPI_Recv(&c, 1, MPI_BYTE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &estado);
			end_time = MPI_Wtime() * 1000000;
			times += (end_time - start_time);
		}
		beta = times/n;
		printf("La latencia es de %f seg\n", beta);
    }
	else if (myrank==1){
		for (int i=0;i<n;++i){
			MPI_Recv(&c, 1, MPI_BYTE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &estado);
			MPI_Send(&c, 1, MPI_BYTE, 0, 1, MPI_COMM_WORLD);
		}
	}
	MPI_Finalize();
}