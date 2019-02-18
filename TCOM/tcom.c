#include <stdio.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char **argv)
{
	int n = 1000; 
	int min_size = 5;
	int max_size = 20; 
	int size;
	int myrank=0, numproc = 0;
    char c = 'a';
	double tcom[max_size];
	double tau[max_size];
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
			times += (end_time - start_time) / 2;
		}
		beta = times/n;
		printf("La latencia es de %.4f ms\n", beta);
		for (int i=min_size;i<max_size;++i){
			size = pow(2, i);
			double array[size];
			tau[i] = 0;
			tcom[i] = 0;
			for (int j=0;j<n;++j){
				start_time = MPI_Wtime() * 1000000;
				MPI_Send(&array, size, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
				MPI_Recv(&array, size, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &estado);
				end_time = MPI_Wtime() * 1000000;
				times += (end_time - start_time) / 2;
				tcom[i] += times;
			}
			tcom[i] = tcom[i]/n;
			tau[i] = (tcom[i] - beta) / (8*size);
			printf("T%d\t\tTcom=%.2f\tms\t\tTau=%.4f\tms/bytes\n", i, tcom[i], tau[i]);
		}

    }
	else if (myrank==1){
		for (int i=0;i<n;++i){
			MPI_Recv(&c, 1, MPI_BYTE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &estado);
			MPI_Send(&c, 1, MPI_BYTE, 0, 1, MPI_COMM_WORLD);
		}
		for (int i=min_size;i<max_size;++i){
			size = pow(2, i);
			double array[size];
			for (int j=0;j<n;++j){
				MPI_Recv(&array, size, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &estado);
				MPI_Send(&array, size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
			}
		}
	}
	MPI_Finalize();
}