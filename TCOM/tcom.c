#include <stdio.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char **argv)
{
	// Inicialización variables 
	int min_size = 5;
	int max_size = 20; 
	int size;
	int myrank=0, numproc = 0;
	int n=1000;
	char c = 'a';
	double tcom[max_size];
	double tau[max_size];
	double beta, start_time, end_time, times;
	MPI_Status estado;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  if (myrank == 0){
		// Estimación de beta
		for (int i =0;i<n;++i){
			start_time = MPI_Wtime() * 1000000; // Inicio
			MPI_Send(&c, 1, MPI_BYTE, 1, 0, MPI_COMM_WORLD); // Envío 
			MPI_Recv(&c, 1, MPI_BYTE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &estado); // Respuesta
			end_time = MPI_Wtime() * 1000000; // Final 
			times += (end_time - start_time) / 2; // Tiempo acumulado
		}
		beta = times/n; // Media aritmética
		printf("La latencia es de %.4f micro s\n", beta);
		
		// Estimación de tau
		printf("Bytes\t\tTiempo(micro s)\t\tTau(micro s/byte)\n");
		for (int i=min_size;i<max_size;++i){ // Tamaños de los mensajes
			size = pow(2, i); // Tamaño del mensaje (en doubles)
			double array[size]; // Datos a enviar
			tau[i] = 0;
			tcom[i] = 0;
			for (int j=0;j<n;++j){ // 1000 repeticiones
				start_time = MPI_Wtime() * 1000000; // Inicio
				MPI_Send(&array, size, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD); // Envío
				MPI_Recv(&array, size, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &estado); // Respuesta
				end_time = MPI_Wtime() * 1000000; // Final 
				tcom[i] += (end_time - start_time) / 2; // Tiempo acumulado para el mensaje de tamaño pow(2, i)*sizeof(double)
			}
			tcom[i] = tcom[i]/n; // Tcom para el tamaño pow(2, i)*sizeof(double)
			tau[i] = (tcom[i] - beta) / (8*size); // Tau para el tamaño pow(2, i)*sizeof(double)
			printf("%d,\t\t%.4f,\t\t%.6f\n", size*8, tcom[i], tau[i]);
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