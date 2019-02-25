#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define maxnprocs 8
#define maxn 100000000

int main(int argc, char** argv){
    int n = 0, myrank, np, k;
    MPI_Status estado;
    double *x, *y;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (np > maxnprocs){
        printf("Se supera el número máximo de procesos\n");
        MPI_Finalize();
        return 0;
    }
    if (myrank == 0){
        // leer el tamaño del vector 
        printf("Indtroduce el tamaño del vector: \n");
        scanf("%u", &n);   
        if (n>maxn){
            printf("Se supera el tamaño máximo del vector\n");
            MPI_Finalize();
            return 0;
        }
    }
    // enviar n a todos los procesos
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    k = n/np; // número de elementos para cada proceso 
    printf("k: %d\n", k);
    if (myrank == 0){
		// crear los vectores en el proceso 0 
		x = malloc(n * sizeof(double));  // asignamos memoria dinámicamente con la funcion malloc perteneciente a stdlib.h
    	y = malloc(n * sizeof(double));
		for (int i=0;i<n;++i){
			x[i] = (double)1/(double)(i+1);
			y[i] = i+1;
            printf("x: \t %2.2f \t\t y: \t %2.2f\n", x[i], y[i]);
		}
		int indice = k+n%np; // número de elementos para el proceso 0 
		// enviar partes correspondientes del vector a cada proceso 
		for (int i=1;i<np;++i){
			MPI_Send(&x[indice], k, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
			MPI_Send(&y[indice], k, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
			indice += k;
		}
	} else {
		x = malloc(k*sizeof(double));
        y = malloc(k*sizeof(double));
        // recibir 
        MPI_Recv(&x[0], k, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD, &estado);
        MPI_Recv(&y[0], k, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD, &estado);

        printf("Soy el proceso %d, los elementos que he recibido son: \n", myrank);
        for (int i=0;i<k;i++){
            printf("myrank: %d \t x: \t %2.2f \t\t y: \t %2.2f\n", myrank, x[i], y[i]);
        }
	}
    MPI_Finalize();
}