#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define maxnprocs 8
#define maxn 100000000

int main(int argc, char** argv){
    int n = 0, myrank, np, k;
    MPI_Status estado;
    double *x, *y;
    double dotproduct = 0.0;
    double final = 0.0;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    if (np > maxnprocs){
        if (myrank == 0) printf("Se supera el número máximo de procesos\n");
        MPI_Finalize();
        return 0;
    }
    if (myrank == 0){
        // leer el tamaño del vector 
        printf("Indtroduce el tamaño del vector: \n");
        scanf("%u", &n);
        if (n>maxn){
            printf("Se supera el tamaño máximo del vector\n");
            k = -1;
        }   
        else k = n/np;
    }
    // enviar n a todos los procesos
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (k==-1){
        MPI_Finalize();
        return 0;
    }
    if (myrank == 0){
		// crear los vectores en el proceso 0 
		x = malloc(n * sizeof(double)); 
    	y = malloc(n * sizeof(double));
		for (int i=0;i<n;++i){
			x[i] = (double)1/(double)(i+1);
			y[i] = i+1;
		}
        // enviar partes correspondientes del vector a cada proceso 
		int indice = k+n%np; 
		for (int i=1;i<np;++i){
			MPI_Send(&x[indice], k, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
			MPI_Send(&y[indice], k, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
			indice += k;
		}
        k = k+n%np; // número de elementos para el proceso 0 
	} else {
		x = malloc(k*sizeof(double));
        y = malloc(k*sizeof(double));
        // recibir 
        MPI_Recv(&x[0], k, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD, &estado);
        MPI_Recv(&y[0], k, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD, &estado);
	}
    // calcular el producto escalar parcial en cada proceso 
    for (int i=0;i<k; ++i){
        dotproduct += x[i] * y[i];
    }
    printf("Soy el proceso %d y mi producto escalar es %.2f\n", myrank, dotproduct);
    // acumular la suma de todos los productos y enviarla al proceso 0
    MPI_Reduce(&dotproduct, &final, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myrank == 0){
        // printf("\nSoy el proceso 0 y el resultado final es %.2f\n\n", final);
        printf("\nEl producto escalar del vector es %.1f\n\n", final);
    }
    free(x);
    free(y);
    MPI_Finalize();
}