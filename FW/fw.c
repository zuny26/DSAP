#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

#define maxnumprocs 8
#define maxn 1000


int main(int argc, char** argv){
	float** crearMatrizPesos(int, int);
	int** crearMatrizCaminos(int, int);
	void definirGrafo(int, float**, int**);
	void printMatrizPesos(float**, int, int);
	void printMatrizCaminos(int**, int, int);
	void destruirMatrizPesos(float**, int);
	void destruirMatrizCaminos(int**, int);

	int myrank = 0, numproc = 0, n = 0, lm = 0, resto = 0;
	int *chunk_sizes, *despl;
	float **dist;
	int **caminos;

	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &numproc);
	if (myrank == 0){
		if (numproc > maxnumprocs){
			printf("ERROR: superado el número máximo de procesos\n");
			MPI_Finalize();
			return 0;
		}
		if (argc > 1) sscanf(argv[1], "%i", &n);
		else{
			// TODO: argumento o input por teclado?? 
			printf("Usage: se debe proporcionar el número de vértices como el primer argumento\n");
			MPI_Finalize();
			return 0;
		}
		if (n>maxn){
			printf("ERROR: supero el número máximo de vértices\n");
			MPI_Finalize();
			return 0;
		}
		else printf("Número de vértices en el grafo: %d\n\n",n); 	
		resto = n % numproc;
		dist = crearMatrizPesos(n, n);
		caminos = crearMatrizCaminos(n, n);
		definirGrafo(n, dist, caminos);
		if (n<10) {
			printMatrizPesos(dist, n, n);
			// printMatrizCaminos(caminos, n, n);
		}
	} 
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	lm = n / numproc; 
	if(myrank != 0) {
		dist = crearMatrizPesos(n, n);
		caminos = crearMatrizCaminos(n, n);
		MPI_Scatter(&dist[resto+myrank+n%numproc][0], lm*n, MPI_FLOAT, &dist[resto+myrank+n%numproc][0], lm*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
	} else{
		MPI_Scatter(&dist[resto][0], lm*n, MPI_FLOAT, MPI_IN_PLACE, lm*n, MPI_FLOAT, 0,MPI_COMM_WORLD);
	}
	destruirMatrizPesos(dist, n);
	destruirMatrizCaminos(caminos, n);
	MPI_Finalize();
}

void definirGrafo(int n,float **dist,int **caminos)
{
	int i,j;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			if (i==j) dist[i][j]=0;
			else {
				dist[i][j]= (11.0 * rand() / ( RAND_MAX + 1.0 )); // aleatorios 0 <= dist < 11
				dist[i][j] = ((double)((int)(dist[i][j]*10)))/10; // truncamos a 1 decimal
			if (dist[i][j] < 2) dist[i][j]=0; // establecemos algunos a 0 
		}
		if (dist[i][j] != 0) caminos[i][j] = i+1;
		else caminos[i][j] = 0;
		}
	}
}

void printMatrizCaminos(int **a, int fila, int col) {
	int i, j;
	char buffer[10];
	printf("     ");
	for (i = 0; i < col; ++i){
		j=sprintf(buffer, "%c%d",'V',i+1 );
		printf("%5s", buffer);
	}
	printf("\n");
	for (i = 0; i < fila; ++i) {
		j=sprintf(buffer, "%c%d",'V',i+1 );
		printf("%5s", buffer);
		for (j = 0; j < col; ++j)
			printf("%5d", a[i][j]);
		printf("\n");
	}
	printf("\n");
}

void printMatrizPesos(float **a, int fila, int col) {
	int i, j;
	char buffer[10];
	printf("     ");
	for (i = 0; i < col; ++i){
		j=sprintf(buffer, "%c%d",'V',i+1 );
		printf("%5s", buffer);
	}
	printf("\n");
	for (i = 0; i < fila; ++i) {
		j=sprintf(buffer, "%c%d",'V',i+1 );
		printf("%5s", buffer);
		for (j = 0; j < col; ++j)
				printf("%5.1f", a[i][j]);
		printf("\n");
	}
	printf("\n");
}

float** crearMatrizPesos(int filas, int columnas){
	if (filas <= 0 || columnas <=0 ){
		printf("ERROR: incorrectas dimensiones de la matriz de pesos\n");
		return 0;
	}
	float** matriz = (float**)malloc(filas*sizeof(float*));
	matriz[0] = (float*)malloc(filas*columnas*sizeof(float));
	if (matriz==NULL){
		printf("ERROR: insuficiente espacio de memoria\n");
		return 0;
	}
	for (int i=1;i<filas;++i){
		matriz[i] = matriz[i-1] + columnas;
	}
	return matriz;
}

int** crearMatrizCaminos(int filas, int columnas){
	if (filas <= 0 || columnas <=0 ){
		printf("ERROR: incorrectas dimensiones de la matriz de pesos\n");
		return 0;
	}
	int** matriz = (int**)malloc(filas*sizeof(int*));
	matriz[0] = (int*)malloc(filas*columnas*sizeof(int));
	if (matriz==NULL){
		printf("ERROR: insuficiente espacio de memoria\n");
		return 0;
	}
	for (int i=1;i<filas;++i){
		matriz[i] = matriz[i-1] + columnas;
	}
	return matriz;
}

void destruirMatrizCaminos(int** caminos, int n){
	for (int i=1;i<n; i++){
		caminos[i] = NULL;
	}
	free(caminos[0]);
	free(caminos);
}

void destruirMatrizPesos(float** dist, int n){
	for (int i=1;i<n; i++){
		dist[i] = NULL;
	}
	free(dist[0]);
	free(dist);
}
