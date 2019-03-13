#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

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
	void calculaCamino(float**, int**, int);

	int myrank = 0, numproc = 0, n = 0, lm = 0, resto = 0;
	int *chunk_sizes, *despl;
	float **dist;
	int **caminos;

	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &numproc);
	if (numproc > maxnumprocs){
		if(myrank==0)printf("ERROR: superado el número máximo de procesos\n");
		MPI_Finalize();
		return 0;
	}
	if (myrank == 0){
		if (argc > 1) sscanf(argv[1], "%i", &n);
		else{
			printf("Introduce el número de vertices:\n");
			scanf("%i", &n);
		}
	} 
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (n>maxn){
		if (myrank == 0)printf("ERROR: supero el número máximo de vértices\n");
		MPI_Finalize();
		exit(0);
	} else if(n<=0){
		if (myrank == 0)printf("ERROR: incorrecto número de vértices\n");
		MPI_Finalize();
		exit(0);
	}
	else if (myrank==0) printf("Número de vértices en el grafo: %d\n\n",n); 	
	lm = n / numproc; 
	dist = crearMatrizPesos(n, n);
	caminos = crearMatrizCaminos(n, n);
	if (myrank == 0){
		definirGrafo(n, dist, caminos);
		if (n<10){
			printMatrizPesos(dist, n,n);
			printMatrizCaminos(caminos, n, n);
		}
		resto = n % numproc;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank != 0) {
		MPI_Scatter(&dist[0][0], lm*n, MPI_FLOAT, &dist[0][0], lm*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Scatter(&caminos[0][0], lm*n, MPI_INT, &caminos[0][0], lm*n, MPI_INT, 0, MPI_COMM_WORLD);
	} else{
		MPI_Scatter(&dist[resto][0], lm*n, MPI_FLOAT, MPI_IN_PLACE, lm*n, MPI_FLOAT, 0,MPI_COMM_WORLD);
		MPI_Scatter(&caminos[resto][0], lm*n, MPI_INT, MPI_IN_PLACE, lm*n, MPI_INT, 0,MPI_COMM_WORLD);
	}
	int sender = 0, fila_k = 0;
	float* auxd = (float*)malloc(n*sizeof(float));
	int* auxc = (int*)malloc(n*sizeof(int));
	for (int k=0;k<n;k++){
		sender = (k<lm+n%numproc) ? 0 : numproc-1-((n-k-1)/lm);
		fila_k = (sender==0) ? k : k - (lm*sender) - n % numproc;
		if (myrank==sender){
			memcpy(&auxd[0], &dist[fila_k][0], n*sizeof(float));
			memcpy(&auxc[0], &caminos[fila_k][0], n*sizeof(int));
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&auxd[0], n, MPI_FLOAT, sender, MPI_COMM_WORLD);
		MPI_Bcast(&auxc[0], n, MPI_INT, sender, MPI_COMM_WORLD);
		for (int i=0;i<lm+resto; i++){
			for (int j=0;j<n;j++){
				if (dist[i][k] * auxd[j] != 0){
					if ( (dist[i][k]+auxd[j] < dist[i][j]) || (dist[i][j]==0)){
						dist[i][j] = dist[i][k] + auxd[j];
						caminos[i][j] = auxc[j];
					}
				}
			}
		}
	}
	auxc = NULL;
	auxd = NULL;
	free(auxc);
	free(auxd);
	if (myrank != 0){
		MPI_Gather(&dist[0][0], lm*n, MPI_FLOAT, &dist[0][0], lm*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Gather(&caminos[resto][0], lm*n, MPI_INT, &caminos[resto][0], lm*n, MPI_INT, 0, MPI_COMM_WORLD);
	} else {
		MPI_Gather(MPI_IN_PLACE, lm*n, MPI_FLOAT, &dist[resto][0], lm*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Gather(MPI_IN_PLACE, lm*n, MPI_INT, &caminos[resto][0], lm*n, MPI_INT, 0, MPI_COMM_WORLD);
	}
	if (myrank == 0){
		if (n<10){
			printMatrizPesos(dist, n,n);
			printMatrizCaminos(caminos, n, n);
		}
		calculaCamino(dist, caminos, n);
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

void calculaCamino(float **a, int **b, int n)
{
	int i,count=2, count2;
	int anterior; 
	int *camino;
	int inicio=-1, fin=-1;
	while ((inicio < 0) || (inicio >n) || (fin < 0) || (fin > n)) {
		printf("Vertices inicio y final: (0 0 para salir)\n");
		scanf("%d %d",&inicio, &fin);
	}
	while ((inicio != 0) && (fin != 0)) {
		anterior = fin;
		while (b[inicio-1][anterior-1] != inicio) {
			anterior = b[inicio-1][anterior-1];
			count = count + 1;
		}
		count2 = count;
		camino = malloc(count * sizeof(int));
		anterior = fin;
		camino[count-1]=fin;
		while (b[inicio-1][anterior-1] != inicio) {
			anterior = b[inicio-1][anterior-1];
			count = count - 1;
			camino[count-1]=anterior;
		}
		camino[0] = inicio;
		printf("\nCamino mas corto de %d a %d:\n", inicio, fin);
		printf("\tPeso: %5.1f\n", a[inicio-1][fin-1]);
		printf("\tCamino: ");
		for (i=0; i<count2; i++) printf("%d ",camino[i]);
		printf("\n");
		free(camino);
		inicio = -1;
		while ((inicio < 0) || (inicio >n) || (fin < 0) || (fin > n)) {
			printf("Vertices inicio y final: (0 0 para salir)\n");
			scanf("%d %d",&inicio, &fin);
		}
	}
}