#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define START_TAG 1
#define STOP_TAG 99

typedef struct Min_Max{
	double Min;
	double Max;
} Min_Max;

int main(int argc, char** argv){
	int** crearMatriz(int, int);
	void destruirMatriz(int**, int);
	Min_Max buscarMin_Max(int**, int, int);
	double ctimer(void);

	int pixelX, pixelY;
	int pixelXmax = 1024, pixelYmax = 1024;
	double Creal, Cimg;
	const double RealMinArray[] = {-2, -1.023438, -1.017579, -1.017523, -1.0190739863281251, -1.0184326403503419, -1.0175024721407624, -1.0176950869990224, 0.2720};
	const double RealMaxArray[] = { 1, -0.992188, -1.016968, -1.017493, -1.0178532832031251, -1.0184266798858643, -1.0175010058651026, -1.0173896129032258, 0.3720};
	const double ImMinArray[] = {-1.5, -0.285156, -0.274444, -0.274065, -0.2672003476562500, -0.2667537630310058,  -0.2740544516129032, -0.2772175483870968, 0.4805};
	const double ImMaxArray[] = { 1.5, -0.25,     -0.273758, -0.274032, -0.2658270664062500, -0.2667336466064453,  -0.2740528387096774, -0.2768738924731183, 0.5805}; 
	double RealMin=-2.0; 
	double RealMax=1.0;  
	double ImMin=-1.5; 
	double ImMax=1.5; 
	double AnchoPixel;
	double AltoPixel;
	double Tinicial, Tfinal, Ttotal;

	const int MaxValorTonos=255; 
	FILE *ImgFile, *ImgFile2;
	char ArchivoImagen[]="imgA.pgm";
	char ArchivoImagen2[]="imgB.pgm";
	char *comentario="# ";
	int bn, bn2;
	double Zx, Zy;             // Z=Zx+Zy*i  
	double Zx2, Zy2;           // Zx2=Zx*Zx, Zy2=Zy*Zy  
	int Iter;
	int IterMax=1000;
	const double Salida=2;    // valor de escape
	double Salida2=Salida*Salida, SumaExponencial;
	int dominio=0;
	int **matriz, **matriz2;
	Min_Max  min_max_img;

	int *buffer; 
	int sizeBuffer;


	int myrank = 0, numproc = 0;
	MPI_Status status;
	MPI_Request request;

	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &numproc);

	if (myrank == 0){


		printf("\n  *************** POSIBLES PARÁMETROS ***************\n");
		printf("  * Ancho imagen: pixelXmax                         *\n");
		printf("  * Alto imagen: pixelYmax                          *\n");
		printf("  * Numero de iteraciones: IterMax                  *\n");
		printf("  * Nombre archivo imagen de salida: ArchivoImagen  *\n");
		printf("  * Intervalo a considerar: dominio                 *\n");  
		printf("  ***************************************************\n\n");

		// Lectura de datos
		switch (argc){
			case 6: // dominio
				sscanf(argv[5], "%i", &dominio);
				if (8*(dominio + 1) > sizeof(RealMinArray)){
					printf("El array de dominios no es tan grande\n");
					return 0;
				}
				RealMin = RealMinArray[dominio];
				RealMax = RealMaxArray[dominio];
				ImMin = ImMinArray[dominio];
				ImMax = ImMaxArray[dominio];
			case 5: // nombre de los archivos
				strcpy(ArchivoImagen, argv[4]);
				strcat(ArchivoImagen, "A.pgm");
				strcpy(ArchivoImagen2, argv[4]);
				strcat(ArchivoImagen2, "B.pgm");
			case 4: // numero maximo de iteraciones
				sscanf(argv[3], "%i", &IterMax);
			case 3: // alto de la imagen
				sscanf(argv[2], "%i", &pixelYmax);
			case 2: // ancho de la imagen
				sscanf(argv[1], "%i", &pixelXmax);
			case 1:	
				break;
			default: 
				printf("Demasiados parametros\n");
				return 0;
		}

		printf("  *************** DATOS DE LA EJECUCIÓN ****************************\n");
		printf("  * pixelXmax = %4d, pixelYmax = %4d                             *\n", pixelXmax, pixelYmax);
		printf("  * IterMax = %5d                                                *\n", IterMax);
		printf("  * ArchivoImagen = %15s %15s                *\n", ArchivoImagen, ArchivoImagen2);
		printf("  * Dominio = %2d                                                   *\n", dominio);
		printf("  * Intervalo para X: [%20.15f,%20.15f]  *\n  * Intervalo para Y: [%20.15f,%20.15f]  *\n", RealMin, RealMax, ImMin, ImMax);
		printf("  ******************************************************************\n\n");

		/* Se crea un nuevo archivo y se abre en binario  */
		ImgFile= fopen(ArchivoImagen,"wb"); 
		ImgFile2= fopen(ArchivoImagen2,"wb"); 
		/* Se escribe la cabecera - ASCII  */
		fprintf(ImgFile,"P5\n %s\n %d\n %d\n %d\n",comentario,pixelXmax,pixelYmax,MaxValorTonos);
		fprintf(ImgFile2,"P5\n %s\n %d\n %d\n %d\n",comentario,pixelXmax,pixelYmax,MaxValorTonos);

		// las matrices donde se realizaran los calculos
		matriz = crearMatriz(pixelYmax, pixelXmax);
		matriz2 = crearMatriz(pixelYmax, pixelXmax);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// bcast de las variables necesarias para los calculos

	MPI_Bcast(&RealMin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&RealMax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ImMin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ImMax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&pixelXmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&pixelYmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&IterMax, 1, MPI_INT, 0, MPI_COMM_WORLD);

	AnchoPixel=(RealMax-RealMin)/(pixelXmax-1); 
	AltoPixel=(ImMax-ImMin)/(pixelYmax-1);

	// matriz de datos con la que van a trabajar los hijos:
	int* datos = (int*) malloc(pixelXmax*2*sizeof(int));

	// int sizeBuffer;
	// int* buffer;
	MPI_Pack_size(pixelXmax*2, MPI_INT, MPI_COMM_WORLD, &sizeBuffer);
	sizeBuffer = numproc*(sizeBuffer + MPI_BSEND_OVERHEAD);
	buffer = (int*) malloc(sizeBuffer);
	if (buffer == NULL)
		printf("Error al reservar la memoria del buffer\n");
	MPI_Buffer_attach(buffer, sizeBuffer);

	if(myrank==0){
		double Tinicial, Tfinal, Ttotal;
		int fila = 0;


		// PROCESO 0 
		Tinicial = MPI_Wtime();

		// enviar una fila a cada proceso
		for (int i=1;i<numproc; ++i){
			MPI_Send(&fila, 1, MPI_INT, i, START_TAG, MPI_COMM_WORLD);
			fila++;
		}

		// asignar nuevas tareas a procesos inactivos hasta que se acaben
		while(fila < pixelYmax){
			MPI_Recv(&datos[0], 2*pixelXmax, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			pixelY = status.MPI_TAG;
			int source = status.MPI_SOURCE;
			MPI_Send(&fila, 1, MPI_INT, source, START_TAG, MPI_COMM_WORLD);
			fila++;
			memcpy(&matriz[pixelY][0], &datos[0], pixelXmax*sizeof(int));
			memcpy(&matriz2[pixelY][0], &datos[pixelXmax], pixelXmax*sizeof(int));
		}

		// recibir las filas restantes
		for (int i=1;i<numproc; ++i){
			MPI_Recv(&datos[0], 2*pixelXmax, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			pixelY = status.MPI_TAG;
			memcpy(&matriz[pixelY][0], &datos[0], pixelXmax*sizeof(int));
			memcpy(&matriz2[pixelY][0], &datos[pixelXmax], pixelXmax*sizeof(int));
		}

		// indicar a los hijos que finalicen la ejecución
		for (int i=1;i<numproc; ++i){
			MPI_Send(&fila, 1, MPI_INT, i, STOP_TAG, MPI_COMM_WORLD);
		}
		
		Tfinal = MPI_Wtime();
		Ttotal = Tfinal - Tinicial;
		printf("\nTiempo: %f segundos\n", Ttotal);


		min_max_img = buscarMin_Max(matriz2,pixelYmax,pixelXmax);
		for (int i=0;i<pixelYmax;i++) {
			for (int j=0;j<pixelXmax;j++) {
				matriz2[i][j] = matriz2[i][j] - min_max_img.Min;
				matriz2[i][j] = matriz2[i][j] * (255.0/(min_max_img.Max-min_max_img.Min));
			}
		}
		min_max_img = buscarMin_Max(matriz2,pixelYmax,pixelXmax);

		for (pixelY=0;pixelY<pixelYmax;pixelY++) {
			for(pixelX=0;pixelX<pixelXmax;pixelX++){
				fwrite(&matriz[pixelY][pixelX],1,1,ImgFile);
				fwrite(&matriz2[pixelY][pixelX],1,1,ImgFile2);
			}
		} 
		fclose(ImgFile);
		fclose(ImgFile2);
	} else {

		// PROCESOS HIJOS
		while(1){
			// recibir fila y comprobar la etiqueta
			MPI_Recv(&pixelY, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if (status.MPI_TAG == STOP_TAG){
				MPI_Finalize();
				return 0;
			}

			// realizar calculos
			Cimg = ImMin + pixelY*AltoPixel;
			for (pixelX=0;pixelX<pixelXmax;pixelX++){
				Creal = RealMin + pixelX*AnchoPixel;
				Zx=0.0;         // Valor inicial
				Zy=0.0;
				Zx2=Zx*Zx;
				Zy2=Zy*Zy;
				SumaExponencial=0;
				for (Iter=0;Iter<IterMax && ((Zx2+Zy2)<Salida2);Iter++)
				{
					SumaExponencial += exp( -sqrt(Zx2+Zy2) )/IterMax;   // se mantiene siempre entre (0,1)            
					Zy=2*Zx*Zy + Cimg;
					Zx=Zx2-Zy2 + Creal;
					Zx2=Zx*Zx;
					Zy2=Zy*Zy;
				}
				if (Iter==IterMax) { /*  interior del conjunto Mandelbrot = negro */    
					bn=0;          
					bn2=0;            
				}
				else { /* exterior del conjunto Mandelbrot = blanco modificado con exp */
					bn = MaxValorTonos - SumaExponencial * 255; 
					bn2 = (Iter +1 - (int)(log(log(sqrt(Zx2+Zy2)) / log(2)) / log(2)))*255;
				}
				datos[pixelX] = bn;
				datos[pixelXmax+pixelX] = bn2;
			}

			// enviar resultado al padre
			MPI_Bsend(&datos[0], 2*pixelXmax, MPI_INT, 0, pixelY, MPI_COMM_WORLD);
		}

	}
	MPI_Buffer_detach(buffer, &sizeBuffer);
	if(myrank == 0){
		destruirMatriz(matriz, pixelYmax);
		destruirMatriz(matriz2, pixelYmax);
	}
	free(datos);
	free(buffer);
	MPI_Finalize();
}

int** crearMatriz(int filas, int columnas){
	if (filas <= 0 || columnas <=0 ){
		printf("ERROR: incorrectas dimensiones de la matriz\n");
		return 0;
	}
	int** matriz = (int**)malloc(filas*sizeof(int*));
	if (matriz==NULL){
		printf("ERROR: problemas al dimensionar\n");
		return 0;
	}
	for (int i=0;i<filas;++i){
		matriz[i] = (int*)malloc(columnas*sizeof(int));
		if (matriz==NULL){
			printf("ERROR: problemas al dimensionar\n");
			return 0;
		}
	}
	return matriz;
}

void destruirMatriz(int** matriz, int filas){
	for (int i=0;i<filas;i++)
		free(matriz[i]);
	free(matriz);
}

Min_Max buscarMin_Max(int** array, int fil, int col){
	Min_Max  min_max;
	int index,fila,columna,fila_next,columna_next;
	int n = fil*col; 
	if ( n%2 != 0 ){

		min_max.Min = array[0][0];
		min_max.Max = array[0][0];

		index = 1;
	}
	else{
		if ( array[0][0] < array[0][1] ){
			min_max.Min = array[0][0];
			min_max.Max = array[0][1];
		}
		else {
			min_max.Min = array[0][1];
			min_max.Max = array[0][0];
		}
		index = 2;
	}

	int big, small,i;
	for ( i = index; i < n-1; i = i+2 ){
		fila = i / col;
		columna = i % col;  
		fila_next = fila;
		columna_next = columna + 1;
		if (columna_next == col) {
			fila_next = fila +1;
			columna_next = 0;
		}
		if ( array[fila][columna] < array[fila_next][columna_next] ){
			small = array[fila][columna];
			big = array[fila_next][columna_next];
		}
		else{
			small = array[fila_next][columna_next];
			big = array[fila][columna];
		}
		if ( min_max.Min > small ){
			min_max.Min = small;
		}
		if ( min_max.Max < big ){ 
			min_max.Max = big;
		}
	}
	printf("Minimo = %f, Maximo = %f\n", min_max.Min, min_max.Max);
	return min_max;
}