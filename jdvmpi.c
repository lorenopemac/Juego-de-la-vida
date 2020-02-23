#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "emmintrin.h"
#include <mpi.h>

#ifdef _OPENMP
	#include <omp.h>
	#define TRUE 1
	#define FALSE 0
#else
	#define omp_get_thread_num() 0
	#define omp_get_num_threads() 1
#endif

int main(int argc, char *argv[]){

	// Comienza Sección MPI -----------------------------------------------------------------------------

		int n, size, rank,filas,columnas,cols,rows,steps,index,index2,i,j,k,largo,perd,vivas;

		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Status status;
		MPI_Request request,request2;
		MPI_Comm_size(MPI_COMM_WORLD, &size);
                request = malloc(2 * (size-1) * sizeof(MPI_Request));
                request2 = malloc(2 * (size-1) * sizeof(MPI_Request));
                n = size;
                status  = malloc(2 * (size-1) * sizeof(MPI_Status));
                __m128i v1, v2, v3, v4,v5,aux1,aux2,aux3,aux4,final, a,b;
                char recvdimsubgrilla[2];
                char **ptr = NULL;
                char **ptrAux = NULL;	
                char **aux = NULL;
		if(rank == 0){
			
			int indiceFilas[n], indiceCol[n], arr1[n], arr2[n], arrResta[n], menor, scols, sfilas

			// Lectura de archivo
			FILE *fp;
			fp = fopen("entrada.cells","r+");  
                        char dimsubgrilla[2];
			fscanf(fp,"cols %d\n rows %d\n steps %d\n",&cols,&rows,&steps);
				printf("cols %d\nrows %d\nsteps %d\n",cols,rows,steps);	
			//	printf("-------------------------------------------------------------------\n");
			//asignacion de memoria	
			i=0;
			//cols += 17; 	//agrego 16 columnas mas (al principio) para poder trabajar con instrucciones vectoriales (por el 1er shitf-left)
							// y una mas a la izquierda. Con esto quedarian incluidos los dos bordes laterales que se suban en la version serial.
                        
                        for(i=2; i<=n/2; i++)
                        {
				if(n % i == 0){
					indiceFilas[cont] = i;
                                        indiceCol[cont] = n/indiceFilas[cont];
                                        cont=cont+1;
                                }
			}
			for(i=0; i<cont; i++)
                        {   
                                arr1[i] =rows/indiceFilas[i];
				arr2[i] = cols/indiceCol[i];
				arrResta[i] =fabs(arr1[i] - arr2[i]);
                        }
			menor = 0;	//menor almacena indices no valores, por lo tanto este valor inicial corresponde a la primera posicion de arrResta;

 			for(i=1; i<cont; i++){
				if(arrResta[i] < arrResta[menor])
					menor = i;
			}
                        sfilas = indiceFilas[menor];	//CANTIDAD DE FILAS EN LA SUB GRILLA
			scols = indiceCol[menor];    //CANTIDAD DE COLUMNAS EN LA SUB GRILLA
                        
			dimsubgrilla[0]=(char)sfilas;
                        dimsubgrilla[1]=(char)scols;
			for(i=1;i<n ;i++)
                            MPI_ISend(dimsubgrilla,2,MPI_CHAR,i,0,MPI_COMM_WORLD,&request2[i-1]);
                            MPI_Wait(&request2,&status);
			//asigno memoria al contigua para las "rows+2" filas de la matriz    
                        filas=sfilas;//establesco el mismo nombre para todos filas
                        columnas=scols;
                        /*CREACION MATRIZ PROPIA */
                        if(posix_memalign((void**)&ptr0, 16, sizeof(char*)*(filas)) != 0)
			  return 1;
                        for(index=0; index<sfilas; index++)
                        if(posix_memalign((void**)&ptr0[index], 16, columnas) != 0)
				  return 1;
                        /*FIN CREACION MATRIZ PROPIA*/ 
			
                       /* No es necesario reservar espacio para la matriz
                        * if(posix_memalign((void**)&ptr, 16, sizeof(char*)*(rows+2)) != 0)//NO ES NECESARIO RESERVAR TANTO ESPACIO
			  return 1;
			if(posix_memalign((void**)&ptrAux, 16, sizeof(char*)*(rows+2)) != 0)
			  return 1;

			if (!ptr){
				printf ("\nFALLO EN LA ASIGNACION DE MEMORIA !\n\n");
				exit (EXIT_FAILURE);
			}
			//asigno memoria contigua a cada columna, del tamaño de "cols+2"
			for(index=0 ; index<rows+2 ; index++){
				if(posix_memalign((void**)&ptr[index], 16, cols) != 0)
				  return 1;
				if(posix_memalign((void**)&ptrAux[index], 16, cols) != 0)
				  return 1;

				if (!*(ptr+index) || !*(ptrAux+index)) {	
						printf ("\nFALLO EN LA ASIGNACION DE MEMORIA !\n");
						exit (EXIT_FAILURE);
				}
			}*/

			//Se carga la matriz central con los valores del archivo

			for(index=1; index<(rows+1); index++)
					fgets((*(ptr+index)+16), cols, fp);	//comienza a cargarse la matriz desde la columna 16

			//Se imprime la matriz central cargada (testing)
		/*	printf("\nCarga inicial \n\n");
	
			for(index=1 ; index<(rows+1) ; index++){
				printf("%d: %s",index,*(ptr+index)+16);
			}	
		*/

			// Reemplazo Puntos y Ceros
			int indiceCore=0,core=0,cantColsCore;
                        cantColsCore=cols/scols;
			for(index=1; index<(rows+1); index++){//deberia cargar un buffer
				largo = 0;
				while(ptr[index][largo+16] != '\n')	//calcula el largo de la cadena
					largo++;
		
				for(index2=0; index2<largo; index2++){
			
					if(ptr[index][index2] == '.' ){
						ptr[index][index2] = 0;
					}
					else{
						ptr[index][index2] = 1;			
					}
				}	
				//completo las lineas faltantes/incompletas
				largo = 0;
				while(ptr[index][largo] != '\n')	//calcula el largo de la cadena
					largo++;

				if(largo+16 < cols){
					for(index2=largo; (index2)<cols-1; index2++){
						ptr[index][index2] = 0;
					}
				}
				for(i=0; i<cantColsCore; i++){//falta otro for por la cantidad de filas
                                   if(indiceCore+core==0){
                                            ptr0[index-1][j]=ptr[index][j];//ingresa valores en la matriz de pc 0
                                            core++;//aumenta el rank del core al que se le envia el msj
                                    }
                                    else
                                    {
                                            MPI_ISend(ptr[index][(scols*i)],scols,MPI_CHAR,core,0,MPI_COMM_WORLD,&request[i]);//ENVIA FILAS A LOS PROCESOS
                                            MPI_Wait(&request,&status);
                                            core++;//aumenta el rank del core al que se le envia el msj
                                    }
                                    if(core==indiceCore+cantColsore){
                                            core=indiceCore;//incrementa o reset el valor del core
                                    }
                                }
                                if((index+1) % sfilas == 0){
                                    indiceCore+= cantColsCore;
                                    core= indiceCore;
                                }
				// SEND -----------
                        }

			// imprimo matriz completada con ceros
		/*	printf("\nMatriz Completa \n\n");

			for(index=1; index<(rows+1); index++){
				printf("%d: ",index);		
	malloc			for(index2=16; index2<cols-1; index2++){
					printf("%d",(char)ptr[index][index2]);
				}
				printf("\n");
			}
		*/
		// FIN CARGA DE ARCHIVO -----------------------------------------------------
                }
                //MPI_Recv(MATRIZ, PREGUNTAR, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
	
	// FIN DE TAREA DEL PROCESO CON RANK 0 -----------------------------------------------------------		
        if(rank != 0){
            MPI_Recv(recvdimsubgrilla,2,MPI_CHAR,0,0,MPI_COMM_WORLD,&status); //RECEPCION DE LAS DIMENSIONES
            /*CREACION MATRIZ PROPIA */
            filas=recvdimsubgrilla[0];
            columnas=recvdimsubgrilla[1];
            if(posix_memalign((void**)&ptr0, 16, sizeof(char*)*(filas)) != 0)//CREACION MATRIZ PROPIA (ver sfilas+2)
                return 1;
            
            for(index=0; index<filas; index++){
                if(posix_memalign((void**)&ptr0[index], 16, sizeof(char*)*(columnas)) != 0)
                return 1;
            }
            
            /*FIN CREACION MATRIZ PROPIA */
            char *bufferfilas;// BUFFER AUXILIAR PARA ALMACENAR LAS FILAS ENVIADAS
            for(i=0; i<filas, i++){
                MPI_Recv(bufferfilas,columnas,MPI_CHAR,0,0,MPI_COMM_WORLD,&status);//RECEPCION DE CADA FILA ENVIADA POR EL PC 0
                for(j=0; j<columnas; j++)
                    ptr0[i][j]=bufferfilas[j];// INGRESO DE VALORES A LA MATRIZ PROPIA DE CADA PROCESO
            }
        }
	// Comienzan a correr las etapas del juego    
	/*
  	for (i=0; i<steps; i++) {
	
		//copio los bordes de la matriz
		//bordes laterales 
		/*#pragma omp for nowait
		for(index=1; index<(rows+1) ; index++){
			*(*(ptr+index)+15) = *(*(ptr+index)+cols-2);
			*(*(ptr+index)+(cols-1)) = *(*(ptr+index)+16);
		}
		
		//bordes superior e inferior
		#pragma omp parallel for
		for(index=16; index<(cols-1); index++){
			*(*(ptr)+index) = *(*(ptr+rows)+index);
			*(*(ptr+rows+1)+index) = *(*(ptr+1)+index);
		}

		// Intercambio valores de puntas
		*(*(ptr)+15) = *(*(ptr+rows)+(cols-2));
		*(*(ptr)+(cols-1)) = *(*(ptr+rows)+16);
		
		*(*(ptr+rows+1)+15) = *(*(ptr+1)+(cols-2));
		*(*(ptr+rows+1)+(cols-1)) = *(*(ptr+1)+16);

		printf("\nMatriz con bordes Completa \n\n");

		for(index=0; index<(rows+2); index++){
			printf("%d: ",index);		
			for(index2=15; index2<cols; index2++){
				printf("%d",(char)ptr[index][index2]);
			}
			printf("\n");
		}


		// Recorro filas
		#pragma omp parallel for private (k,v1,v2,v3,v4,v5,perd,aux2,aux3,aux4,a,b,final) schedule(static)
		for(j=0; j<rows ; j++){		

			// Recorro columnas
			for(k=16; k<cols ; k+=16){

				v1 = _mm_load_si128((__m128i*)&ptr[j][k]); 
				v2 = _mm_load_si128((__m128i*)&ptr[j+1][k]);
				v3 = _mm_load_si128((__m128i*)&ptr[j+2][k]);
			
				v4 = _mm_add_epi8(v1,v2);
				v4 = _mm_add_epi8(v3,v4); 	//v4 almacena la suma de cada columna de v1,v2,v3
				
				v1 = _mm_slli_si128(v4,1); //v1 almacena shift-left de v4
				v3 = _mm_srli_si128(v4,1); //v3 almacena shift-right de v4
				
				perd = ptr[j][k-1] +  ptr[j+1][k-1] +  ptr[j+2][k-1];	//suma almacena la suma vertical del borde izquierdo perdido
				v5 = _mm_set_epi8 (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,perd);
				v1 = _mm_add_epi8(v1,v5);
			
				perd = ptr[j][k+16] +  ptr[j+1][k+16] +  ptr[j+2][k+16]; //suma almacena la suma vertical del borde derecho perdido
				v5 = _mm_set_epi8 (perd,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
				v3 = _mm_add_epi8(v3,v5); 

				v5= _mm_add_epi8(v1,v4);
				v1= _mm_add_epi8(v3,v5); // ahora v1 almacena v1+v3+v4, ie, la suma de los tres vectores SIN PERDIDA

				v1= _mm_sub_epi8(v1,v2); //resto el casillero evaluado para que en v1 quede la suma de los 8 vecinos

				aux4=_mm_set1_epi8(1);

				a=_mm_set1_epi8(3);
				aux2=_mm_cmpeq_epi8(v1,a); 	//para c/elto si v1=='3' ...
				aux2=_mm_and_si128(aux2,aux4);	// ...entonces aux2 es 1 (viva) en esa celda

				a=_mm_set1_epi8(2);
				b=_mm_cmpeq_epi8(v1,a); 	//para c/elto si v1=='2' ...
				b=_mm_and_si128(b,aux4);	// ...entonces en r2 el valor se mantiene
				aux3=_mm_and_si128(v2,b); 
			
				final =_mm_or_si128(aux2,aux3);
				_mm_store_si128((__m128i*)&(ptrAux[j+1][k]), final);			
			}
		}
		//intercambio las matrices
		aux = ptr;
		ptr = ptrAux;	
		ptrAux = aux;
	}

	//calcula las celulas vivas en la ultima etapa
	vivas = 0;
	#pragma omp parallel for private(j) reduction(+:vivas)
	for(i=1;i<rows+1;i++){
		for(j=16;j<cols-1;j++){
			vivas += ptr[i][j];
		} 
	}	
//    printf("\nFIN DEL JUEGO, CELULAS VIVAS: %d\n",vivas);
		
//	printf("MATRIZ FINAL: \n");
	// imprimo matriz final
	for(i=1;i<rows+1;i++){
		for(j=16;j<cols-1;j++){
			if(ptr[i][j] == 0)			
				printf("%c", '.');
			else printf("%c", 'O');
		} 
		printf("\n");
	}
	printf("\n");
	printf("Number of live cells = %d",vivas);
*/
                /*
        if(rank!=0){
            for(i=0;i<recvdimsubgrilla[0];i++){
                MPI_Isend(ptr0[i],recvdimsubgrilla[1],MPI_CHAR,0,0,MPI_COMM_WORLD,&request[i]);
                MPI_Wait(&request,&status);
            }
        }
        if(rank==0)
        {
            for(i=1;i<rows;i++){//FILAS TOTALES MATRIZ ORIGINAL
                for(j=1;j<scols; j++)
                    printf("%c", '.');
                    
                    for(index=1; index<n; index++){
			recive
			printf
                    }
		printf("\n");
                //FALTA OTRO  PARA CUANDO EL 0 NO IMPRIME. 
            }
        }*/
        MPI_Finalize();
	return 0;
}
