/*FALTA ENCONTRAR A LOS VECINOS DE UN PROCESO Y TERMINAR EL PROCESAMIENTO INTERNO*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "emmintrin.h"
#include <mpi.h>

int main(int argc, char *argv[]){

	// Comienza Sección MPI -----------------------------------------------------------------------------

        int n,q,rank,filas,size,columnas,cols,rows,steps,index,index2,i,j,k,largo,perd,vivas,indiceCore=0,core=0,cantColsCore,scols, sfilas;
        int aIzq,aMed,aDer,der,izq,bIzq,bMed,bDer,sumaAnte,sumaPost,total;

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD,&size);
        MPI_Request *request,reqs2[16];
        MPI_Status *status;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        n = size;
        char hostname[20];
        status  = malloc((size) * sizeof(MPI_Status));
        request  = malloc((size) * sizeof(MPI_Request));
        gethostname(hostname, sizeof(hostname));
        __m128i v1, v2, v3, v4,v5,v6,v7,v8,v9,aux1,aux2,aux3,aux4,final, a,b;
        int *dimsubgrilla=NULL, *recvdimsubgrilla=NULL;
        char *ptr = NULL;
        char **ptr0=NULL;
        char *ptr1=NULL;
        char *bufferfilas=NULL;// BUFFER AUXILIAR PARA ALMACENAR LAS FILAS ENVIADAS
        char **ptrAux = NULL;	
        char **aux = NULL;
        char *buffer=NULL;
        char *buffer1=NULL;
        char *resultado=NULL;
        //printf("proceso: %d, pid: %d, host: %s\n",rank, (int)getpid(), hostname);
        recvdimsubgrilla= malloc(sizeof(int)*5);    //buffer donde cada proceso recibe datos de las subgrillas
        //sleep(10);
        if(rank==0){
                        char indiceFilas[n],cont=0, indiceCol[n], arr1[n], arr2[n], arrResta[n], menor;

			// Lectura de archivo
			FILE *fp;
			fp = fopen("salida.cells","r+");  

			dimsubgrilla= malloc(sizeof(int)*5);
			fscanf(fp,"cols %d\n rows %d\n steps %d\n",&cols,&rows,&steps);
				printf("cols %d\nrows %d\nsteps %d\n",cols,rows,steps);
			i=0;
                        ptr= malloc(sizeof(char)*cols);
                        //printf("   {%d}   ",sizeof(ptr));
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
                                arr1[i] = rows/indiceFilas[i];
				arr2[i] = cols/indiceCol[i];
				arrResta[i] = fabs(arr1[i] - arr2[i]);
                                //printf("[%d _ %d] %d",indiceFilas[i],indiceCol[i],arrResta[i]);        
			}
			menor = 0;	//menor almacena indices no valores, por lo tanto este valor inicial corresponde a la primera posicion de arrResta;

 			for(i=1; i<cont; i++){
				if(arrResta[i] <= arrResta[menor])
					menor = i;
			}
                        sfilas = indiceFilas[menor];	//CANTIDAD DE FILAS EN LA SUB GRILLA
			scols = indiceCol[menor];    //CANTIDAD DE COLUMNAS EN LA SUB GRILLA
                        filas=rows/sfilas;//establesco el mismo nombre para todos los procesos 'filas' 
                        columnas=cols/scols;
			dimsubgrilla[0]=filas;//cantidad de filas por subgrilla
                        dimsubgrilla[1]=columnas;//cantidad de columnas por subgrilla
                        dimsubgrilla[2]=steps;
			dimsubgrilla[3]=sfilas;
			dimsubgrilla[4]=scols;
                        
			//printf("filas: %d columnas: %d",filas,columnas);
                        ptr1=malloc(sizeof(char)*columnas);
			for(i=1;i<n ;i++){
				  MPI_Send(dimsubgrilla,5,MPI_INT,i,0,MPI_COMM_WORLD);
                        }
                        /*CREACION MATRIZ PROPIA */
                        if(posix_memalign((void**)&ptr0, 16, sizeof(char*)*(filas)) != 0)
			  		return 1;
                        //printf("Tamaño: %d",rows/filas);
                        //printf("\n");
                        for(index=0; index<(rows/sfilas); index++)
                            if(posix_memalign((void**)&ptr0[index], 16, columnas+16) != 0)
                                        return 1;
                        /*FIN CREACION MATRIZ PROPIA*/
                        //se lee una linea de archivo, se traduce a 0 o 1 y se la envia a los demas procesos
                        printf("\n");
			buffer=malloc((cols+1)*sizeof(char));
			resultado=malloc(cols*sizeof(char));
			// Ingreso de valores a la matriz
                        resultado=fgets(buffer,cols,fp);
			index=0;
                        while(resultado!=NULL){
                            for(q=0;q<cols;q++){
                                    if(buffer[q]=='O'){
                                            ptr[q]=1;	
                                    }
                                    else{
                                            ptr[q]=0;
                                    }
                            }
                            for(i=0; i<scols; i++){
                                if(indiceCore+core==0){
				    for(j=0;j<columnas;j++){
                                        ptr0[index][j+16]=ptr[j];//ingresa valores en la matriz de pc
                                    }
                                    core++;//aumenta el rank del core al que se le envia el msj
                                }
                                else{
                                     if(core < size){
                                          for(k=0;k<columnas;k++){
                                                ptr1[k]=ptr[columnas*i+k];//buffer auxiliar
                                          }
                                           MPI_Send(ptr1,sizeof(char)*columnas,MPI_CHAR,core,0,MPI_COMM_WORLD);//ENVIA FILAS A LOS PROCESOS 
                                      }
                                      core++;//aumenta el rank del core al que se le envia el ms
                                }
		 		
                                if(core==indiceCore+scols){
                                    core=indiceCore;//incrementa o reset el valor del core
                                }
                            }
                             if((index+1) % (filas) == 0){
                                indiceCore+= scols;
                                core = indiceCore;
                             }
                            if(!feof(fp)){
                                memset(buffer,'\0',cols);
                                resultado=fgets(buffer,cols,fp);
                            }
                            index=index+1;
                           
                        }
                        fclose(fp);
                   //printf("Entro FIN LECTURA");
    }
                
    if(rank!=0){
        MPI_Recv(recvdimsubgrilla,5,MPI_INT,0,0,MPI_COMM_WORLD,status); //RECEPCION DE LAS DIMENSIONES
        //printf("Entro 2 RANK:%d",rank);
        //CREACION MATRIZ PROPIA 
        filas=recvdimsubgrilla[0];
        columnas=recvdimsubgrilla[1];
        steps=recvdimsubgrilla[2];
        sfilas=recvdimsubgrilla[3];
        scols=recvdimsubgrilla[4];
        
       //Se reserva espacio para almacenar cada subgrilla
        if(posix_memalign((void**)&ptr0, 16, sizeof(char*)*(filas)) != 0)//CREACION MATRIZ PROPIA (ver sfilas+2)
                return 1;
        for(index=0; index<filas; index++){
            if(posix_memalign((void**)&ptr0[index], 16, sizeof(char*)*(columnas+17)) != 0)
                return 1;
        }
        bufferfilas=malloc(columnas*sizeof(char));
        for(i=0; i<filas; i++){
            MPI_Recv(bufferfilas,columnas,MPI_CHAR,0,0,MPI_COMM_WORLD,status);//RECEPCION DE CADA FILA ENVIADA POR EL PC 0
            for(j=0; j<columnas; j++){
                ptr0[i][j+16]=bufferfilas[j];// INGRESO DE VALORES A LA MATRIZ PROPIA DE CADA PROCESO
            }
        }
                //FIN CREACION MATRIZ PROPIA 
    }
    
    if(rank>=(sfilas-1)*scols){
        bMed= rank % scols;
        if((rank % scols)== 0){
            bIzq= bMed+(scols-1);
            bDer=bMed+1;
        }
        else{
            if(((rank+1) % scols)==0){
                bIzq= bMed-1;
                bDer= bMed-(scols-1);
            }
            else{
                bIzq = bMed-1;
                bDer = bMed+1;
            }
        }
    }
    else{
        bMed= rank+scols;
        if(((rank+1)%scols)==0){
            bDer= bMed-(scols-1);
            bIzq= bMed-1; 
        }
        else{
            if(rank % scols == 0){
                bDer=bMed+1;
                bIzq=bMed+(scols-1);
            }
            else{
                bDer= bMed+1;
                bIzq= bMed-1;
            }   
        }
    }
				
//LOS DE ARRIBA
    if(rank<scols){
        aMed= scols*(sfilas-1)+rank;
        if(rank == 0){
            aDer= aMed+1;
            aIzq= aMed+(scols-1);
        }
        else{
            if(rank == scols-1){
                aDer=aMed-(scols-1);
                aIzq=aMed-1;
            }
            else{
                aIzq = aMed-1;
                aDer = aMed+1;
            }
        }
    }
    else{
        aMed= rank-scols;
        if(((rank+1)%scols)==0){
            aDer= aMed-(scols-1);
            aIzq= aMed-1;
        }
        else{
            if(rank%scols == 0){
                aDer=aMed+1;
                aIzq=aMed+(scols-1);
            }
            else{
                aDer=aMed+1;
                aIzq=aMed-1;
            }
        }
    }
    //FIN LOS DE ARRIBA

    //Izquierda  Derecha
    if(((rank+1)%scols)==0){
        izq=rank-1;
        der=rank-(scols-1);
    }
    else{
        if(rank%scols == 0){
            izq= rank+(scols-1);
            der= rank+1;
        }
        else{
            izq=rank-1;
            der=rank+1;
        }
    }
                //Fin Izquierda Derecha
//DEFINCION DE BUFFERS            
    /*if(rank==1){
        printf("\n");
        printf("rank=%d, izq=%d, der=%d, arribaDer=%d, arribaIzq=%d, arriba=%d,  abajoIzq=%d, abajoDer=%d, abajo=%d ",rank,izq,der,aDer,aIzq,aMed,bIzq,bDer,bMed);
        printf("\n");
    }
    else{sleep(2);}*/
    char *bordeInfEnvia=NULL;
    char *bordeSupEnvia=NULL;
    char *bordeInfIzqEnvia=NULL;
    char *bordeSupDerEnvia=NULL;
    char *bordeInfDerEnvia=NULL;
    char *bordeSupIzqEnvia=NULL;
    char *bordeIzqEnvia= NULL;
    char *bordeDerEnvia= NULL;
    char *bordeIzq=NULL;
    char *bordeDer=NULL;
    char *bordeInf=NULL;
    char *bordeSup=NULL;
    char *bordeSupIzq=NULL;
    char *bordeSupDer=NULL;
    char *bordeInfIzq=NULL;
    char *bordeInfDer=NULL;
                
    bordeSupEnvia= malloc(sizeof(char)*columnas);
    bordeInfEnvia= malloc(sizeof(char)*columnas);	
    bordeIzqEnvia= malloc(sizeof(char)*filas);
    bordeDerEnvia= malloc(sizeof(char)*filas);
    bordeIzq= malloc(sizeof(char)*filas);
    bordeDer= malloc(sizeof(char)*filas);
    bordeSup= malloc(sizeof(char)*columnas);
    bordeInf= malloc(sizeof(char)*columnas);
    bordeSupIzq= malloc(sizeof(char));
    bordeSupDer= malloc(sizeof(char));
    bordeInfIzq= malloc(sizeof(char));
    bordeInfDer= malloc(sizeof(char));
    bordeInfIzqEnvia= malloc(sizeof(char));
    bordeSupDerEnvia= malloc(sizeof(char));
    bordeInfDerEnvia= malloc(sizeof(char));
    bordeSupIzqEnvia= malloc(sizeof(char));
                
    if(posix_memalign((void**)&ptrAux, 16, sizeof(char*)*(filas)) != 0)//CREACION MATRIZ PROPIA (ver sfilas+2)
            return 1;
            
    for(index=0; index<filas; index++){
        if(posix_memalign((void**)&ptrAux[index], 16, sizeof(char*)*(columnas+17)) != 0)
            return 1;
    }
//FIN DEFINCION DE BUFFERS                
    
           
//EMPIEZA PROCESAMIENTO INTERNO
    for (i=0; i<steps; i++) {
            
//eviar bordes
        //borde superior
        for(index=0;index<columnas;index++)
            bordeSupEnvia[index]=ptr0[0][index+16];
        MPI_Isend(bordeSupEnvia,sizeof(char)*columnas,MPI_CHAR,aMed,1,MPI_COMM_WORLD,&reqs2[0]);
        
        //borde inferior
        for(index=0;index<columnas;index++)
            bordeInfEnvia[index]=ptr0[filas-1][index+16];
        MPI_Isend(bordeInfEnvia,sizeof(char)*columnas,MPI_CHAR,bMed,2,MPI_COMM_WORLD,&reqs2[1]);
        
        //borde izquierdo
        for(index=0;index<filas;index++)
            bordeIzqEnvia[index]=ptr0[index][16+columnas-1];
        MPI_Isend(bordeIzqEnvia,sizeof(char)*filas,MPI_CHAR,der,3,MPI_COMM_WORLD,&reqs2[2]);
        
        //borde derecho
        for(index=0;index<filas;index++)
            bordeDerEnvia[index]=ptr0[index][16];
        MPI_Isend(bordeDerEnvia,sizeof(char)*filas,MPI_CHAR,izq,4,MPI_COMM_WORLD,&reqs2[3]);
        
        //borde superior Izquierdo
        bordeSupIzqEnvia[0]=ptr0[filas-1][15+columnas];
        MPI_Isend(bordeSupIzqEnvia,sizeof(char),MPI_CHAR,aIzq,5,MPI_COMM_WORLD,&reqs2[4]);
        //borde superior Derecho
        bordeSupDerEnvia[0]=ptr0[filas-1][16];
        MPI_Isend(bordeSupDerEnvia,sizeof(char),MPI_CHAR,aDer,6,MPI_COMM_WORLD,&reqs2[5]);
        //borde inferior Izquierdo
        bordeInfIzqEnvia[0]=ptr0[0][16];
        MPI_Isend(bordeInfIzqEnvia,sizeof(char),MPI_CHAR,bIzq,7,MPI_COMM_WORLD,&reqs2[6]);
        //borde inferior Derecho
        bordeInfDerEnvia[0]=ptr0[0][15+columnas];
        MPI_Isend(bordeInfDerEnvia,sizeof(char),MPI_CHAR,bDer,8,MPI_COMM_WORLD,&reqs2[7]);
//Fin enviar bordes
        MPI_Irecv(bordeInf,sizeof(char)*columnas,MPI_CHAR,bMed,1,MPI_COMM_WORLD,&reqs2[8]);
    //    if(rank==0){printf("\n");printf("Envio borde bordeSupEnvia");printf("\n");}
        //borde Superior
        MPI_Irecv(bordeSup,sizeof(char)*columnas,MPI_CHAR,aMed,2,MPI_COMM_WORLD,&reqs2[9]);
    //    if(rank==0){printf("Envio borde bordeInfEnvia");printf("\n");}
        //borde Izquierdo
        MPI_Irecv(bordeIzq,sizeof(char)*filas,MPI_CHAR,izq,3,MPI_COMM_WORLD,&reqs2[10]);
    //    if(rank==0){printf("Envio borde bordeIzqEnvia");printf("\n");}
        //borde Derecho
        MPI_Irecv(bordeDer,sizeof(char)*filas,MPI_CHAR,der,4,MPI_COMM_WORLD,&reqs2[11]);
    //    if(rank==0){printf("Envio borde bordeDerEnvia");printf("\n");}
        //borde Superior Izquierdo
        MPI_Irecv(bordeInfDer,sizeof(char),MPI_CHAR,bDer,5,MPI_COMM_WORLD,&reqs2[12]);
    //    if(rank==0){printf("Envio borde bordeSupIzqEnvia");printf("\n");}
        //borde Superior Derecho
        MPI_Irecv(bordeInfIzq,sizeof(char),MPI_CHAR,bIzq,6,MPI_COMM_WORLD,&reqs2[13]);
    //    if(rank==0){printf("Envio borde bordeSupDerEnvia");printf("\n");}
        //borde Inferior Izquierdo
        MPI_Irecv(bordeSupDer,sizeof(char),MPI_CHAR,aDer,7,MPI_COMM_WORLD,&reqs2[14]);
    //    if(rank==0){printf("Envio borde bordeInfIzqEnvia");printf("\n");}
        //borde Inferior Derecho
        MPI_Irecv(bordeSupIzq,sizeof(char),MPI_CHAR,aIzq,8,MPI_COMM_WORLD,&reqs2[15]);
    //    if(rank==0){printf("Envio borde bordeInfDerEnvia");printf("\n");}
        

        // Recorro filas
        for(j=1; j<filas-1 ; j++){	
            // Recorro columnas
            for(k=16; k<columnas+16 ; k+=16){
                sumaAnte=ptr0[j][k-1] +  ptr0[j+1][k-1] +  ptr0[j-1][k-1];
                sumaPost=ptr0[j][k+16] +  ptr0[j+1][k+16] +  ptr0[j-1][k+16];
                v7=_mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,sumaAnte);
                v8=_mm_set_epi8(sumaPost,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
                //fin construccion marcos interiores
                v1=_mm_load_si128((__m128i*)&ptr0[j][k]);// fila del medio
                v3=_mm_load_si128((__m128i*)&ptr0[j-1][k]);//fila de arriba
                v4=_mm_load_si128((__m128i*)&ptr0[j+1][k]);//fila de abajo
                v2=_mm_adds_epi8(v1,v3);
                v2=_mm_adds_epi8(v2,v4);//suma del medio 'original'
                v5=_mm_srli_si128(v2,1);//desplaza izq
                v6=_mm_slli_si128(v2,1);//desplaza der
                v6=_mm_or_si128(v7,v6);//or entre la shift der y sumaAnte
                v5=_mm_or_si128(v5,v8);//or entre la shift izq y sumaPost
                v2=_mm_adds_epi8(v2,v5);//d+d'
                v2=_mm_adds_epi8(v2,v6);//d+d'+d''
                v2=_mm_sub_epi8(v2,v1);//resta la fila del medio (suma)
                //COMPARACION
                v2=_mm_or_si128(v2,v1);//OR entre la suma y el valor de la fila del medio, si debe vivir da como resultado 3 
                v9=_mm_set1_epi8(3);
                v2=_mm_cmpeq_epi8(v2,v9);//verifico si el resultado es igual a 3
                v9=_mm_set1_epi8(1);
                v2=_mm_and_si128(v2,v9);//AND entre el resultado de comparar y el registro con 1, me van a quedar con 1 solo los que tienen 0xFF
                //FIN COMPARACION
                _mm_store_si128((__m128i*)&ptrAux[j][k], v2);//almaceno la fila calculada en la matriz auxiliar
            }
    }
        //intercambio las matrices
   //     printf("ENTRA WAIT RANK:%d \n ",rank);
        MPI_Waitall(16,reqs2,MPI_STATUS_IGNORE);
        
//RECIVO BORDES
        //borde Inferior

        
        
        for(index=1;index<filas-1;index++){//ARMO UN BUFFER CON LAS FILA 0 Y LA ULTIMA FILA
            //ULTIMA ULTIMA COLUMNA
            total=bordeDer[index]+ptr0[index][columnas+14]+ptr0[index+1][columnas+15]+ptr0[index-1][columnas+15]+bordeDer[index-1]+ptr0[index-1][columnas+14]+bordeDer[index+1]+ptr0[index+1][columnas+14];
            if( ( ptr0[index][columnas+15]==1 && (total==3 || total==2) ) || ( ptr0[index][columnas+15]==0 && total==3) ){
                ptrAux[index][columnas+15]=1;
            }
            else{
                ptrAux[index][columnas+15]=0;
            }
            //FIN ULTIMA COLUMNA
            
            //PRIMERA COLUMNA
            total=bordeIzq[index]+ptr0[index][17]+ptr0[index+1][16]+ptr0[index-1][16]+bordeIzq[index-1]+ptr0[index-1][17]+bordeIzq[index+1]+ptr0[index+1][17];
            if( ( ptr0[index][16]==1 && (total==3 || total==2) ) || ( ptr0[index][16]==0 && total==3) ){
                ptrAux[index][16]=1;
            }
            else{
                ptrAux[index][16]=0;
            }
            //FIN ULTIMA COLUMNA
        }
        
        //Celda Primera Fila Primera Columna
        total=bordeIzq[0]+bordeSupIzq[0]+bordeIzq[1]+bordeSup[0]+ptr0[1][16]+bordeSup[1]+ptr0[0][17]+ptr0[1][17];
        if( (ptr0[0][16]==1 && (total==3 || total==2) ) || ( ptr0[0][16]==0 && total==3)){
            ptrAux[0][16]=1;
        }
        else{
            ptrAux[0][16]=0;
        }
        //Fin Celda Primera Fila Primera Columna
        
        //Celda Ultima Fila Primera Columna
        total=bordeIzq[filas-1]+bordeInfIzq[0]+bordeIzq[filas-2]+bordeInf[0]+ptr0[filas-1][17]+bordeInf[1]+ptr0[filas-2][17]+ptr0[filas-2][16];
        if( (ptr0[filas-1][16]==1 && (total==3 || total==2) ) || ( ptr0[filas-1][16]==0 && total==3)){
            ptrAux[filas-1][16]=1;
        }
        else{
            ptrAux[filas-1][16]=0;
        }
        //Fin Celda Ultima Fila Primera Columnaizq
        
        //Celda Primera Fila Ultima Columna
        total=bordeSup[columnas-1]+bordeSup[columnas-2]+ptr0[0][14+columnas]+ptr0[1][14+columnas]+ptr0[1][15+columnas]+bordeSupDer[0]+bordeDer[0]+bordeDer[1];
        if( (ptr0[0][15+columnas]==1 && (total==3 || total==2) ) || ( ptr0[0][15+columnas]==0 && total==3)){
            ptrAux[0][15+columnas]=1;
        }
        else{
            ptrAux[0][15+columnas]=0;
        }
        //Fin Celda Primera Fila Ultima Columna
        
        //Celda Ultima Fila Ultima Columna
        total=bordeInf[columnas-1]+bordeInf[columnas-2]+ptr0[filas-1][columnas-2]+ptr0[filas-2][columnas-2]+ptr0[filas-2][columnas-1]+bordeInfDer[0]+bordeDer[filas-1]+bordeDer[filas-2];
        if( (ptr0[filas-1][15+columnas]==1 && (total==3 || total==2) ) || ( ptr0[filas-1][15+columnas]==0 && total==3)){
            ptrAux[filas-1][15+columnas]=1;
        }
        else{
            ptrAux[filas-1][15+columnas]=0;
        }
        //Fin Celda Ultima Fila Ultima Columna
        
        //FILAS
        for(index=1;index<columnas-1;index++){
            //FILA ARRIBA
            total=ptr0[0][index+15]+ptr0[0][index+17]+ptr0[1][index+16]+bordeSup[index]+bordeSup[index-1]+bordeSup[index+1]+ptr0[1][index+15]+ptr0[1][index+17];
            if( ( ptr0[0][16+index]==1 && (total==3 || total==2) ) || ( ptr0[0][16+index]==0 && total==3) ){
                ptrAux[0][16+index]=1;
            }
            else{
                ptrAux[0][16+index]=0;
            }
            //FIN FILA ARRIBA
            
            //FILA ABAJO
            total=ptr0[filas-1][index+15]+ptr0[filas-1][index+17]+ptr0[filas-2][index+15]+ptr0[filas-2][index+16]+ptr0[filas-2][index+17]+bordeInf[index]+bordeInf[index-1]+bordeInf[index+1];
            if( ( ptr0[filas-1][16+index]==1 && (total==3 || total==2) ) || ( ptr0[filas-1][16+index]==0 && total==3) ){
                ptrAux[filas-1][16+index]=1;
            }
            else{
                ptrAux[filas-1][16+index]=0;
            }
            //FIN FILA ABAJO
        }
        //FIN FILAS
        
        

        aux = ptr0;
        ptr0 = ptrAux;	
        ptrAux = aux;
        }        
        
        
		//FIN PROCESAMIENTO INTERNO
        if(rank!=0){
            for(i=0;i<filas;i++){
                for(j=0; j<columnas; j++){
                    bufferfilas[j]= ptr0[i][j+16];// INGRESO DE VALORES A LA MATRIZ PROPIA DE CADA PROCESO}
                }
                MPI_Send(bufferfilas,columnas,MPI_CHAR,0,0,MPI_COMM_WORLD);
            }
        }
        if(rank==0){
            core=indiceCore=0;                           
            buffer1= malloc(sizeof(char)*columnas);
            //Proceso cero recibe las filas de otros procesos y las imprime
            for(index=0; index<rows; index++){
                for(i=0; i<scols; i++){//falta otro for por la cantidad de filas
                    if(indiceCore+core==0){//escribe proceso cero
                        for(k=0;k<columnas;k++){
                            if(ptr0[index][k+16]==1){
                                printf("0");//ingresa valores en la matriz de pc 0
                            }
                            else{
                                printf(".");
                            }
                        }    
                            core++;//aumenta el rank del core al que se le envia el msj
                    }
                    else{
                        if(core<size){
                            MPI_Recv(buffer1,columnas,MPI_CHAR,core,0,MPI_COMM_WORLD,status);//ENVIA FILAS A LOS PROCESOS
                            for(j=0;j<columnas;j++){
                                    if(buffer1[j]==1){
                                        printf("0");//ingresa valores en la matriz de pc 0
                                    }
                                    else{
                                        printf(".");
                                    }
                            }
                            core++;//aumenta el rank del core al que se le envia el msj
                        }
                    }
                    if(core==indiceCore+scols)
                        core=indiceCore;//incrementa o reset el valor del core
                }
                if((index+1) % filas == 0){
                    indiceCore+= scols;
                    core = indiceCore;
                }
                printf("\n");
            }
        }
    MPI_Finalize();
    return 0;
}
