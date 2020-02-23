// Peña Lorenzo
#include <stdio.h>
#include <stdlib.h>
#include <emmintrin.h>
int main()
{
char *grilla,*grilla2,*grilla3,*buffer,*resultado,*temp1,*temp2;
int cols,rows,juegos,k,j,l,n,r,s,p,t,q,w,z,i,res,bordeFin,sumaAnte,sumaPost;
__m128i v1,v2,v3,v4,v5,v6,v7,v8,v9;
// Lectura del Archivo
FILE *fp;
fp=fopen("salida.cells","r");
fscanf(fp,"cols %i\n rows %i\n steps %i",&cols,&rows,&juegos);
cols=cols+17;
rows=rows+2;
res= cols%16;//Buscar mod
bordeFin=cols;
if(res!=0)
{
 cols=cols+(16-res);
}
if(posix_memalign((void**)&grilla, 16, cols*rows*sizeof(char)) != 0)
      return 1;
if(posix_memalign((void**)&grilla2, 16, cols*rows*sizeof(char)) != 0)
      return 1;
buffer=malloc((cols+1)*sizeof(char));
resultado=malloc(cols*sizeof(char));
// Ingreso de valores a la matriz
  resultado=fgets(buffer,cols,fp);
  memset(buffer,'\0',cols);
  w=0;
  while(resultado!=NULL)
  {
   for(q=16;q<cols+1;q++)
   {	  
	  if(buffer[q-16]=='O')
	  {
	 	grilla[w*cols+q]=1;	
	  }
	  else
	  {
		grilla[w*cols+q]=0;
	  }
   }
   if(!feof(fp))
   {
       memset(buffer,'\0',cols);
       resultado=fgets(buffer,cols,fp);
   }
   w=w+1;
}
fclose(fp);
for(z=0;z<juegos;z++)
{
/*Creacion de marco de la matriz*/
  grilla[0*cols+15]=grilla[(rows-2)*cols+(bordeFin-2)];// bordeFin son la cantidad de columnas
  grilla[(rows-1)*cols+(bordeFin-1)]=grilla[1*cols+16];
  grilla[(rows-1)*cols+15]=grilla[1*cols+(bordeFin-2)];
  grilla[0*cols+(bordeFin-1)]=grilla[(rows-2)*cols+16];
  for(l=1;l<(rows-1);l++) 
   {
	grilla[l*cols+15]=grilla[l*cols+(bordeFin-2)];
	grilla[l*cols+(bordeFin-1)]=grilla[l*cols+16];   
   }
  for(n=16;n<(bordeFin-1);n++)
   { 
	grilla[0*cols+n]=grilla[(rows-2)*cols+n];
	grilla[(rows-1)*cols+n]=grilla[1*cols+n];
   }
//Fin creacion de marco
for(r=1;r<rows-1;r++)//inicio J vida(rows-1 por las 2 filas del marco)
 {
  for(i=1;i<(cols/16);i++)//INICIO RECORRIDO COLUMNAS 
  {
//inicio construccion marcos interiores
    sumaAnte=grilla[r*cols+i*16-1] + grilla[(r-1)*cols+i*16-1] + grilla[(r+1)*cols+i*16-1];
    sumaPost=grilla[r*cols+i*16+16] + grilla[(r-1)*cols+i*16+16] + grilla[(r+1)*cols+i*16+16];
    v7=_mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,sumaAnte);
    v8=_mm_set_epi8(sumaPost,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
 //fin construccion marcos interiores
    v1=_mm_load_si128(((__m128i*)&grilla[r*cols])+i);// fila del medio
    v3=_mm_load_si128(((__m128i*)&grilla[(r-1)*cols])+i);//fila de arriba
    v4=_mm_load_si128(((__m128i*)&grilla[(r+1)*cols])+i);//fila de abajo
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
    _mm_store_si128(((__m128i*)&grilla2[r*cols])+i, v2);//almaceno la fila calculada en la matriz auxiliar
  }			//FIN RECORRIDO COLUMNAS
 }			//FIN RECORRIDO FILAS
/*Intercambio punteros*/
  grilla3=grilla;
  grilla=grilla2;
  grilla2=grilla3;
/*Fin intercambio punteros*/
}
//imprimir Matriz
  res=0;
  for(s=1;s<(rows-1);s++)
  {
 	for(r=16;r<bordeFin-1;r++)
	{
            if(grilla[s*cols+r])
            {
                res++;    
                printf("%c",'O');
            }
           else
           {
                printf("%c",'.');
           }
        }
	printf("%s\n","");
  }
  printf("Number of live cells = %d",res);
//fin imprimir
}
