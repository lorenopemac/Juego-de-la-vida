 
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
const int RMAX=100;

void main (int argc, char* argv[])
{
    int phase,i,n,*a,tmp;
    n=5000;
    a=malloc(n*sizeof(int));
    srandom(1);
    for(i=0;i<n;i++)
        a[i]=random() % RMAX;

    for(phase=0;phase<n;phase++)
        if(phase % 2 == 0)
            for(i = 1;i<n; i+=2)
            {
                if(a[i-1]>a[i])
                {
                    tmp= a[i - 1];
                    a[i-1]=a[i];
                    a[i]=tmp;
                }
            }
        else
            for(i = 1;i < n-1; i +=2)
            {
                if(a[i]>a[i+1])
                {
                    tmp=a[i+1];
                    a[i+1]=a[i];
                    a[i]=tmp;
                }
            }
}