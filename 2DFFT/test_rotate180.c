#include<stdio.h>
#include"fft.h"

void main(){
  int n=5;
  int i,j;
  int r=3,c=2;
  int **a;
  a=(int**)malloc(sizeof(int*)*n);
  for(i=0;i<n;i++)
    a[i]=(int*)malloc(sizeof(int)*n);

  for(i=0;i<r;i++){
    for(j=0; j<c;j++){
       a[i][j]=i*n+j;
	printf("%d ",a[i][j]);
    }
    printf("\n");
  }

  printf("\n");

  rotate180(a, r, c);

  for(i=0;i<r;i++){
    for(j=0; j<c;j++){
        printf("%d ",a[i][j]);
    }
    printf("\n");
  }

}
