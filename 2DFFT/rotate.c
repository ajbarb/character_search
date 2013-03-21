#include<stdio.h>
#include"fft.h"

void rotate180_square(int **a, int n) {

  int i,j,temp;
  for(i=0; i<=n/2; i++){
    for(j=0; j<n&&(i*n+j)<(int)(n*n/2); j++){
      temp = a[i][j];
      a[i][j] = a[n-i-1][n-j-1];
      a[n-i-1][n-j-1] = temp;
    }
  }

}

void rotate180(mycomplex **a, int r, int c) {
  int i,j,temp;
  for(i=0; i<=r/2; i++){
    for(j=0; j<c&&(i*c+j)<(int)(r*c/2); j++){
      temp = a[i][j].r;
      a[i][j].r = a[r-i-1][c-j-1].r;
      a[r-i-1][c-j-1].r = temp;
    }
  }
}
