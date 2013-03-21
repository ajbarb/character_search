/*
 * n x n complex 2 dimensional fast Fourier transform (2DFFT)
 * David O'Hallaron, Carnegie Mellon University
 */

#include "fft.h"

#include <stdio.h>
#include <math.h>
#include <sys/time.h>

//static mycomplex a[MAXN][MAXN];     /* input matrix */
//static mycomplex w_common[MAXN/2];  /* twiddle factors */
mycomplex *w_common;

void dfft(mycomplex **a, int n, int inv)
{
  int logn;        /* log base 2 of n */
  int nx;          /* used to compute logn */
  int i,j;         /* index variables */

  if (inv == 0 || inv >1 || inv <-1) {
    printf("ERROR: Provide the correct value for inv \n");
    abort();
  }

  nx=n;
  logn = 0;
  while(( nx >>= 1) > 0)
    logn++;

  /* compute logn and ensure that n is a power of two */
  
  w_common=(mycomplex *)malloc((n/2)*sizeof(mycomplex));

  /* precompute the complex constants (twiddle factors) for the 1D FFTs */
  for (i=0;i<n/2;i++) {
    w_common[i].r = (float) cos((double)(inv*(2.0*PI*i)/(float)n));
    w_common[i].i = (float) sin((double)(inv*(2.0*PI*i)/(float)n));
  }
 
  /* 
   * now do the actual 2D FFT
   */

  /* first do a set of row FFTs */
  for (i=0; i<n; i++)
    (void)fft(&a[i][0], w_common, n, logn);

  /* then transpose the matrix */
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++) 
      SWAP(a[i][j], a[j][i]); 
      
  /* then do another set of row FFTs */
  for (i=0; i<n; i++)
    (void)fft(&a[i][0], w_common, n, logn);

  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++) 
      SWAP(a[i][j], a[j][i]); 
  
  if(inv==1){
    for (i=0; i<n; i++){
	for (j=0; j<n; j++){
	   a[i][j].r=(a[i][j].r+1)/(n*n); 
	}
    }
  }
}
