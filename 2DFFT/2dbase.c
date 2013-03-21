/*
 * n x n complex 2 dimensional fast Fourier transform (2DFFT)
 * David O'Hallaron, Carnegie Mellon University
 */

#include"man.h"
#include"image_utils/img_write.h"

#include "fft.h"

#include <stdio.h>
#include <math.h>
#include <sys/time.h>

//static mycomplex a[MAXN][MAXN];     /* input matrix */
mycomplex **a;
//static mycomplex w_common[MAXN/2];  /* twiddle factors */
mycomplex *w_common;

main(argc,argv)
int argc;
char **argv;
{
  int n;           /* FFT size */
  int nin;
  int r,c;
  int k,l;
  int logn;        /* log base 2 of n */
  int errors,sign; /* used for error checking */
  int nx;          /* used to compute logn */
  int flops;       /* total number of floating point ops */
  float mflops;    /* Mflops/s */
  double fsecs;    /* time spent doing 1D FFT's */
  double tsecs;    /* time spent doing the transpose */
  double secs;     /* total time of the 2D FFT */
  int i,j;         /* index variables */

  struct timeval start, finish; /* gettimeofday structures */
  void print_cmat(), print_cvec(); /* forward definitions */

  /* get the FFT size and ensure that it's in range */
  if (argc < 2) {
    (void)fprintf(stderr, "usage: %s <n>\n", argv[0]);
    exit(0);
  }
  nin = r = c = atoi(argv[1]);

  if (argc>2) {
	c = atoi(argv[2]);
	if (c>r) nin=c; 
  }
  
  if (nin < 2) {
    (void)fprintf(stderr, "%s: bad fft size\n", argv[0]);
    exit(0);
  }

// AJ Hard Code
  c=1024;
  r=1024;
  nin=1024;
//

  /* compute logn and ensure that n is a power of two */
  nx = nin;
  logn = 0;
  while(( nx >>= 1) > 0) 
    logn++; 
  nx = 1;
  for (i=0; i<logn; i++)
    nx = nx*2;

//AJ -s
  n = nx<nin?nx*2:nin;
  printf("n=%d, nx=%d, r=%d, c=%d\n", nin, n, r, c);
//AJ -e

//Allocate memory for the Matrix and twiddle factor
//AJ -s
  a = (mycomplex **)malloc(n*sizeof(mycomplex *));
  for(i=0;i<n;i++){
    a[i]=(mycomplex *)malloc(n*sizeof(mycomplex));
  }
  w_common=(mycomplex *)malloc((n/2)*sizeof(mycomplex));
//AJ -e

//AJ -s 
  /* initialize the input matrix with a centered point source */
 /*
 for (i=0; i<r; i++)
    for (j=0; j<c; j++) 
      a[i][j].r = a[i][j].i = 0.0;

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Padding all the extra columns already accessed
//////////////////////////////////////////////////////////////////////////////////////////////////////

  for (k=0; k<i; k++)
    for (l=j; l<n; l++)
	a[k][l].r = a[k][l].i = 0.0;
	
// Print the last column accessed
  printf("i=%d, j=%d\n",i,j);
 
// Padding the last rows 
  for (i; i<n; i++)
    for (j; j<n; j++)
      a[i][j].r = a[i][j].i = 0.0;

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Done Padding
//////////////////////////////////////////////////////////////////////////////////////////////////////

a[n/2][n/2].r =  a[n/2][n/2].i = (float)n;
*/

  for(i=0; i<n; i++)
    for (j=0; j<n; j++){
	a[i][j].r=in[i*n+j];
	a[i][j].i = 0.0;
    }

  /* precompute the complex constants (twiddle factors) for the 1D FFTs */
  for (i=0;i<n/2;i++) {
    w_common[i].r = (float) cos((double)((2.0*PI*i)/(float)n));
    w_common[i].i = (float) -sin((double)((2.0*PI*i)/(float)n));
  }
  
  /* 
   * now do the actual 2D FFT
   */

  /* first do a set of row FFTs */
  (void)gettimeofday(&start, 0);
  for (i=0; i<n; i++)
    (void)fft(&a[i][0], w_common, n, logn);
  (void)gettimeofday(&finish, 0);
  fsecs = ((((finish.tv_sec - start.tv_sec) * 1000000.0) +
	    finish.tv_usec) - start.tv_usec) / 1000000.0;

  /* then transpose the matrix */
  (void)gettimeofday(&start, 0);
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++) 
      SWAP(a[i][j], a[j][i]); 
  (void)gettimeofday(&finish, 0);
  tsecs = ((((finish.tv_sec - start.tv_sec) * 1000000.0) +
	    finish.tv_usec) - start.tv_usec) / 1000000.0;
      
  /* then do another set of row FFTs */
  (void)gettimeofday(&start, 0);
  for (i=0; i<n; i++)
    (void)fft(&a[i][0], w_common, n, logn);
  (void)gettimeofday(&finish, 0);
  fsecs += ((((finish.tv_sec - start.tv_sec) * 1000000.0) +
	    finish.tv_usec) - start.tv_usec) / 1000000.0;

  /* transpose the matrix back */
  (void)gettimeofday(&start, 0);
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++) 
      SWAP(a[i][j], a[j][i]); 
  (void)gettimeofday(&finish, 0);
  tsecs += ((((finish.tv_sec - start.tv_sec) * 1000000.0) +
	    finish.tv_usec) - start.tv_usec) / 1000000.0;
      
  /* check the answers for an alternating sequence of (+-n,+-n) */
  errors = 0;
  for (i=0; i<n; i++) {
    if (((i+1)/2)*2 == i) 
      sign = 1;
    else
      sign = -1;
    for (j=0; j<n; j++) {
      if (a[i][j].r > n*sign+EPSILON ||
	  a[i][j].r < n*sign-EPSILON ||
	  a[i][j].i > n*sign+EPSILON ||
	  a[i][j].i < n*sign-EPSILON) {
//	printf("%d,%d, %f + j %f\n",i,j,sign*a[i][j].r, sign*a[i][j].i);
	errors++;
      }
      sign *= -1;
    }
  }
  if (errors) { 
    printf("%d errors!\n", errors);
    //exit(0);
  }

  /* summarize the 2d FFT performance */
  printf("%d x %d 2D FFT\n", n, n);
  secs = fsecs + tsecs;
  flops = (n*n*logn)*10;
  mflops = ((float)flops/1000000.0);
  mflops = mflops/(float)secs;
  printf("\n");
  printf("1D FFTs   : %10.6f secs (%2d%%)\n", fsecs, (int)(fsecs/secs*100));
  printf("transpose : %10.6f secs (%2d%%)\n", tsecs, (int)(tsecs/secs*100));
  printf("total     : %10.6f secs\n", secs);
  printf("\n");
  printf("%10.6f Mflop/s\n", mflops);

//////////////////////////////////////
/// INV FFT
/////////////////////////////////////
  /* precompute the complex constants (twiddle factors) for the 1D FFTs */
  for (i=0;i<n/2;i++) {
    w_common[i].r = (float) cos((double)((2.0*PI*i)/(float)n));
    w_common[i].i = (float) sin((double)((2.0*PI*i)/(float)n));
  }

  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++) 
      SWAP(a[i][j], a[j][i]);

  for (i=0; i<n; i++)
    (void)fft(&a[i][0], w_common, n, logn); 

  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++) 
      SWAP(a[i][j], a[j][i]);

  for (i=0; i<n; i++)
    (void)fft(&a[i][0], w_common, n, logn);



////////////////////////////////////

// AJ -s
  int *out = (int *)malloc(sizeof(int)*r*c);
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
	out[i*n+j]=(unsigned int)(a[i][j].r/(n*n));

  if(memcmp(in, out, r*c*sizeof(int))!=0) printf("ERROR: Input and Output not Same\n");;
  img_write(out, 1024, 1024);
  for(i=0;i<n;i++)
  	free(a[i]);
  free(a);
  free(w_common);
  free(out);
// AJ -e
  exit(0);
}

/*
 * routines to print complex matrices and vectors
 */
void print_cvec(a, n)
mycomplex *a;
int n;
{
  int i;

  for (i=0; i<n; i++) {
    if (i%4 == 0 && i)  
      printf("\n"); 
    printf("(%.3f,%.3f) ", a[i].r, a[i].i);   
  }
  printf("\n");
}

void print_cmat(a, m, n)
mycomplex a[][MAXN];
int m,n;
{
  int i,j;

  for (i=0; i<m; i++) {
    for (j=0; j<n; j++)
      printf("(%6.2f,%6.2f) ",a[i][j].r, a[i][j].i);
    printf("\n");
  }
}





