/*
 * fft -  n-point in-place decimation-in-time FFT of complex vector 
 *        "data" using the n/2 complex twiddle factors in "w_common".
 */
#include "fft.h"

void fft(data,w_common,n,logn)
mycomplex *data, *w_common;
int n,logn;
{
  int incrvec, i0, i1, i2, nx;
  float f0, f1;

  void bit_reverse();

  /* bit-reverse the input vector */
  (void)bit_reverse(data,n);

  /* do the first logn-1 stages of the fft */
  i2 = logn;
  for (incrvec=2;incrvec<n;incrvec<<=1) {
    i2--;
    for (i0 = 0; i0 < incrvec >> 1; i0++) {
      for (i1 = 0; i1 < n; i1 += incrvec) {
        f0 = data[i0+i1 + incrvec/2].r * w_common[i0<<i2].r - 
	  data[i0+i1 + incrvec/2].i * w_common[i0<<i2].i;
        f1 = data[i0+i1 + incrvec/2].r * w_common[i0<<i2].i + 
	  data[i0+i1 + incrvec/2].i * w_common[i0<<i2].r;
        data[i0+i1 + incrvec/2].r = data[i0+i1].r - f0;
        data[i0+i1 + incrvec/2].i = data[i0+i1].i - f1;
        data[i0+i1].r = data[i0+i1].r + f0;
        data[i0+i1].i = data[i0+i1].i + f1;
      }
    }
  }

  /* do the last stage of the fft */
  for (i0 = 0; i0 < n/2; i0++) {
    f0 = data[i0 + n/2].r * w_common[i0].r - 
      data[i0 + n/2].i * w_common[i0].i;
    f1 = data[i0 + n/2].r * w_common[i0].i + 
      data[i0 + n/2].i * w_common[i0].r;
    data[i0 + n/2].r = data[i0].r - f0;
    data[i0 + n/2].i = data[i0].i - f1;
    data[i0].r = data[i0].r + f0;
    data[i0].i = data[i0].i + f1;
  }
}

/* 
 * bit_reverse - simple (but somewhat inefficient) bit reverse 
 */
void bit_reverse(a,n)
mycomplex *a;
int n;
{
  int i,j,k;

  j = 0;
  for (i=0; i<n-2; i++){
    if (i < j)
      SWAP(a[j],a[i]);
    k = n>>1;
    while (k <= j) {
      j -= k; 
      k >>= 1;
    }
    j += k;
  }
}
