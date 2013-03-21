#include<stdio.h>
#include"template_source.h"
#include"fft.h"
#include"source_image.h"

mycomplex multiplyMatrix(mycomplex x,mycomplex y);

void main(){
  mycomplex **template, **in, **out;
  int n, in_r, in_c, template_r, template_c;
  int nx;
  int logn;
  int i,j;
  int threshold=256;

  in_r=512;
  in_c=1024;
  template_r=26;
  template_c=46;

  nx = in_r>in_c?in_r:in_c;
  n = nx;
  logn = 0;
  while(( nx >>= 1) > 0)
    logn++;
  nx = 1;
  for (i=0; i<logn; i++)
    nx = nx*2;

  n = nx<n?nx*2:n;

  template = (mycomplex **)malloc(n*sizeof(mycomplex *));
  for(i=0;i<n;i++){
    template[i]=(mycomplex *)malloc(n*sizeof(mycomplex));
  }

  for (i=0; i<template_r; i++)
    for (j=0; j<template_c; j++){
      template[i][j].r=temp_s[i*template_c+j];
      template[i][j].i = 0.0;
    }

  in = (mycomplex **)malloc(n*sizeof(mycomplex *));
  for(i=0;i<n;i++){
    in[i]=(mycomplex *)malloc(n*sizeof(mycomplex));
  }

  for (i=0; i<in_r; i++)
    for (j=0; j<in_c; j++){
      in[i][j].r = in_s[i*in_c+j];
      in[i][j].i = 0.0;
  }

  rotate180(template, template_r, template_c);
 
  out = (mycomplex **)malloc(n*sizeof(mycomplex *));
  for(i=0;i<n;i++){
    out[i]=(mycomplex *)malloc(n*sizeof(mycomplex));
  }


for (i=0; i<template_r; i++)
    for (j=template_c; j<n; j++)
        template[i][j].r = template[i][j].i = 0.0;

for (i=template_r; i<n; i++)
    for (j=template_c; j<n; j++)
      template[i][j].r = template[i][j].i = 0.0;


for (i=0; i<in_r; i++)
    for (j=in_c; j<n; j++)
      in[i][j].r = in[i][j].i = 0.0;

for (i=in_r; i<n; i++)
    for (j=in_c; j<n; j++)
      in[i][j].r = in[i][j].i = 0.0;


dfft(template, n, -1);
dfft(in, n, -1);

for(i=0;i<n;i++)
  for(j=0;j<n;j++){
    out[i][j]=multiplyMatrix(in[i][j], template[i][j]);    
  }

  dfft(out, n, 1);

for(i=0;i<n;i++)
  for(j=0;j<n;j++){
    if(out[i][j].r<threshold) out[i][j].r=0.0;
  }


  img_write(out, in_c, in_r);

}

mycomplex multiplyMatrix(mycomplex x,mycomplex y){
    mycomplex z;
    z.r = x.r*y.r - x.i*y.i;
    z.i = x.r*y.i + x.i*y.r;
    return z;
}
