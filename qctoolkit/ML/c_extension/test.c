#include "kernels_capi.h"
#include <stdio.h>
#include <stdlib.h>

int main(){
  double *vec1 = (double*) malloc(5*sizeof(double));
  double *vec2 = (double*) malloc(5*sizeof(double));
  double *inp = (double*) malloc(sizeof(double));
  double out;
  int i;

  inp[0] = 3.0;

  for(i=0;i<5;i++){
    vec1[i] = i;
    vec2[i] = 2*i;
  }
  kernelpt kernel = kernel_create("Gaussian", inp);
  out = kernel_evaluate("Gaussian",kernel,vec1,vec2,5);
  printf("%lf\n",out);
  //out = g.evaluate(vec1,vec2,5);

  kernel_free("Gaussian", kernel);
}
