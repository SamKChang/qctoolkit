#include<math.h>
#include<string.h>
#include"kernels.h"

/**********************
** interface wrapper **
**********************/
double kerEval(char *type, double *input,
               double *vec1, double *vec2, int size) {
  if(strcmp(type,"Gaussian")==0)
    return GaussianKernel(input, vec1, vec2, size);
}

/**********************
**  Gaussian kernel  **
**********************/
double GaussianKernel(double *input, 
                      double *vec1, double *vec2, int size){
  /* the rest of input array has no effect */
  double sigma = input[0];
  double norm = 0;
  double diff;
  int i;
  for(i=0;i<size;i++){
    diff = vec1[i] - vec2[i];
    norm += diff*diff;
  }
  norm = sqrt(norm);
  return exp(-0.5 * pow((norm / sigma), 2));
}

