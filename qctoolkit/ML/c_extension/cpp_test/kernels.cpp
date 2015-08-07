#include<iostream>
#include<math.h>
#include"kernels.h"
#include<string.h>

/**********************
**   C wrapper API   **
***********************/
// allocate memory
extern "C" void* kernel_create(char *type, double *input){
  if(strcmp(type,"Gaussian")==0)
    return new Gaussian(input[0]);
}
// free memory
extern "C" void kernel_free(char *type, void *kernel) {
  if(strcmp(type,"Gaussian")==0)
    delete static_cast<Gaussian*>(kernel);
}
extern "C" double kernel_evaluate(char *type,
                                void *kernel, 
                                double *vec1, 
                                double *vec2, 
                                int size) {
  if(strcmp(type,"Gaussian")==0)
    static_cast<Gaussian*>(kernel)->evaluate(vec1,vec2,size);
}

/**********************
**  Gaussian kernel  **
**********************/
// constructor, take one input
Gaussian::Gaussian(double input){
  sigma = input;
}
// interface, evaluate distance
double Gaussian::evaluate(double *vec1, double *vec2, int size){
  double norm = 0;
  double diff;
  for(int i=0;i<size;i++){
    diff = vec1[i] - vec2[i];
    norm += diff*diff;
  }
  norm = sqrt(norm);
  return exp(-0.5 * pow((norm / sigma), 2));
}

