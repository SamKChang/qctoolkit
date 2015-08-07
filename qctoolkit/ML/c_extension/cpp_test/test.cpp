#include "kernels.h"
#include <iostream>
#include <stdio.h>

int main(){
  double *vec1 = new double[5];
  double *vec2 = new double[5];
  double out;
  for(int i=0;i<5;i++){
    vec1[i] = i;
    vec2[i] = 2*i;
  }
  Gaussian g(3.0);
  out = g.evaluate(vec1,vec2,5);
  std::cout << out << std::endl;
}
