#include <iostream>
#include "descriptors.h"

CoulombMatrix::CoulombMatrix(int size){
  base_dim = size;
  matrix = new double[size*size];
}

double* CoulombMatrix::getVector(){
  return matrix;
}
