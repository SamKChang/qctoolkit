// C routine to read CUBE file

#include <Python.h>
#include <numpy/arrayobject.h>
#include "utilities.h"

///////////////////////////////////////////////////
// main function for constructing Coulomb matrix //
///////////////////////////////////////////////////
// or what ever computation need to be done
// only minor changes should be necessary
// to interface other type of calculation
void coulombmatrix(char* inp, 
                    int dim, 
                    double** coulombMatrix,
                    int dims[2])
{
  FILE *fr;
  int Na;
  double read;
  double *coord;
  double Ri[3], Rj[3], dR;
  int *Z;
  int i=0, j=0, k=0, s1, s2;
  char *string = (char *) malloc(80);
  size_t len=0;


  // read xyz file into coord and Z
  fr = fopen(inp, "r");
  fscanf(fr, "%d", &Na);
  coord = (double *) malloc(Na * 3 * sizeof(double));
  Z = (int *) malloc(Na * sizeof(int));
  for(i=0;i<2;i++){
    getline(&string, &len, fr);
  }
  for(i=0;i<Na;i++){
    fscanf(fr, "%s", string);
    Z[i] = n2Z(string);
    for(j=0;j<3;j++) fscanf(fr, "%lf", &coord[3*i+j]);
  }

  // construct Coulomb matrix
  *coulombMatrix = (double *) malloc(dim * dim * sizeof(double));
  // initialize Coulomb matrix
  for(i=0;i<dim*dim;i++) (*coulombMatrix)[i] = 0;
  // fill in matrix elements
  for(i=0;i<Na;i++){
    for(k=0;k<3;k++) Ri[k] = coord[i*3 + k];
    for(j=0;j<Na;j++){
      s1 = dim*i + j;
      if(i==j)(*coulombMatrix)[s1] = 0.5 * pow(Z[i], 2.4);
      else{
        s2 = dim*j + i;
        dR = 0;
        for(k=0;k<3;k++) dR += sqr(Rj[k] - coord[j*3 + k]);
        (*coulombMatrix)[s1] = Z[i]*Z[j]/dR;
        (*coulombMatrix)[s2] = Z[i]*Z[j]/dR;
      }
    }
  }

  fclose(fr);
  free(string);
  free(coord);
  free(Z);
}

///////////////////////////////////////
// !!! PYTHON INTERFACE FROM HERE!!! //
///////////////////////////////////////

//////////////////////////////
// python callable function //
//////////////////////////////
// input: string, for file name
// output: PyArrayObject, as volumetric data
static PyObject * coulomb_matrix(PyObject * self, PyObject * args){
  char *input; // filename string as input
  int dim; // base dimension of coulomb matrix
  // C-data need to be converted to numpy object
  double *coulombMatrix; 
  // dimensions of numpy array output 
  int dims[2];
  // numpy object pass back to python
  PyObject *npCoulombMatrix; 

  // parse arguments check and/or error handling
  if(!PyArg_ParseTuple(args, "sn", &input, &dim))
    return NULL;
  dims[0] = dim;
  dims[1] = dim;

  // run the actual function
  // NOTE the data is passed as datatype double**
  // to allocate memory in the function
  coulombmatrix(input, dim, &coulombMatrix, dims);

  // numpy C-API: allocate memory to copy C-array data
  // Assume 1-D C-array, N-D numpy array will be filled as
  // [data[0]         : data[dims[0]], 
  //  data[dims[0]+1] : data[dims[1]],
  //  data[dims[1]+1] : data[dims[2]],
  //                 .......             ,
  // ]
  npCoulombMatrix = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, 
                                              coulombMatrix);

  return Py_BuildValue("O", npCoulombMatrix);
  //return PyArray_Return(npCoulombMatrix);
  //return npCoulombMatrix;
  printf("yo\n");

  // NOTE data CANNOT be freed before returning output to python
  // otherwise the first element of output will be overwritten!
  free(coulombMatrix);
}

/////////////////////////////////////////
// register python module symbol table //
/////////////////////////////////////////
// PyMethodDef: struct of four field
//   string: method_name
//   PyCFunction: method_function
//   int: flag
//   string: documentation
static PyMethodDef CoulombMatrixXYZ[] = {
  {"coulomb_matrix", 
   coulomb_matrix, 
   METH_VARARGS, 
   "read xyz file and construct Coulomb matrix"
  },
  {NULL, NULL, 0, NULL} // sentinel?
};

//////////////////////////
//  method constructor  //
//////////////////////////
// the name MUST be init{name} otherwise python cannot find it
// depends on numpy C-API, import_array()
// and/or import_ufunc() are necessary
// otherwise the code return segfault
PyMODINIT_FUNC initcoulomb_matrix(void){
  Py_InitModule("coulomb_matrix", CoulombMatrixXYZ);
  import_array(); // necessary!
}
