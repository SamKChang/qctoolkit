#include<Python.h>
#include <numpy/arrayobject.h>

/////////////////////////////////////////
// main function for reading CUBE file //
/////////////////////////////////////////
// or what ever computation need to be done
// only minor changes should be necessary
// to interface other type of calculation
void readcube_c(char *inp, double **out, int dims[3]){
  FILE *fp;
  char *string = (char *) malloc(80);
  size_t len=0;
  ssize_t read;
  int i;

  //double *data = (double *) malloc(40 * sizeof(double));
  *out = (double *) malloc(40 * sizeof(double));
  fp = fopen(inp, "r");
  read = getline(&string, &len, fp);
  printf("%s", string);

  dims[0] = 2;
  dims[1] = 2;
  dims[2] = 2;
  
  for(i=0;i<8;i++){
    (*out)[i] = (double) i;
    printf("%d, %lf\n", i, (*out)[i]);
  }

  //*out = (char *) malloc(100);
  fclose(fp);
  free(string);
}

///////////////////////////////////////
// !!! PYTHON INTERFACE FROM HERE!!! //
///////////////////////////////////////

//////////////////////////////
// python callable function //
//////////////////////////////
// input: string, for file name
// output: PyArrayObject, as volumetric data
static PyObject * read_cube(PyObject * self, PyObject * args){
  char *input; // filename string as input
  double *result; // C-data need to be converted to numpy object
  int dims[3];
  PyObject *output; // numpy object pass back to python

  // parse arguments check and/or error handling
  if(!PyArg_ParseTuple(args, "s", &input))
    return NULL;

  // run the actual function
  // NOTE the result is passed as datatype double**
  // to allocate memory in the function
  readcube_c(input, &result, dims);

  // numpy C-API: allocate memory to copy C-array data
  // Assume 1-D C-array, N-D numpy array will be filled as
  // [result[0]         : result[dims[0]], 
  //  result[dims[0]+1] : result[dims[1]],
  //  result[dims[1]+1] : result[dims[2]],
  //                 .......             ,
  // ]
  output = PyArray_SimpleNewFromData(3, dims, NPY_DOUBLE, result);

  return output;
  // NOTE result CANNOT be freed before returning output to python
  // otherwise the first element of output will be overwritten!
  free(result);
}

/////////////////////////////////////////
// register python module symbol table //
/////////////////////////////////////////
// PyMethodDef: struct of four field
//   string: method_name
//   PyCFunction: method_function
//   int: flag
//   string: documentation
static PyMethodDef ReadCubeMethods[] = {
  {"read_cube", read_cube, METH_VARARGS, "read CUBE file"},
  {NULL, NULL, 0, NULL} // sentinel?
};

//////////////////////////
//  method constructor  //
//////////////////////////
// the name MUST be init{name} otherwise python cannot find it
// depends on numpy C-API, import_array()
// and/or import_ufunc() are necessary
// otherwise the code return segfault
PyMODINIT_FUNC initread_cube(void){
  Py_InitModule("read_cube", ReadCubeMethods);
  import_array(); // necessary!
}
