// C routine to read CUBE file

#include<Python.h>
#include <numpy/arrayobject.h>

///////////////////////////////////////////////////
// read CUBE file to determin needed memory size //
///////////////////////////////////////////////////
// open CUBE file but read ONLY the header part
// for main python code to allocate Numpy array
void readCubeHeader(char *inp, int dims[3], int size[2]){
  FILE *fr;
  double read;
  int i, j=0;
  double grid[16];
  char *string = (char *) malloc(80);
  size_t len=0;

  // CUBE file contains fixed size grid specification
  fr = fopen(inp, "r");
  for(i=0;i<2;i++) getline(&string, &len, fr);
  for(i=0;i<16;i++){
    fscanf(fr, "%lf", &read);
    grid[i] = read;
  }
  getline(&string, &len, fr);

  size[0] = grid[0];
  size[1] = 4;
  for(i=0;i<3;i++){
    dims[i] = (int) grid[4 * (i+1)];
  }
  fclose(fr);
  free(string);
}

/////////////////////////////////////////
// main function for reading CUBE file //
/////////////////////////////////////////
// or what ever computation need to be done
// only minor changes should be necessary
// to interface other type of calculation
void readcube_c(char *inp, 
                double *cube, 
                double *structure, 
                double *grid, 
                int dims[3],
                int size[2])
{
  FILE *fr;
  int Na, N3=1;
  double read;
  int i, j=0;
  char *string = (char *) malloc(80);
  size_t len=0;

  // CUBE file contains fixed size grid specification
  fr = fopen(inp, "r");
  for(i=0;i<2;i++) getline(&string, &len, fr);
  for(i=0;i<16;i++){
    fscanf(fr, "%lf", &read);
    grid[i] = read;
  }
  getline(&string, &len, fr);

  // assign proper value to related parameters
  Na = size[0];
  for(i=0;i<3;i++) N3 *= dims[i];

  for(i=0;i<Na*5;i++){
    fscanf(fr, "%lf", &read);
    if((i%5) != 0){
      structure[j++] = read;
    }
  }

  // check extra line for molecular orbital index
  // signiture of corner wavefunction/density 
  // value is always a small
  fscanf(fr, "%lf", &read);
  if(read*read<1.0E-4){
    cube[0] = read; 
    j=1; // without MO index
  }else{
    getline(&string, &len, fr);
    printf("skipping line: (double)% le (string)%s\n", 
           read, string);
    j=0; // with MO index
  }

  for(i=j;i<N3;i++){
    if(!feof(fr)){
      fscanf(fr,"%le", &read);
      cube[i] = read;
    }else{
      printf("===== ERROR REPORT =====");
      printf("last point read: % le\n", cube[i-1]);
      printf("end at i=%d, while N=%d\n", i, N3);
      // check for incomplete CUBE file
      printf("filename: %s\n", inp);
      printf("====== END REPORT ======");
      break;
    }
  }
  fclose(fr);
  free(string);
}

//////////////////////////////////////////////
// redirect PyArray data contant to C-array //
//////////////////////////////////////////////
double *pyvector_to_Carrayptrs(PyArrayObject *arrayin){
  int n=arrayin->dimensions[0];
  /* pointer to arrayin data as double */
  return (double *) arrayin->data;
}

// !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*! //
///////////////////////////////////////
// !!! PYTHON INTERFACE FROM HERE!!! //
///////////////////////////////////////
// !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*! //

//////////////////////////////
// python callable function //
//////////////////////////////
// input: string, for file name
// output: PyArrayObject, as volumetric data
static PyObject * read_cube(PyObject * self, PyObject * args){
  char *input; // filename string as input
  double *data; // C-data need to be converted to numpy object
  double *structure; // C-data need to be converted to numpy object
  double *grid; // C-data need to be converted to numpy object
  int dims[3], size[2], gsize[2]={4,4};
  int N3=1, i;
  // numpy object pass back to python
  PyArrayObject *npdata, *npstructure, *npgrid; 

  // parse arguments check and/or error handling
  if(!PyArg_ParseTuple(args, "s", &input))
    return NULL;

  // read CUBE file header for memory size
  readCubeHeader(input, dims, size);

  // to allocate memory in the function
  for(i=0;i<3;i++) N3 *= dims[i];
  data = (double*) malloc(N3*sizeof(double));
  structure = (double*) malloc(4*size[0]*sizeof(double));
  grid = (double*) malloc(16*sizeof(double));
  // set up numpy array data structure
  npdata = (PyArrayObject*) PyArray_FromDims(3, dims, NPY_DOUBLE);
  npstructure = (PyArrayObject*) PyArray_FromDims(2, 
                                   size, NPY_DOUBLE);
  npgrid = (PyArrayObject*) PyArray_FromDims(2, gsize, NPY_DOUBLE);
  // set numpy data point to C-array data
  data = pyvector_to_Carrayptrs(npdata);
  structure = pyvector_to_Carrayptrs(npstructure);
  grid = pyvector_to_Carrayptrs(npgrid);
  // run the actual function
  // NOTE the data is passed as datatype double**
  readcube_c(input, data, structure, grid, dims, size);

  // build return tuple of PyObjects
  // NOTE: data which Numpy array points to CANNOT be freeed
  return Py_BuildValue("OOO", npdata, npstructure, npgrid);
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
