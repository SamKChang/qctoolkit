// C routine to calculate ESP for a point

#include<Python.h>
#include <numpy/arrayobject.h>
#include <omp.h>

/////////////////////////////////////////
// main function for writing CUBE file //
/////////////////////////////////////////
// or what ever computation need to be done
// only minor changes should be necessary
// to interface other type of calculation
static esp_cube_c(double *grid, 
                  double *structure, 
                  double *data,
                  int Ndim[4],
                  int N3,
                  double *data_out)
{
  int i, j, k, s, I;
  int a, b, c, t;
  double x, y, z;
  double X, Y, Z, r, Q, V = 0;
  double XI, YI, ZI, QI;
  double dV;

  dV = grid[5] * grid[10] * grid[15];
#pragma omp parallel private(a,b,c,x,y,z,t,i,j,k,X,Y,Z,s,Q,r,I) \
shared(data_out)
{
  #pragma omp for
  for(a=0;a<Ndim[1];a++){
    x = grid[1] + a * grid[5];
    for(b=0;b<Ndim[2];b++){
      y = grid[2] + b * grid[10];
      for(c=0;c<Ndim[3];c++){
        // return cube file initialized at new point
        z = grid[3] + c * grid[15];
        t = c + b*Ndim[3] + a*Ndim[2]*Ndim[3];
        data_out[t] = 0;
        // loop through all data point
        for(i=0;i<Ndim[1];i++){
          X = grid[1] + i * grid[5];
          for(j=0;j<Ndim[2];j++){
            Y = grid[2] + j * grid[10];
            for(k=0;k<Ndim[3];k++){
              Z = grid[3] + k * grid[15];
              s = k + j*Ndim[3] + i*Ndim[2]*Ndim[3];
              Q = data[s] * dV;
              r = sqrt(pow(X-x, 2) + pow(Y-y, 2) + pow(Z-z, 2));
              if(r > 1) data_out[t] += Q / r;
            }
          }
        }
        // loop through all atoms
        for(I=0;I<grid[0];I++){
          Q = -structure[4*I];
          X = structure[4*I + 1];
          Y = structure[4*I + 2];
          Z = structure[4*I + 3];
          r = sqrt(pow(X-x, 2) + pow(Y-y, 2) + pow(Z-z, 2));
          if(r > 1) data_out[t] += Q / r;
        }
      }
    }
  }
} 
// end omp for
}

static void nplist_to_C_double_array(PyArrayObject *in_array, 
                              double *out_list,
                              int N)
{
  NpyIter *in_iter;
  NpyIter_IterNextFunc *in_iternext;
  double **in_ptr;
  int itr;

  in_iter = NpyIter_New(in_array, NPY_ITER_READONLY,
                        NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  if (in_iter == NULL) goto fail;
  in_iternext = NpyIter_GetIterNext(in_iter, NULL);
  if (in_iternext == NULL){
    NpyIter_Deallocate(in_iter);
    goto fail;
  }
  /* interator pointer to actual numpy data */
  in_ptr = (double **) NpyIter_GetDataPtrArray(in_iter);
  itr=0;
  do {
    out_list[itr++] = **in_ptr;
  } while(in_iternext(in_iter));
  NpyIter_Deallocate(in_iter);

  fail:
    return NULL;
}

//////////////////////////////////////////////
// redirect PyArray data contant to C-array //
//////////////////////////////////////////////
static double *pyvector_to_Carrayptrs(PyArrayObject *arrayin){
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
// input: grid, structure, data 
// output: none
static PyObject * esp_cube(PyObject * self, PyObject * args){
  PyObject *py_list_holder; // variable to pass python to C
  double x, y, z; // coordinate from python
  PyArrayObject *py_data; // numpy data from python
  PyArrayObject *py_structure; // numpy data from python
  PyArrayObject *py_grid; // numpy data from python
  PyArrayObject *py_data_out; // numpy data to return
  double *data; // C-data need to be converted from python object
  double *structure; // C-data need to be converted from python object
  double *grid; // C-data need to be converted from python object
  double *data_out; // C-data for output
  double V; // resaulting ESP at (x, y, z)
  int Ndim[4], N3, Ns;
  int i;
  int dims[3];

  // parse arguments check and/or error handling
  if(!PyArg_ParseTuple(args, "O!O!O!", 
                       &PyArray_Type, &py_grid,
                       &PyArray_Type, &py_structure,
                       &PyArray_Type, &py_data
                      )) return NULL;
  if((py_grid == NULL)||(py_structure == NULL)||(py_data == NULL))
    return NULL;

  /*** access numpy array data ***/
  grid = (double *) malloc(16 * sizeof(double));
  nplist_to_C_double_array(py_grid, grid, 16);
  for(i=0;i<4;i++) Ndim[i] = grid[i*4];
  for(i=1;i<4;i++) dims[i-1] = Ndim[i];
  N3 = Ndim[1] * Ndim[2] * Ndim[3];
  Ns = Ndim[0] * 4;
  structure = (double *) malloc(Ns * sizeof(double));
  data = (double *) malloc(N3 * sizeof(double));
  data_out = (double *) malloc(N3 * sizeof(double));
  nplist_to_C_double_array(py_structure, structure, Ns);
  nplist_to_C_double_array(py_data, data, N3);

  py_data_out = (PyArrayObject*) PyArray_FromDims(3, dims, NPY_DOUBLE);
  data_out = pyvector_to_Carrayptrs(py_data_out);
  /*** end of list input data ***/

  esp_cube_c(grid, structure, data, Ndim, N3, data_out);

  //free(grid);
  //free(structure);
  //free(data);

  return Py_BuildValue("O", py_data_out);
}

/////////////////////////////////////////
// register python module symbol table //
/////////////////////////////////////////
// PyMethodDef: struct of four field
//   string: method_name
//   PyCFunction: method_function
//   int: flag
//   string: documentation
static PyMethodDef ESPCubeMethods[] = {
  {"esp_cube", esp_cube, METH_VARARGS, "ESP cube from density"},
  {NULL, NULL, 0, NULL} // sentinel?
};

//////////////////////////
//  method constructor  //
//////////////////////////
// the name MUST be init{name} otherwise python cannot find it
// depends on numpy C-API, import_array()
// and/or import_ufunc() are necessary
// otherwise the code return segfault
PyMODINIT_FUNC initesp_cube(void){
  Py_InitModule("esp_cube", ESPCubeMethods);
  import_array(); // necessary!
}
