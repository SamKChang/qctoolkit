// C routine to calculate ESP for a point

#include<Python.h>
#include <numpy/arrayobject.h>

/////////////////////////////////////////
// main function for writing CUBE file //
/////////////////////////////////////////
// or what ever computation need to be done
// only minor changes should be necessary
// to interface other type of calculation
static double esp_point_c(double *grid, 
                 double *structure, 
                 double *data,
                 int Ndim[4],
                 int N3,
                 double x,
                 double y,
                 double z)
{
  FILE *fr;
  int i, j, k, s, I;
  double X, Y, Z, r, Q, V = 0;
  double XI, YI, ZI, QI;
  double dV;
  int num;


  dV = grid[5] * grid[10] * grid[15];
  for(i=0;i<Ndim[1];i++){
    X = grid[1] + i * grid[5];
    for(j=0;j<Ndim[2];j++){
      Y = grid[2] + i * grid[10];
      for(k=0;k<Ndim[3];k++){
        Z = grid[3] + i * grid[15];
        s = k + j*Ndim[3] + i*Ndim[2]*Ndim[3];
        Q = data[s] * dV;
        r = sqrt(pow(X-x, 2) + pow(Y-y, 2) + pow(Z-z, 2));
        V += Q / r;
      }
    }
  }
  for(I=0;I<grid[0];I++){
    Q = -structure[4*I];
    X = structure[4*I + 1];
    Y = structure[4*I + 2];
    Z = structure[4*I + 3];
    r = sqrt(pow(X-x, 2) + pow(Y-y, 2) + pow(Z-z, 2));
    V += Q / r;
  }
  return V;
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
static PyObject * esp_point(PyObject * self, PyObject * args){
  PyObject *py_list_holder; // variable to pass python to C
  double x, y, z; // coordinate from python
  PyArrayObject *py_data; // numpy data from python
  PyArrayObject *py_structure; // numpy data from python
  PyArrayObject *py_grid; // numpy data from python
  double *data; // C-data need to be converted from python object
  double *structure; // C-data need to be converted from python object
  double *grid; // C-data need to be converted from python object
  double V; // resaulting ESP at (x, y, z)
  PyObject *pyV; // resaulting ESP to python
  int Ndim[4], N3, Ns;
  int i;

  // parse arguments check and/or error handling
  if(!PyArg_ParseTuple(args, "O!O!O!ddd", 
                       &PyArray_Type, &py_grid,
                       &PyArray_Type, &py_structure,
                       &PyArray_Type, &py_data,
                       &x, &y, &z
                      )) return NULL;
  if((py_grid == NULL)||(py_structure == NULL)||(py_data == NULL))
    return NULL;

  /*** access numpy array data ***/
  grid = (double *) malloc(16 * sizeof(double));
  nplist_to_C_double_array(py_grid, grid, 16);
  for(i=0;i<4;i++) Ndim[i] = grid[i*4];
  N3 = Ndim[1] * Ndim[2] * Ndim[3];
  Ns = Ndim[0] * 4;
  structure = (double *) malloc(Ns * sizeof(double));
  data = (double *) malloc(N3 * sizeof(double));
  nplist_to_C_double_array(py_structure, structure, Ns);
  nplist_to_C_double_array(py_data, data, N3);
  /*** end of list input data ***/

  V = esp_point_c(grid, structure, data, Ndim, N3, x, y, z);
  free(grid);
  free(structure);
  free(data);

  return Py_BuildValue("d", V);
}

/////////////////////////////////////////
// register python module symbol table //
/////////////////////////////////////////
// PyMethodDef: struct of four field
//   string: method_name
//   PyCFunction: method_function
//   int: flag
//   string: documentation
static PyMethodDef ESPPointMethods[] = {
  {"esp_point", esp_point, METH_VARARGS, "ESP at a point"},
  {NULL, NULL, 0, NULL} // sentinel?
};

//////////////////////////
//  method constructor  //
//////////////////////////
// the name MUST be init{name} otherwise python cannot find it
// depends on numpy C-API, import_array()
// and/or import_ufunc() are necessary
// otherwise the code return segfault
PyMODINIT_FUNC initesp_point(void){
  Py_InitModule("esp_point", ESPPointMethods);
  import_array(); // necessary!
}
