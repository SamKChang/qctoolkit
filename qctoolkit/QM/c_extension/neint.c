#include <Python.h>
#include <numpy/arrayobject.h>
#include <omp.h>
#include "gaussian.h"

/**********************************
*  necessary Numpy functionality  *
**********************************/
// return 1D memory array pointer to allocated numpy array 
double* pyvector_to_Carrayptrs(PyArrayObject *arrayin){
  int n = arrayin->dimensions[0];
  return (double*) arrayin->data;
}

/*********************
*  python interface  *
*********************/
static PyObject* neint(PyObject* self, PyObject* args){

  /*** input variables ***/
  /* dictionary */
  PyObject *in_dict1, *in_dict2;
  PyObject *dictList;
  PyObject *dictList_fast;
  PyObject *pyStr;
  PyObject *item;
  /* numpy array */
  PyArrayObject *in_array1, *in_array2, *in_array3, *in_array4;
  NpyIter *in_iter;
  NpyIter_IterNextFunc *in_iternext;
  /* C data type for input */
  int Ngo, Nao;    // number of gaussian and atomic orbiatls
  int fNgo, fNao;  // number of gaussian and atomic densities
  int N, itr;
  /* basis function variables */
  double **in_ptr; // address to numpy pointer
  double **fin_ptr;// address to numpy pointer
  double *center;  // gaussian center
  double *fcenter; // gaussian center, for density-fitting
  double *cef;     // contraction coefficients
  double *fcef;    // contraction coefficients, for density-fitting
  int *ng;         // number of gaussians per AO
  int *fng;        // number of gaussians per atomic density
  double *exp;     // gaussian exponents
  double *fexp;    // gaussian exponents
  int *lm_xyz;     // angular momentum
  int *flm_xyz;    // angular momentum, for density-fitting

  /* python output variables */
  PyObject *py_out;
  //PyObject *py_out2;
  double *data, *overlap;
  double element;
  int s, a, i, j;
  int mat_dim[3];

  /*  parse numpy array and two integers as argument */
  if (!PyArg_ParseTuple(args, "OO!O!OO!O!",
                        &in_dict1, 
                        &PyArray_Type, &in_array1,
                        &PyArray_Type, &in_array2,
                        &in_dict2, 
                        &PyArray_Type, &in_array3,
                        &PyArray_Type, &in_array4
                       )) return NULL;
  if(in_array1 == NULL) return NULL;

  /***********************
  * Construct input data *
  ***********************/
  /*** access orbital dict_list1 data ***/
  /* exponents */
  pyStr = PyString_FromString("exponents");
  dictList = PyDict_GetItem(in_dict1, pyStr);
  dictList_fast = PySequence_Fast(
                    dictList, "expected a sequence");
  Ngo = PySequence_Size(dictList);
  exp = (double*) malloc(Ngo * sizeof(double));
  for(i=0;i<Ngo;i++){
    item = PySequence_Fast_GET_ITEM(dictList_fast, i);
    exp[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(dictList_fast);
  /* coefficients */
  pyStr = PyString_FromString("coefficients");
  dictList = PyDict_GetItem(in_dict1, pyStr);
  dictList_fast = PySequence_Fast(
                    dictList, "expected a sequence");
  cef = (double*) malloc(Ngo * sizeof(double));
  for(i=0;i<Ngo;i++){
    item = PySequence_Fast_GET_ITEM(dictList_fast, i);
    cef[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(dictList_fast);
  /* number of gaussians per shell */
  pyStr = PyString_FromString("n_gaussians");
  dictList = PyDict_GetItem(in_dict1, pyStr);
  dictList_fast = PySequence_Fast(
                    dictList, "expected a sequence");
  Nao = PySequence_Size(dictList);
  ng = (int*) malloc(Nao * sizeof(int));
  for(i=0;i<Nao;i++){
    item = PySequence_Fast_GET_ITEM(dictList_fast, i);
    ng[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(dictList_fast);
  /*** end of orbital dict_list1 data ***/

  /*** access density dict_list2 data ***/
  /* exponents */
  pyStr = PyString_FromString("exponents");
  dictList = PyDict_GetItem(in_dict2, pyStr);
  dictList_fast = PySequence_Fast(
                    dictList, "expected a sequence");
  fNgo = PySequence_Size(dictList);
  fexp = (double*) malloc(fNgo * sizeof(double));
  for(i=0;i<fNgo;i++){
    item = PySequence_Fast_GET_ITEM(dictList_fast, i);
    fexp[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(dictList_fast);
  /* coefficients */
  pyStr = PyString_FromString("coefficients");
  dictList = PyDict_GetItem(in_dict2, pyStr);
  dictList_fast = PySequence_Fast(
                    dictList, "expected a sequence");
  fcef = (double*) malloc(fNgo * sizeof(double));
  for(i=0;i<fNgo;i++){
    item = PySequence_Fast_GET_ITEM(dictList_fast, i);
    fcef[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(dictList_fast);
  /* number of gaussians per shell */
  pyStr = PyString_FromString("n_gaussians");
  dictList = PyDict_GetItem(in_dict2, pyStr);
  dictList_fast = PySequence_Fast(
                    dictList, "expected a sequence");
  fNao = PySequence_Size(dictList);
  fng = (int*) malloc(fNao * sizeof(int));
  for(i=0;i<fNao;i++){
    item = PySequence_Fast_GET_ITEM(dictList_fast, i);
    fng[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(dictList_fast);
  /*** end of density dict_list2 data ***/

  /*** Access numpy array for orbital ***/
  /* center */
  // create the iterators, necessary to access numpy array
  in_iter = NpyIter_New(in_array1, NPY_ITER_READONLY,
                         NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  if (in_iter == NULL) goto fail;
  in_iternext = NpyIter_GetIterNext(in_iter, NULL);
  if (in_iternext == NULL){
    NpyIter_Deallocate(in_iter);
    goto fail;
  }
  /* interator pointer to actual numpy data */
  in_ptr = (double **) NpyIter_GetDataPtrArray(in_iter);
  center = (double*) malloc(3 * Nao * sizeof(double));
  itr=0;
  do {
    center[itr++] = **in_ptr;
  } while(in_iternext(in_iter));
  NpyIter_Deallocate(in_iter);

  /* lm_xyz */
  in_iter = NpyIter_New(in_array2, NPY_ITER_READONLY,
                         NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  if (in_iter == NULL) goto fail;
  in_iternext = NpyIter_GetIterNext(in_iter, NULL);
  if (in_iternext == NULL){
    NpyIter_Deallocate(in_iter);
    goto fail;
  }
  /* interator pointer to actual numpy data */
  int **in_ptr2 = (int **) NpyIter_GetDataPtrArray(in_iter);
  lm_xyz = (int*) malloc(3 * Nao * sizeof(int));
  itr=0;
  do {
    lm_xyz[itr++] = **in_ptr2;
  } while(in_iternext(in_iter));
  NpyIter_Deallocate(in_iter);
  /*** end of orbital numpy input data ***/

  /*** Access numpy array for orbital ***/
  /* fcenter */
  // create the iterators, necessary to access numpy array
  in_iter = NpyIter_New(in_array3, NPY_ITER_READONLY,
                         NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  if (in_iter == NULL) goto fail;
  in_iternext = NpyIter_GetIterNext(in_iter, NULL);
  if (in_iternext == NULL){
    NpyIter_Deallocate(in_iter);
    goto fail;
  }
  /* interator pointer to actual numpy data */
  in_ptr = (double **) NpyIter_GetDataPtrArray(in_iter);
  fcenter = (double*) malloc(3 * Nao * sizeof(double));
  itr=0;
  do {
    fcenter[itr++] = **in_ptr;
  } while(in_iternext(in_iter));
  NpyIter_Deallocate(in_iter);

  /* flm_xyz */
  in_iter = NpyIter_New(in_array4, NPY_ITER_READONLY,
                         NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  if (in_iter == NULL) goto fail;
  in_iternext = NpyIter_GetIterNext(in_iter, NULL);
  if (in_iternext == NULL){
    NpyIter_Deallocate(in_iter);
    goto fail;
  }
  /* interator pointer to actual numpy data */
  int **fin_ptr2 = (int **) NpyIter_GetDataPtrArray(in_iter);
  flm_xyz = (int*) malloc(3 * Nao * sizeof(int));
  itr=0;
  do {
    flm_xyz[itr++] = **fin_ptr2;
  } while(in_iternext(in_iter));
  NpyIter_Deallocate(in_iter);
  /*** end of orbital numpy input data ***/

  /***** end of input data construction *****/


  /*******************
  * construct output *
  *******************/
  for(i=1;i<3;i++) mat_dim[i] = Nao;
  mat_dim[0] = fNao;
  py_out = (PyArrayObject*) 
           PyArray_FromDims(3, mat_dim, NPY_DOUBLE);
  data = pyvector_to_Carrayptrs(py_out);

  /* renormalization of orbital */
  renormalize(center, exp, cef, ng, lm_xyz, Nao);
  /* renormalization of density */
  densityRenormalize(fcenter, fexp, fcef, fng, flm_xyz, fNao);

  /* orthogonalize */
  // no effect at the moment, for debug purpose
  //py_out2 = (PyArrayObject*) 
  //          PyArray_FromDims(2, mat_dim, NPY_DOUBLE);
  //overlap = pyvector_to_Carrayptrs(py_out2);
  //orthogonalize(overlap, center, exp, cef, ng, lm_xyz, Nao);

#pragma omp parallel private(a, i, j, s, element) shared(data)
{
  #pragma omp for schedule(dynamic)
  for(a=0;a<fNao;a++){
    for(i=0;i<Nao;i++){
      for(j=i;j<Nao;j++){
        s = j + i*Nao + a*Nao*Nao;
        element = neMatrix(center, exp, cef, ng, lm_xyz, Nao,
                           fcenter, fexp, fcef, fng, flm_xyz, fNao,
                           a, i, j);
        data[s] = element;
        if(j>i) data[i + j*Nao + a*Nao*Nao] = element;
      }
    }
  }
} // end of omp loop

  /*********************************
  * clean up and return the result *
  *********************************/
  free(exp);
  free(cef);
  free(ng);
  free(center);
  free(lm_xyz);
  free(fexp);
  free(fcef);
  free(fng);
  free(fcenter);
  free(flm_xyz);
  //free(overlap);

  Py_INCREF(py_out);
  return Py_BuildValue("O", py_out);

  /*  in case bad things happen */
  fail:
    Py_XDECREF(py_out);
    return NULL;
} // end of neint function

/*  define functions in module */
static PyMethodDef NEInt[] ={
  {"neint", neint, METH_VARARGS,
   "analytic Gaussian-Coulomb integral for n-psi^2"},
  {NULL, NULL, 0, NULL}
};

/* module initialization */
PyMODINIT_FUNC initneint(void){
  (void) Py_InitModule("neint", NEInt);
  /* IMPORTANT: this must be called */
  import_array();
}
