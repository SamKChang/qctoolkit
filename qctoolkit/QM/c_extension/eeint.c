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
static PyObject* eeint(PyObject* self, PyObject* args){

  /*** input variables ***/
  /* dictionary */
  PyObject *in_dict;
  PyObject *dictList;
  PyObject *dictList_fast;
  PyObject *pyStr;
  PyObject *item;
  /* numpy array */
  PyArrayObject *in_array1, *in_array2;
  NpyIter *in_iter;
  NpyIter_IterNextFunc *in_iternext;
  /* C data type for input */
  int Ngo, Nao;
  int N, itr;
  /* basis function variables */
  double **in_ptr; // address to numpy pointer
  double *center;  // gaussian center
  double *cef;     // contraction coefficients
  int *ng;         // number of gaussians per AO
  double *exp;     // gaussian exponents
  int *lm_xyz;     // angular momentum

  /* python output variables */
  PyObject *py_out;
  PyObject *py_out2;
  double *data, *overlap;
  double element;
  int i, j, k, l;
  int s, t, u, v, w;
  int mat_dim[4];

  /*  parse numpy array and two integers as argument */
  if (!PyArg_ParseTuple(args, "OO!O!",
                        &in_dict, 
                        &PyArray_Type, &in_array1,
                        &PyArray_Type, &in_array2
                       )) return NULL;
  if(in_array1 == NULL) return NULL;

  /***********************
  * Construct input data *
  ***********************/
  /*** access dict_list data ***/
  /* exponents */
  pyStr = PyString_FromString("exponents");
  dictList = PyDict_GetItem(in_dict, pyStr);
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
  dictList = PyDict_GetItem(in_dict, pyStr);
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
  dictList = PyDict_GetItem(in_dict, pyStr);
  dictList_fast = PySequence_Fast(
                    dictList, "expected a sequence");
  Nao = PySequence_Size(dictList);
  ng = (int*) malloc(Nao * sizeof(int));
  for(i=0;i<Nao;i++){
    item = PySequence_Fast_GET_ITEM(dictList_fast, i);
    ng[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(dictList_fast);
  /*** end of dict_list data ***/

  /*** create the iterators, necessary to access numpy array ***/
  /* center */
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

  /*** end of numpy input data ***/

  /***** end of input data construction *****/


  /*******************
  * construct output *
  *******************/
  for(i=0;i<4;i++) mat_dim[i] = Nao;
  py_out = (PyArrayObject*) 
           PyArray_FromDims(4, mat_dim, NPY_DOUBLE);
  data = pyvector_to_Carrayptrs(py_out);

  /* renormalization */
  renormalize(center, exp, cef, ng, lm_xyz, Nao);

  /* orthogonalize */
  // no effect at the moment, for debug purpose
  //py_out2 = (PyArrayObject*) 
  //          PyArray_FromDims(2, mat_dim, NPY_DOUBLE);
  //overlap = pyvector_to_Carrayptrs(py_out2);
  //orthogonalize(overlap, center, exp, cef, ng, lm_xyz, Nao);

#pragma omp parallel private(i, j, k, l, s, t, u, v, w)\
shared(data)
{
  #pragma omp for schedule(dynamic)
  for(i=0;i<Nao;i++){
    for(j=i;j<Nao;j++){
      for(k=0;k<Nao;k++){
        for(l=k;l<Nao;l++){
          if(l+k >= i+j){
            s = l + k*Nao + j*Nao*Nao + i*Nao*Nao*Nao;
            element = eeMatrix(center, exp, cef, ng, lm_xyz,
                                Nao, i, j, k, l);
            data[s] = element;

            // symmetry for (ij|kl)=(ij|lk)=(ji|kl)=(ji|lk)
            t = k + l*Nao + j*Nao*Nao + i*Nao*Nao*Nao; //(ij|lk)
            u = l + k*Nao + i*Nao*Nao + j*Nao*Nao*Nao; //(ji|kl)
            v = k + l*Nao + i*Nao*Nao + j*Nao*Nao*Nao; //(ji|lk)
            data[t] = data[s];
            data[u] = data[s];
            data[v] = data[s];
    
            // symmetry for (ij|kl)=(kl|ij)=(kl|ji)=(lk|ij)=(lk|ji)
            t = j + i*Nao + l*Nao*Nao + k*Nao*Nao*Nao; //(kl|ij)
            u = i + j*Nao + l*Nao*Nao + k*Nao*Nao*Nao; //(kl|ji)
            v = j + i*Nao + k*Nao*Nao + l*Nao*Nao*Nao; //(lk|ij)
            w = i + j*Nao + k*Nao*Nao + l*Nao*Nao*Nao; //(lk|ji)
            data[t] = data[s];
            data[u] = data[s];
            data[v] = data[s];
            data[w] = data[s];
          }
        }
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
  //free(overlap);

  Py_INCREF(py_out);
  return Py_BuildValue("O", py_out);

  /*  in case bad things happen */
  fail:
    Py_XDECREF(py_out);
    return NULL;
} // end of eeint function

/*  define functions in module */
static PyMethodDef EEInt[] ={
  {"eeint", eeint, METH_VARARGS,
      "analytic two-center Gaussian-Coulomb integral"},
  {NULL, NULL, 0, NULL}
};

/* module initialization */
PyMODINIT_FUNC initeeint(void){
  (void) Py_InitModule("eeint", EEInt);
  /* IMPORTANT: this must be called */
  import_array();
}
