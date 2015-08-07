#include <Python.h>
#include <numpy/arrayobject.h>
#include "kernels.h"
#include <omp.h>

/* extract row vector with index='index' in 'data' matrix */
void getVector(double *vector, double *data, int index, int max){
  int itr;
  for(itr=0;itr<max;itr++){
    vector[itr] = data[index * max + itr];
  }
}

/*  python interface */
static PyObject* kernel_vectors(PyObject* self, PyObject* args){

  /* input variables */
  PyArrayObject *ref_array;
  NpyIter *ref_iter;
  NpyIter_IterNextFunc *ref_iternext;
  PyArrayObject *tar_array;
  NpyIter *tar_iter;
  NpyIter_IterNextFunc *tar_iternext;
  PyObject *klargs;
  PyObject *seq;
  PyObject *item;
  double *rdata, *tdata;
  char *kernel = (char*) malloc(20);
  double *kerargs;
  int rrows, rcolumns, trows, tcolumns, len;

  /* python output variables */
  PyObject *np_matrix;
  double *matrix;
  double *xi, *xj;
  int vec_dim[1];
  int mat_dim[2];
  int i, j, itr = 0;

  /*  parse numpy array and two integers as argument */
  if (!PyArg_ParseTuple(args, "O!iiO!iis|O", 
                        &PyArray_Type, &ref_array,// O!, NpArray
                        &rrows,                   // i
                        &rcolumns,                // i
                        &PyArray_Type, &tar_array,// O!, NpArray
                        &trows,                   // i
                        &tcolumns,                // i
                        &kernel,                  // s, kernel
      /* kernel args */ &klargs)) return NULL;
  
  /***********************
  * Construct input data *
  ***********************/
  rdata = (double*) malloc(rrows * rcolumns * sizeof(double));
  tdata = (double*) malloc(trows * tcolumns * sizeof(double));
  /*  create the iterators, necessary to access numpy array */
  ref_iter = NpyIter_New(ref_array, NPY_ITER_READONLY, 
                         NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  tar_iter = NpyIter_New(tar_array, NPY_ITER_READONLY, 
                         NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  /* check input status */
  if (ref_iter == NULL || tar_iter == NULL)
      goto fail;
  ref_iternext = NpyIter_GetIterNext(ref_iter, NULL);
  tar_iternext = NpyIter_GetIterNext(tar_iter, NULL);
  if (ref_iternext == NULL || tar_iternext == NULL) {
      NpyIter_Deallocate(ref_iter);
      NpyIter_Deallocate(tar_iter);
      goto fail;
  }
  /* interator pointer to actual numpy data */
  double **ref_dataptr = (double **) 
                           NpyIter_GetDataPtrArray(ref_iter);
  double **tar_dataptr = (double **) 
                           NpyIter_GetDataPtrArray(tar_iter);
  /*  iterate over the arrays */
  itr=0;
  do {
      rdata[itr++] = **ref_dataptr;
  } while(ref_iternext(ref_iter));
  itr=0;
  do {
      tdata[itr++] = **tar_dataptr;
  } while(tar_iternext(tar_iter));

  /* access python list data */
  seq = PySequence_Fast(klargs, "expected a sequence");
  len = PySequence_Size(klargs);
  kerargs = (double*) malloc(len * sizeof(double));
  for(i=0;i<len;i++){
    item = PySequence_Fast_GET_ITEM(seq, i);
    kerargs[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(seq);
  /***** end of numpy array construction *****/

  /**************************
  * construct output matrix *
  **************************/
  matrix = (double *) malloc(rrows * trows * sizeof(double));
  mat_dim[0] = trows;
  mat_dim[1] = rrows;
  if (rcolumns != tcolumns){
    printf("data dimension not fit\n");
    return NULL;
  }
  vec_dim[0] = rcolumns;

#pragma omp parallel private(i,j,xi,xj) shared(matrix)
{
  #pragma omp for
  for(i=0;i<trows;i++){
    for(j=0;j<rrows;j++){
      xi = (double*) malloc(tcolumns * sizeof(double));
      xj = (double*) malloc(rcolumns * sizeof(double));
      getVector(xi, tdata, i, tcolumns);
      getVector(xj, rdata, j, rcolumns);
      matrix[j+rrows*i] = kerEval(kernel,kerargs,xi,xj,rcolumns);
      free(xi);
      free(xj);
    }
  }
}
  
  np_matrix = PyArray_SimpleNewFromData(2, mat_dim, 
                                        NPY_DOUBLE, matrix);
  /***** end of output matrix construction *****/

  /*********************************
  * clean up and return the result *
  *********************************/
  Py_INCREF(np_matrix);
  return np_matrix;
  NpyIter_Deallocate(ref_iter);
  NpyIter_Deallocate(tar_iter);

  /*  in case bad things happen */
  fail:
      Py_XDECREF(np_matrix);
      return NULL;

  free(rdata);
  free(tdata);
  free(matrix);
}

/*  define functions in module */
static PyMethodDef KernelVectors[] ={
  {"kernel_vectors", kernel_vectors, METH_VARARGS,
      "evaluate the cosine on a numpy array"},
  {NULL, NULL, 0, NULL}
};

/* module initialization */
PyMODINIT_FUNC initkernel_vectors(void){
  (void) Py_InitModule("kernel_vectors", KernelVectors);
  /* IMPORTANT: this must be called */
  import_array();
}
