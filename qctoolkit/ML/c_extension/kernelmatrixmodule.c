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
static PyObject* kernel_matrix(PyObject* self, PyObject* args){

  /* input variables */
  PyArrayObject *in_array;
  NpyIter *in_iter;
  NpyIter_IterNextFunc *in_iternext;
  PyObject *klargs;
  PyObject *seq;
  PyObject *item;
  double *data;
  char *kernel = (char*) malloc(20);
  double *kerargs;
  int rows, columns, len;

  /* python output variables */
  PyObject *np_matrix;
  double *matrix;
  double *xi, *xj;
  int vec_dim[1];
  int mat_dim[2];
  int i, j, itr = 0;

  /*  parse numpy array and two integers as argument */
  if (!PyArg_ParseTuple(args, "O!iis|O", 
                        &PyArray_Type, &in_array, // O!, NpArray
                        &rows,                    // i
                        &columns,                 // i
                        &kernel,                  // s, kernel
     /* kernel args */  &klargs)) return NULL;

  /***********************
  * Construct input data *
  ***********************/
  data = (double*) malloc(rows * columns * sizeof(double));
  /*  create the iterators, necessary to access numpy array */
  in_iter = NpyIter_New(in_array, NPY_ITER_READONLY, 
                        NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  /* check input status */
  if (in_iter == NULL)
      goto fail;
  in_iternext = NpyIter_GetIterNext(in_iter, NULL);
  if (in_iternext == NULL) {
      NpyIter_Deallocate(in_iter);
      goto fail;
  }
  /* interator pointer to actual numpy data */
  double **in_dataptr = (double **) 
                           NpyIter_GetDataPtrArray(in_iter);
  /*  iterate over the arrays */
  do {
      data[itr++] = **in_dataptr;
  } while(in_iternext(in_iter));

  /* access python list data */
  seq = PySequence_Fast(klargs, "expected a sequence");
  len = PySequence_Size(klargs);
  kerargs = (double*) malloc(len * sizeof(double));
  for(i=0;i<len;i++){
    item = PySequence_Fast_GET_ITEM(seq, i);
    kerargs[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(seq);
  /***** end of input data construction *****/

  /**************************
  * construct output matrix *
  **************************/
  matrix = (double *) malloc(rows * rows * sizeof(double));
  mat_dim[0] = rows;
  mat_dim[1] = rows;
  vec_dim[0] = columns;

#pragma omp parallel private(i,j,xi,xj) shared(matrix)
{
  #pragma omp for
  for(i=0;i<rows;i++){
    for(j=i;j<rows;j++){
      xi = (double*) malloc(columns * sizeof(double));
      xj = (double*) malloc(columns * sizeof(double));
      /* construct callback input numpy arrays*/
      getVector(xi, data, i, columns);
      getVector(xj, data, j, columns);
      matrix[j+rows*i] = kerEval(kernel,kerargs,xi,xj,columns);
      if(i!=j) matrix[i+rows*j] = matrix[j+rows*i];
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
  NpyIter_Deallocate(in_iter);

  /*  in case bad things happen */
  fail:
      Py_XDECREF(np_matrix);
      return NULL;

  free(data);
  free(matrix);
}

/*  define functions in module */
static PyMethodDef KernelMatrix[] ={
  {"kernel_matrix", kernel_matrix, METH_VARARGS,
      "evaluate the cosine on a numpy array"},
  {NULL, NULL, 0, NULL}
};

/* module initialization */
PyMODINIT_FUNC initkernel_matrix(void){
  (void) Py_InitModule("kernel_matrix", KernelMatrix);
  /* IMPORTANT: this must be called */
  import_array();
}
