/*  Example of wrapping the cos function from math.h using the Numpy-C-API. */

#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
//#include <omp.h>

/* python callable pointer */
static PyObject *my_callback = NULL;

/* extract row vector with index='index' in 'data' matrix */
void getVector(double **vector, double *data, int index, int max){
  int itr;
  free(*vector); /* free previously allocated memory */
  *vector = (double*) malloc(max * sizeof(double));
  for(itr=0;itr<max;itr++){
    (*vector)[itr] = data[index * max + itr];
  }
}

void getMatrixElement(
                      int i,
                      int j,
                      int rrows,
                      int columns,
                      double *tdata,
                      double *rdata,
                      int vec_dim[],
                      double *matrix,
                      PyObject *callback
                     )
{
  int s = j+rrows*i;
  double *xi = (double*) malloc(sizeof(double));
  double *xj = (double*) malloc(sizeof(double));
  PyObject *np_xi, *np_xj;
  PyObject *arglist;
  PyObject *element;

  /* construct callback input numpy arrays*/
  getVector(&xi, tdata, i, columns);
  getVector(&xj, rdata, j, columns);
  
  np_xi = PyArray_SimpleNewFromData(1, vec_dim, 
                                    NPY_DOUBLE,xi);
  np_xj = PyArray_SimpleNewFromData(1, vec_dim, 
                                    NPY_DOUBLE,xj);
  /* construct callback argument list */
  arglist = Py_BuildValue("(OO)", np_xi, np_xj);
  /* simple python call with return PyObject */
  element = PyObject_CallObject(callback, arglist);
  /* extract 'double' data type from PyObject */
  matrix[s] = PyFloat_AsDouble(element);

  Py_DECREF(arglist);
  Py_DECREF(element);
  Py_DECREF(np_xi);
  Py_DECREF(np_xj);
}

/*  python interface */
static PyObject* kernel_vectors(PyObject* self, PyObject* args){

  int itr = 0;

  /* numpy input variables */
  PyArrayObject *ref_array;
  NpyIter *ref_iter;
  NpyIter_IterNextFunc *ref_iternext;
  PyArrayObject *tar_array;
  NpyIter *tar_iter;
  NpyIter_IterNextFunc *tar_iternext;
  double *rdata, *tdata;
  int rrows, rcolumns, trows, tcolumns;

  /* python callback variables */
  PyObject *result = NULL;
  PyObject *py_obj;
  PyObject *py_objMethod;
  PyObject *arglist;
  PyObject *np_xi, *np_xj;
  double *xi = (double *) malloc(sizeof(double));
  double *xj = (double *) malloc(sizeof(double));
  int vec_dim[1];

  /* python output variables */
  PyObject *np_matrix;
  PyObject *element;
  double *matrix;
  int mat_dim[2];
  int i, j;

  /*  parse numpy array and two integers as argument */
  if (!PyArg_ParseTuple(args, "O!iiO!iiO", 
                        &PyArray_Type, &ref_array,// O!, NpArray
                        &rrows,                   // i
                        &rcolumns,                // i
                        &PyArray_Type, &tar_array,// O!, NpArray
                        &trows,                   // i
                        &tcolumns,                // i
                        &py_obj))                 // O, class
      return NULL;
  
  /******************************
  * Construct input numpy array *
  ******************************/
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
  /***** end of numpy array construction *****/

  /*********************************
  * setup python callback function *
  *********************************/
  /* !!!!WARNING!!!! */
  /* instance method 'evaluate' is hard coded here */
  py_objMethod = PyObject_GetAttrString(py_obj, "evaluate");
  if (!PyCallable_Check(py_objMethod)){
    PyErr_SetString(PyExc_TypeError, 
                    "parameter must be callable");
    return NULL;
  }

  Py_XINCREF(py_objMethod);   /* Add a reference to new callback */
  Py_XDECREF(my_callback);    /* Dispose of previous callback */
  my_callback = py_objMethod; /* Remember new callback */
  /***** end of callback setup *****/

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

  for(i=0;i<trows;i++){
    for(j=0;j<rrows;j++){
      getMatrixElement(i,j,rrows,rcolumns,tdata,rdata,
                       vec_dim,matrix,my_callback);
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
  free(xi);
  free(xj);
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
