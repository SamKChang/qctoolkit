/**************************************************************
*  C interface to take a Numpy array of descriptor vectors    *
*  as input. Callback python kernel class methond 'evaluate'  *
*  to obtain kernel matrix                                    *
*  adapted from numpy interface example at:                   *
*  https://scipy-lectures.github.io/advanced/                 *
*    interfacing_with_c/interfacing_with_c.html               *
**************************************************************/
#include <Python.h>
#include <numpy/arrayobject.h>
#include <omp.h>
#include <pthread.h>

/* python callable pointer */
static PyObject *my_callback = NULL;

/* extract row vector with index='index' in 'data' matrix */
void getVector(double *vector, double *data, int index, int max){
  int itr;
  //vector = (double*) malloc(max * sizeof(double));
  for(itr=0;itr<max;itr++){
    vector[itr] = data[index * max + itr];
  }
}

/*  python interface */
static PyObject* kernel_matrix(PyObject* self, PyObject* args){

  int itr = 0;

  /* numpy input variables */
  PyArrayObject *in_array;
  NpyIter *in_iter;
  NpyIter_IterNextFunc *in_iternext;
  double *data;
  int rows, columns;

  /* python callback variables */
  PyObject *result = NULL;
  PyObject *py_obj;
  PyObject *py_objMethod;
  PyObject *arglist;
  PyObject *np_xi, *np_xj;
//  double *xi = (double *) malloc(sizeof(double));
//  double *xj = (double *) malloc(sizeof(double));
  double *xi, *xj;
  int vec_dim[1];

  /* python output variables */
  PyObject *np_matrix;
  PyObject *element;
  double *matrix;
  int mat_dim[2];
  int i, j, s;

  /*  parse numpy array and two integers as argument */
  if (!PyArg_ParseTuple(args, "O!iiO", 
                        &PyArray_Type, &in_array, // O!, NpArray
                        &rows,                    // i
                        &columns,                 // i
                        &py_obj))                 // O, class
      return NULL;
  
  /******************************
  * Construct input numpy array *
  ******************************/
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
  matrix = (double *) malloc(rows * rows * sizeof(double));
  for(i=0;i<rows;i++){
    for(j=i;j<rows;j++){
      matrix[j+i*rows] = 0;
      if(i!=j) matrix[i+j*rows] = 0;
    }
  }
  mat_dim[0] = rows;
  mat_dim[1] = rows;
  vec_dim[0] = columns;

  /* openmp parallelization */
  /* parallel test */
  printf("yo %d\n", rows);

Py_BEGIN_ALLOW_THREADS
  #pragma omp parallel for default(shared) \
   private(s, i, j ,xi, xj, np_xi, np_xj, arglist, element)
//   shared(data, matrix)
//    for(s=0;s<rows*rows;s++){
//      printf("%d\n",s);
//    }
   //private(i, j, xi, xj, np_xi, np_xj, arglist, element) 
     

    for(s=0;s<rows*rows;s++){
//    for(i=0;i<rows;i++){
//      for(j=i;j<rows;j++){
        i = s/rows;
        j = s%rows;
        if(j>=i){
          xi = (double *) malloc(columns * sizeof(double));
          xj = (double *) malloc(columns * sizeof(double));
          /* construct callback input numpy arrays*/
          getVector(xi, data, i, columns);
          getVector(xj, data, j, columns);
          
          np_xi = PyArray_SimpleNewFromData(1, vec_dim, 
                                            NPY_DOUBLE,xi);
          //printf("%d %d %d\n", s, i, j);
          np_xj = PyArray_SimpleNewFromData(1, vec_dim, 
                                            NPY_DOUBLE,xj);
          /* construct callback argument list */
          arglist = Py_BuildValue("(OO)", np_xi, np_xj);
          /* simple python call with return PyObject */
          element = PyObject_CallObject(my_callback, arglist);
          printf("%d %d %d\n", s, i, j);
          /* extract 'double' data type from PyObject */
          matrix[i+rows*j] = PyFloat_AsDouble(element);
          if(i!=j) matrix[j+rows*i] = matrix[i+rows*j];
          //Py_INCREF(element); // ????
          printf("%d %d %d, %lf\n", s, i, j, matrix[i+rows*j]);
          free(xi);
          free(xj);
        }
//      }
    }
  /* end of openmp parallelization */
Py_END_ALLOW_THREADS
  Py_DECREF(arglist);
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
//  free(xi);
//  free(xj);
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
