/*  Example of wrapping the cos function from math.h using the Numpy-C-API. */

#include <iostream>
#include <Python.h>
#include <numpy/arrayobject.h>
#include "kernels.h"
#include <string.h>
#define NAME_SIZE 30

///* python callable pointer */
//static PyObject *my_callback = NULL;

/* extract row vector with index='index' in 'data' matrix */
void getVector(double **vector, double *data, int index, int max){
  int itr;
  delete(*vector); /* free previously allocated memory */
  *vector = new double[max];
  for(itr=0;itr<max;itr++){
    (*vector)[itr] = data[index * max + itr];
  }
}

/*  python interface */
static PyObject* kernel_matrix(PyObject* self, 
                               PyObject* args,
                               PyObject* kwargs){
  int itr = 0;

  /* numpy input variables */
  PyArrayObject *in_array;
  NpyIter *in_iter;
  NpyIter_IterNextFunc *in_iternext;
  double *data;
  int rows, columns;
  char kernel[NAME_SIZE];

  /* optional kernel args and kwargs */
  static char *kwlist[] = {"to", "be", "added"};
  double klarg[5];

//  /* python callback variables */
//  PyObject *result = NULL;
//  PyObject *py_obj;
//  PyObject *py_objMethod;
//  PyObject *arglist;
//  PyObject *np_xi, *np_xj;

  double *xi = new double[1];
  double *xj = new double[1];
  int vec_dim[1];

  /* python output variables */
  PyObject *np_matrix;
//  PyObject *element;
  double *matrix;
  int mat_dim[2];
  int i, j;

  /*  parse numpy array and two integers as argument */
//  if (!PyArg_ParseTuple(args, "O!iiOs", 
//                        &PyArray_Type, &in_array, // O!, NpArray
//                        &rows,                    // i
//                        &columns,                 // i
//                        kernel                    // s, kernelStr
////                        &py_obj))                 // O, class
  if (!PyArg_ParseTupleAndKeywords(
                        args, 
                        kwargs,
                        "O!iis|d", 
                        kwlist,
                        &PyArray_Type, &in_array, // O!, NpArray
                        &rows,                    // i
                        &columns,                 // i
                        kernel,                   // s, kernelStr
            /* optional kernel argument */
                        &klarg[0]                 // d, klarg[0]
                       ))
      return NULL;
  
  if(strcmp(kernel,"Gaussian")==0){
    Gaussian ker(klarg[0]);
  }else{
    std::cout << kernel << " not implemented" << std::endl;
    return NULL;
  }



  /******************************
  * Construct input numpy array *
  ******************************/
  data = new double[rows * columns];
  /*  create the iterators, necessary to access numpy array */
  in_iter = NpyIter_New(in_array, NPY_ITER_READONLY, 
                        NPY_KEEPORDER, NPY_NO_CASTING, NULL);
//  /* check input status */
//  if (in_iter == NULL)
//      goto fail;
//  in_iternext = NpyIter_GetIterNext(in_iter, NULL);
//  if (in_iternext == NULL) {
//      NpyIter_Deallocate(in_iter);
//      goto fail;
//  }
  /* interator pointer to actual numpy data */
  double **in_dataptr = (double **) 
                           NpyIter_GetDataPtrArray(in_iter);
  /*  iterate over the arrays */
  do {
      data[itr++] = **in_dataptr;
  } while(in_iternext(in_iter));
  /***** end of numpy array construction *****/

//  /*********************************
//  * setup python callback function *
//  *********************************/
//  /* !!!!WARNING!!!! */
//  /* instance method 'evaluate' is hard coded here */
//  py_objMethod = PyObject_GetAttrString(py_obj, "evaluate");
//  if (!PyCallable_Check(py_objMethod)){
//    PyErr_SetString(PyExc_TypeError, 
//                    "parameter must be callable");
//    return NULL;
//  }
//
//  Py_XINCREF(py_objMethod);   /* Add a reference to new callback */
//  Py_XDECREF(my_callback);    /* Dispose of previous callback */
//  my_callback = py_objMethod; /* Remember new callback */
//  /***** end of callback setup *****/

  /**************************
  * construct output matrix *
  **************************/
  matrix = new double[rows * rows];
  mat_dim[0] = rows;
  mat_dim[1] = rows;
  vec_dim[0] = columns;

  for(i=0;i<rows;i++){
    for(j=i;j<rows;j++){
      /* construct callback input numpy arrays*/
      getVector(&xi, data, i, columns);
      getVector(&xj, data, j, columns);
      
//      np_xi = PyArray_SimpleNewFromData(1, vec_dim, 
//                                        NPY_DOUBLE,xi);
//      np_xj = PyArray_SimpleNewFromData(1, vec_dim, 
//                                        NPY_DOUBLE,xj);
//      /* construct callback argument list */
//      arglist = Py_BuildValue("(OO)", np_xi, np_xj);
//      /* simple python call with return PyObject */
//      element = PyObject_CallObject(my_callback, arglist);
//      /* extract 'double' data type from PyObject */
//      matrix[i+rows*j] = PyFloat_AsDouble(element);
      if(i!=j) matrix[j+rows*i] = matrix[i+rows*j];
//      Py_DECREF(arglist);
      //Py_INCREF(element); // ????
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
//  fail:
//      Py_XDECREF(np_matrix);
//      return NULL;

  delete(data);
  delete(xi);
  delete(xj);
  delete(matrix);
}

/*  define functions in module */
static PyMethodDef KernelMatrix[] ={
  {"kernel_matrix", kernel_matrix, METH_VARARGS,
      "evaluate the cosine on a numpy array"},
  {NULL, NULL, 0, NULL}
};

/* module initialization */
PyMODINIT_FUNC initkernel_matrix(void){
  extern "C" (void) Py_InitModule("kernel_matrix", KernelMatrix);
  /* IMPORTANT: this must be called */
  import_array();
}
