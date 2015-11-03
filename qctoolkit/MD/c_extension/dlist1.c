#include <Python.h>
#include <numpy/arrayobject.h>
#include <omp.h>

/*  python interface */
static PyObject* dlist_1(PyObject* self, PyObject* args){

  /* input variables */
  PyObject *in_array;
  PyObject *idata;
  PyObject *l1_inp;
  PyObject *l2_inp;
  PyObject *l1;
  PyObject *l2;
  PyObject *item;
  int *atom_list1;
  int *atom_list2;
  int nt, N, len0, len1, len2;
  double *data;

  /* python output variables */
  PyObject *np_matrix;
  double *matrix;
  int mat_dim[1];
  int i, j, k, t;
  int I, J;
  int itr = 0;
  double Rij_t;

  /*  parse numpy array and two integers as argument */
  //                      &PyArray_Type, &in_array, // O!, NpArray
  if (!PyArg_ParseTuple(args, "OiiOO", 
                        &in_array, // O, PyObject
                        &nt,       // i, integer
                        &N,        // i, integer
                        &l1_inp,   // O, PyObject
                        &l2_inp    // O, PyObject
                       )) return NULL;
  if(in_array == NULL) return NULL;

  /***********************
  * Construct input data *
  ***********************/
  /* access python list data */
  idata = PySequence_Fast(in_array, "expected a sequence");
  len0 = PySequence_Size(in_array);
  nt = len0/3/N;
  data = (double*) malloc(nt * N * 3 * sizeof(double));
  for(i=0;i<len0;i++){
    item = PySequence_Fast_GET_ITEM(idata, i);
    data[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(idata);

  l1 = PySequence_Fast(l1_inp, "expected a sequence");
  len1 = PySequence_Size(l1_inp);
  atom_list1 = (int*) malloc(len1 * sizeof(int));
  for(i=0;i<len1;i++){
    item = PySequence_Fast_GET_ITEM(l1, i);
    atom_list1[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(l1);

  l2 = PySequence_Fast(l2_inp, "expected a sequence");
  len2 = PySequence_Size(l2_inp);
  atom_list2 = (int*) malloc(len2 * sizeof(int));
  for(i=0;i<len2;i++){
    item = PySequence_Fast_GET_ITEM(l2, i);
    atom_list2[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(l2);
  /***** end of input data construction *****/

  /**************************
  * construct output matrix *
  **************************/
  matrix = (double *) malloc(len1 * len2 * nt * sizeof(double));

#pragma omp parallel private(i,j,t,I,J,Rij_t) shared(matrix)
{
  #pragma omp for
  for(t=0;t<nt;t++){
    for(i=0;i<len1-1;i++){
      I = i + t*N*3;
      for(j=i+1;j<len2;j++){
        J = j + t*N*3;
        Rij_t = 0.0;
        for(k=0;k<3;k++){
          Rij_t += pow((data[I+k] - data[J+k]), 2);
        }
        matrix[j + i*len2 + t*len1*len2] = pow(Rij_t, 0.5);
      }
    }
  }
}

  mat_dim[0] = len1 * len2;
  np_matrix = PyArray_SimpleNewFromData(1, mat_dim, 
                                        NPY_DOUBLE, matrix);
  /***** end of output matrix construction *****/

  /*********************************
  * clean up and return the result *
  *********************************/
  Py_INCREF(np_matrix);
  return np_matrix;

  /*  in case bad things happen */
  fail:
      Py_XDECREF(np_matrix);
      return NULL;

  free(data);
  free(matrix);
}

/*  define functions in module */
static PyMethodDef DList_1[] ={
  {"dlist_1", dlist_1, METH_VARARGS,
      "loop through all time frame to calculate g or r"},
  {NULL, NULL, 0, NULL}
};

/* module initialization */
PyMODINIT_FUNC initdlist_1(void){
  (void) Py_InitModule("dlist_1", DList_1);
  /* IMPORTANT: this must be called */
  import_array();
}
