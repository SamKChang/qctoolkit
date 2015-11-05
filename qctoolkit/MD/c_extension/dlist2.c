#include <Python.h>
#include <numpy/arrayobject.h>
#include <omp.h>
#define PI 3.14159265358979323846

/*  python interface */
static PyObject* dlist_2(PyObject* self, PyObject* args){

  /* input variables */
  PyObject *in_array;
  PyObject *idata;
  PyObject *l1_inp;
  PyObject *l2_inp;
  PyObject *l1;
  PyObject *l2;
  PyObject *cell_inp;
  PyObject *cell_py;
  PyObject *item;
  int *atom_list1;
  int *atom_list2;
  int nt, N, len0, len1, len2;
  double *data;
  double cell[3];

  /* python output variables */
  PyObject *np_matrix;
  double *matrix;
  int mat_dim[1];
  int i, j, k, t;
  int I, J;
  int itr = 0;
  double Rij_t, dij, surf;

  /*  parse numpy array and two integers as argument */
  //                      &PyArray_Type, &in_array, // O!, NpArray
  if (!PyArg_ParseTuple(args, "OiiOOO", 
                        &in_array, // O, PyObject
                        &nt,       // i, integer
                        &N,        // i, integer
                        &l1_inp,   // O, PyObject
                        &l2_inp,   // O, PyObject
                        &cell_inp  // O, PyObject
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

  cell_py = PySequence_Fast(cell_inp, "expected a sequence");
  for(i=0;i<3;i++){
    item = PySequence_Fast_GET_ITEM(cell_py, i);
    cell[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(cell_py);
  /***** end of input data construction *****/

  /**************************
  * construct output matrix *
  **************************/
  matrix = (double *) malloc(len1 * len2 * nt * sizeof(double));

#pragma omp parallel private(i,j,t,I,J,Rij_t) shared(matrix)
{
  #pragma omp for
  for(t=0;t<nt;t++){
    for(i=0;i<len1;i++){
      I = i + t*N*3;
      for(j=0;j<len2;j++){
        J = j + t*N*3;
        Rij_t = 0.0;
        for(k=0;k<3;k++){
          dij = data[I+k] - data[J+k];
          if(dij>cell[k]) dij = cell[k] - dij;
          Rij_t += pow(dij, 2);
        }
        surf = 4*PI*Rij_t;
        matrix[j + i*len2 + t*len1*len2] = pow(Rij_t, 0.5)/surf;
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
static PyMethodDef DList_2[] ={
  {"dlist_2", dlist_2, METH_VARARGS,
      "loop through all time frame to calculate g or r"},
  {NULL, NULL, 0, NULL}
};

/* module initialization */
PyMODINIT_FUNC initdlist_2(void){
  (void) Py_InitModule("dlist_2", DList_2);
  /* IMPORTANT: this must be called */
  import_array();
}
