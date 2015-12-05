#include <Python.h>
#include <numpy/arrayobject.h>
#include <omp.h>

/*  python interface */
static PyObject* vacf(PyObject* self, PyObject* args){

  /* input variables */
  PyObject *in_array;
  PyObject *idata;
  PyObject *list;
  PyObject *list_inp;
  PyObject *item;
  int *atom_list;
  int nt, N, len0, len1;
  double *data;

  /* python output variables */
  PyObject *np_vax, *np_vay, *np_vaz, *np_va;
  double *vax, *vay, *vaz, *va;
  int *count;
  int i, I_0, I_t, t, dt;
  int mat_dim[1];
  int size;

  /*  parse numpy array and two integers as argument */
  //                      &PyArray_Type, &in_array, // O!, NpArray
  if (!PyArg_ParseTuple(args, "OiiO", 
                        &in_array, // O, PyObject
                        &nt,       // i, integer
                        &N,        // i, integer
                        &list_inp  // O, PyObject
                       )) return NULL;
  if(in_array == NULL) return NULL;

  /***********************
  * Construct input data *
  ***********************/
  /* access python list data */
  idata = PySequence_Fast(in_array, "expected a sequence");
  len0 = PySequence_Size(in_array);
  //nt = len0/3/N;
  data = (double*) malloc(nt * N * 3 * sizeof(double));
  for(i=0;i<len0;i++){
    item = PySequence_Fast_GET_ITEM(idata, i);
    data[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(idata);

  list = PySequence_Fast(list_inp, "expected a sequence");
  len1 = PySequence_Size(list_inp);
  atom_list = (int*) malloc(len1 * sizeof(int));
  for(i=0;i<len1;i++){
    item = PySequence_Fast_GET_ITEM(list, i);
    atom_list[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(list);
  /***** end of input data construction *****/

  /*******************
  * construct output *
  *******************/
  vax = (double *) malloc(nt * sizeof(double));
  vay = (double *) malloc(nt * sizeof(double));
  vaz = (double *) malloc(nt * sizeof(double));
  va = (double *) malloc(nt * sizeof(double));
  count = (int *) malloc(nt * sizeof(int));
  for(i=0;i<nt;i++){
    vax[i] = 0.0;
    vay[i] = 0.0;
    vaz[i] = 0.0;
    va[i] = 0.0;
    count[i] = 0;
  }

#pragma omp parallel private(t, dt, i, I_0, I_t) \
shared(vax, vay, vaz, va)
{
  #pragma omp for schedule(dynamic)
  for(t=0;t<nt-1;t++){
    for(dt=0;(t+dt)<nt;dt++){
      count[dt]++;
      for(i=0;i<len1;i++){
        I_0 = atom_list[i]*3 + t*N*3;
        I_t = atom_list[i]*3 + (dt+t)*N*3;
        vax[dt] += data[I_0] * data[I_t];
        vay[dt] += data[I_0+1] * data[I_t+1];
        vaz[dt] += data[I_0+2] * data[I_t+2];
        va[dt] += vax[dt] + vay[dt] + vaz[dt];
      }
    }
  }
}
  for (t=0;t<nt;t++) {
    vax[t] /= count[t]?(len1*count[t]):1;
    vay[t] /= count[t]?(len1*count[t]):1;
    vaz[t] /= count[t]?(len1*count[t]):1;
    va[t] /= count[t]?(len1*count[t]):1;
  }

  mat_dim[0] = nt;
  np_vax = PyArray_SimpleNewFromData(1, mat_dim, NPY_DOUBLE, vax);
  np_vay = PyArray_SimpleNewFromData(1, mat_dim, NPY_DOUBLE, vay);
  np_vaz = PyArray_SimpleNewFromData(1, mat_dim, NPY_DOUBLE, vaz);
  np_va = PyArray_SimpleNewFromData(1, mat_dim, NPY_DOUBLE, va);
  /***** end of output g construction *****/

  /*********************************
  * clean up and return the result *
  *********************************/
  Py_INCREF(np_vax);
  Py_INCREF(np_vay);
  Py_INCREF(np_vaz);
  Py_INCREF(np_va);
  return Py_BuildValue("OOOO", np_va, np_vax, np_vay, np_vaz);

  /*  in case bad things happen */
  fail:
    Py_XDECREF(np_vax);
    Py_XDECREF(np_vay);
    Py_XDECREF(np_vaz);
    Py_XDECREF(np_va);
    return NULL;

  free(data);
}

/*  define functions in module */
static PyMethodDef Vacf[] ={
  {"vacf", vacf, METH_VARARGS,
      "velocity autocorrelation function"},
  {NULL, NULL, 0, NULL}
};

/* module initialization */
PyMODINIT_FUNC initvacf(void){
  (void) Py_InitModule("vacf", Vacf);
  /* IMPORTANT: this must be called */
  import_array();
}
