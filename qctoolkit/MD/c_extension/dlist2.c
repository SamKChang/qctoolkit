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
  double dr, rho;
  double *data;
  double cell[3];
  double cell_min;

  /* python output variables */
  PyObject *np_g;
  PyObject *np_r;
  double *g, *r;
  int mat_dim[1];
  int i, j, k, t, itr;
  int I, J;
  int size;
  double Rij_t, dij, V;

  /*  parse numpy array and two integers as argument */
  //                      &PyArray_Type, &in_array, // O!, NpArray
  if (!PyArg_ParseTuple(args, "OiiOOOd", 
                        &in_array, // O, PyObject
                        &nt,       // i, integer
                        &N,        // i, integer
                        &l1_inp,   // O, PyObject
                        &l2_inp,   // O, PyObject
                        &cell_inp, // O, PyObject
                        &dr        // f, double
                       )) return NULL;
  if(in_array == NULL) return NULL;

  /***********************
  * Construct input data *
  ***********************/
  /* access python list data */
  idata = PySequence_Fast(in_array, "expected a sequence");
  len0 = PySequence_Size(in_array);
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
  cell_min = cell[0];
  rho = len1 + len2;
  for(i=0;i<3;i++){
    if(cell[i]<cell_min) cell_min = cell[i];
    rho /= cell[i];
  }
  /***** end of input data construction *****/

  /*********************
  * construct output g *
  *********************/
  size = (int)(0.5*cell_min/dr) + 1;
  g = (double *) malloc(size * sizeof(double));
  r = (double *) malloc(size * sizeof(double));
  for(i=0;i<size;i++){
    r[i] = dr * (i+0.5);
    g[i] = 0;
  }

#pragma omp parallel private(itr,i,j,t,I,J,Rij_t,dij) shared(g)
{
  #pragma omp for schedule(dynamic)
  for(t=0;t<nt;t++){
    itr = 0;
    for(i=0;i<len1-1;i++){
      I = atom_list1[i]*3 + t*N*3;
      for(j=i+1;j<len2;j++){
        J = atom_list2[j]*3 + t*N*3;
        Rij_t = 0.0;
        for(k=0;k<3;k++){
          dij = data[I+k] - data[J+k];
          dij = dij - cell[k]*round(dij/cell[k]);
          Rij_t += pow(dij, 2);
        }
        Rij_t = sqrt(Rij_t);
        if(Rij_t < 0.5*cell_min){
          g[(int)(Rij_t/dr)] += 2.0;
        }
        itr++;
      }
    }
  }
}
  for(i=0;i<size;i++){
    V = 4*PI*(pow((i+1)*dr,3) - pow(i*dr,3))/3.0;
    g[i] /= V*nt*rho*len1;
  }

  mat_dim[0] = size;
  np_g = PyArray_SimpleNewFromData(1, mat_dim, NPY_DOUBLE, g);
  np_r = PyArray_SimpleNewFromData(1, mat_dim, NPY_DOUBLE, r);
  /***** end of output g construction *****/

  /*********************************
  * clean up and return the result *
  *********************************/
  Py_INCREF(np_g);
  Py_INCREF(np_r);
  return Py_BuildValue("OO", np_r, np_g);

  /*  in case bad things happen */
  fail:
    Py_XDECREF(np_g);
    Py_XDECREF(np_r);
    return NULL;

  free(data);
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
