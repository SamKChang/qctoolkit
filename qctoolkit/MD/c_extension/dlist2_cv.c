#include <Python.h>
#include <numpy/arrayobject.h>
#include <omp.h>
#define PI 3.14159265358979323846

//////////////////////////////////////////////
// redirect PyArray data contant to C-array //
//////////////////////////////////////////////
double *pyvector_to_Carrayptrs(PyArrayObject *arrayin){
  int n=arrayin->dimensions[0];
  /* pointer to arrayin data as double */
  return (double *) arrayin->data;
}

/*  python interface */
static PyObject* dlist_2_cv(PyObject* self, PyObject* args){

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
  double *data_frac;
  double *data;
  double cell[9];
  double cell_diag;
  double diag_comp;
  double Ri[3];
  double V_cell;

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
  data_frac = (double*) malloc(nt * N * 3 * sizeof(double));
  data = (double*) malloc(nt * N * 3 * sizeof(double));
  for(i=0;i<len0;i++){
    item = PySequence_Fast_GET_ITEM(idata, i);
    data_frac[i] = PyFloat_AsDouble(item);
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
  for(i=0;i<9;i++){
    item = PySequence_Fast_GET_ITEM(cell_py, i);
    cell[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(cell_py);


  /*
     construct liagonal length 
     cell: x1 y1 z1
           x2 y2 z2
           x3 y3 z3
     diag: (x1+x2+x3, y1+y2+y3, z1+z2+z3)
  */

  //periodic image reconstruction is not necessary for crystals
  cell_diag = 0.0;
  for(i=0;i<3;i++){ // loop for x,y,z
    diag_comp = 0.0;
    for(j=0;j<3;j++) diag_comp += cell[i + j*3]; // loop for 1,2,3
    cell_diag += pow(diag_comp, 2);
  }
  cell_diag = pow(cell_diag, 0.5);
  /***** end of input data construction *****/

  /*********************
  * construct output g *
  *********************/
  //size = (int)(0.5*cell_diag/dr) + 1;
  size = (int)(cell_diag/dr) + 1;
  mat_dim[0] = size;
  np_g = (PyArrayObject*) PyArray_FromDims(1, mat_dim, NPY_DOUBLE);
  np_r = (PyArrayObject*) PyArray_FromDims(1, mat_dim, NPY_DOUBLE);
  g = pyvector_to_Carrayptrs(np_g);
  r = pyvector_to_Carrayptrs(np_r);

  for(i=0;i<size;i++){
    r[i] = dr * (i+0.5);
    g[i] = 0;
  }

  // construct real coordinates from fractional coordinates
  for(t=0;t<nt;t++){
    for(i=0;i<N;i++){
      I = i*3 + t*N*3;
      for(j=0;j<3;j++){
        Ri[j] = 0.0;
        for(k=0;k<3;k++){
          Ri[j] += data_frac[I+k] * cell[k*3 + j];
        }
        data[I+j] = Ri[j];
      }
    }
  }
  // calculate cell valume
  V_cell = cell[0] * (cell[4] * cell[8] - cell[7] * cell[5])
          -cell[1] * (cell[3] * cell[8] - cell[6] * cell[5])
          +cell[2] * (cell[3] * cell[7] - cell[6] * cell[4]);
  // calculate density

#pragma omp parallel private(itr,t,i,j,k,I,J,Rij_t,dij) shared(g, data)
{
  #pragma omp for schedule(dynamic)
  for(t=0;t<nt;t++){
    itr = 0;
    for(i=0;i<len1;i++){
      I = atom_list1[i]*3 + t*N*3;
      for(j=0;j<len2;j++){
        J = atom_list2[j]*3 + t*N*3;
        Rij_t = 0.0;
        for(k=0;k<3;k++){
          dij = data[I+k] - data[J+k];
          //dij = dij - cell[k]*round(dij/cell[k]);
          Rij_t += pow(dij, 2);
        }
        Rij_t = sqrt(Rij_t);
        //printf("I=%d J=%d i=%d j=%d, Rij=%f\n", I, J, i, j, Rij_t);
        //if((Rij_t < 0.5*cell_diag)&&(Rij_t > 0.00001)){
        if(Rij_t > 0.00001){
          g[(int)(Rij_t/dr)] += 1.0;
        }
        itr++;
      }
    }
  }
}
  //// devide by radial shell volumn
  //for(i=0;i<size;i++){
  //  V = 4*PI*(pow((i+1)*dr,3) - pow(i*dr,3))/3.0;
  //  //g[i] /= V*nt*rho*(len1 + len2);
  //  g[i+1] /= V*nt;
  //}
  return Py_BuildValue("OO", np_r, np_g);
}

/*  define functions in module */
static PyMethodDef DList_2_cv[] ={
  {"dlist_2_cv", dlist_2_cv, METH_VARARGS,
   "loop through all time frame to calculate g of r, using fractional coordinate"},
  {NULL, NULL, 0, NULL}
};

/* module initialization */
PyMODINIT_FUNC initdlist_2_cv(void){
  (void) Py_InitModule("dlist_2_cv", DList_2_cv);
  /* IMPORTANT: this must be called */
  import_array();
}
