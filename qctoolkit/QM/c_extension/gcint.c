#include <Python.h>
#include <numpy/arrayobject.h>
#include <omp.h>

/************************************************
*  main function for Gaussian-Coulomb integral  *
************************************************/
void gcMatrix(double *data,   //output result
              double *R,      //all the rest are input
              double *Z,
              double *center,
              double *exp,
              double *cef,
              int *lm_xyz,
              int *ng,
              int Nao,
              int N,
              int aoi,
              int aoj
             ){

  /* input related variables */
  int i, j, k, i0 = 0, j0 = 0;
  double ci[3], cj[3];
  int lmi[3], lmj[3];
  int ngi = ng[aoi], ngj = ng[aoj];
  double *expi = (double*) malloc(ngi * sizeof(double));
  double *expj = (double*) malloc(ngj * sizeof(double));
  double *cefi = (double*) malloc(ngi * sizeof(double));
  double *cefj = (double*) malloc(ngj * sizeof(double));

  /* gaussian integral variables */
  double cij[3], P[3], PA[3], PB[3], PC[3], PC2;
  int lmij[3];
  int t, u, v;
  double p, mu;

	for(i=0;i<aoi;i++) i0 += ng[i];
  for(j=0;j<aoj;j++) j0 += ng[j];
  for(i=0;i<ngi;i++){
    expi[i] = exp[i0+i];
    cefi[i] = cef[i0+i];
  }
  for(j=0;j<ngj;j++){
    expj[j] = exp[j0+j];
    cefj[j] = cef[j0+j];
  }
  for(i=0;i<3;i++){
    ci[i] = center[i + 3*aoi];
    lmi[i] = lm_xyz[i + 3*aoi];
  }
  for(j=0;j<3;j++){
    cj[j] = center[j + 3*aoj];
    lmj[j] = lm_xyz[j + 3*aoj];
  }
  for(i=0;i<3;i++){
    cij[i] = ci[i] - cj[i];
    lmij[i] = lmi[i] + lmj[i];
  }

  /* normalization constant */
  /* Boys function */
  /* Hermite-Coulomb coefficient */

  data[aoj+aoi*Nao] = 0.1;
  free(expi);
  free(expj);
  free(cefi);
  free(cefj);
}

/*  python interface */
static PyObject* gcint(PyObject* self, PyObject* args){

  /*** input variables ***/
  /* dictionary */
  PyObject *in_dict;
  PyObject *dictList;
  PyObject *dictList_fast;
  PyObject *pyStr;
  PyObject *item;
  /* list */
  PyObject *in_list;
  PyObject *list_fast;
  /* numpy array */
  PyArrayObject *in_array1, *in_array2, *in_array3;
  NpyIter *in_iter;
  NpyIter_IterNextFunc *in_iternext;
  /* C data type for input */
  int Ngo, Nao;
  int N, itr;
  double *Z;
  double *R;
  /* basis function variables */
  double **in_ptr; // address to numpy pointer
  double *center;  // gaussian center
  double *cef;     // contraction coefficients
  int *ng;         // number of gaussians per AO
  double *exp;     // gaussian exponents
  int *lm_xyz;     // angular momentum

  /* python output variables */
  PyObject *py_out;
  double *data;
  int i, j;
  int mat_dim[2];

  /*  parse numpy array and two integers as argument */
  if (!PyArg_ParseTuple(args, "OO!O!O!O", 
                        &in_dict, 
                        &PyArray_Type, &in_array1,
                        &PyArray_Type, &in_array2,
                        &PyArray_Type, &in_array3,
                        &in_list
                       )) return NULL;
  if(in_array1 == NULL) return NULL;

  /***********************
  * Construct input data *
  ***********************/
  /*** access dict_list data ***/
  /* exponents */
  pyStr = PyString_FromString("exponents");
  dictList = PyDict_GetItem(in_dict, pyStr);
  dictList_fast = PySequence_Fast(
                    dictList, "expected a sequence");
  Ngo = PySequence_Size(dictList);
  exp = (double*) malloc(Ngo * sizeof(double));
  for(i=0;i<Ngo;i++){
    item = PySequence_Fast_GET_ITEM(dictList_fast, i);
    exp[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(dictList_fast);
  /* coefficients */
  pyStr = PyString_FromString("coefficients");
  dictList = PyDict_GetItem(in_dict, pyStr);
  dictList_fast = PySequence_Fast(
                    dictList, "expected a sequence");
  cef = (double*) malloc(Ngo * sizeof(double));
  for(i=0;i<Ngo;i++){
    item = PySequence_Fast_GET_ITEM(dictList_fast, i);
    cef[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(dictList_fast);
  /* number of gaussians per shell */
  pyStr = PyString_FromString("n_gaussians");
  dictList = PyDict_GetItem(in_dict, pyStr);
  dictList_fast = PySequence_Fast(
                    dictList, "expected a sequence");
  Nao = PySequence_Size(dictList);
  ng = (int*) malloc(Nao * sizeof(int));
  for(i=0;i<Nao;i++){
    item = PySequence_Fast_GET_ITEM(dictList_fast, i);
    ng[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(dictList_fast);
  /*** end of dict_list data ***/

  /*** access python list data ***/
  list_fast = PySequence_Fast(in_list, "expected a sequence");
  N = PySequence_Size(in_list);
  Z = (double*) malloc(N * sizeof(double));
  for(i=0;i<N;i++){
    item = PySequence_Fast_GET_ITEM(list_fast, i);
    Z[i] = PyFloat_AsDouble(item);
  }
  Py_DECREF(list_fast);
  /*** end of list input data ***/

  /*** create the iterators, necessary to access numpy array ***/
  /* center */
  in_iter = NpyIter_New(in_array1, NPY_ITER_READONLY,
                         NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  if (in_iter == NULL) goto fail;
  in_iternext = NpyIter_GetIterNext(in_iter, NULL);
  if (in_iternext == NULL){
    NpyIter_Deallocate(in_iter);
    goto fail;
  }
  /* interator pointer to actual numpy data */
  in_ptr = (double **) NpyIter_GetDataPtrArray(in_iter);
  center = (double*) malloc(3 * Nao * sizeof(double));
  itr=0;
  do {
    center[itr++] = **in_ptr;
  } while(in_iternext(in_iter));
  NpyIter_Deallocate(in_iter);

  /* lm_xyz */
  in_iter = NpyIter_New(in_array2, NPY_ITER_READONLY,
                         NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  if (in_iter == NULL) goto fail;
  in_iternext = NpyIter_GetIterNext(in_iter, NULL);
  if (in_iternext == NULL){
    NpyIter_Deallocate(in_iter);
    goto fail;
  }
  /* interator pointer to actual numpy data */
  int **in_ptr2 = (int **) NpyIter_GetDataPtrArray(in_iter);
  lm_xyz = (int*) malloc(3 * Nao * sizeof(int));
  itr=0;
  do {
    lm_xyz[itr++] = **in_ptr2;
  } while(in_iternext(in_iter));
  NpyIter_Deallocate(in_iter);

  /* R */
  in_iter = NpyIter_New(in_array3, NPY_ITER_READONLY,
                         NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  if (in_iter == NULL) goto fail;
  in_iternext = NpyIter_GetIterNext(in_iter, NULL);
  if (in_iternext == NULL){
    NpyIter_Deallocate(in_iter);
    goto fail;
  }
  /* interator pointer to actual numpy data */
  in_ptr = (double **) NpyIter_GetDataPtrArray(in_iter);
  R = (double*) malloc(3 * N * sizeof(double));
  itr=0;
  do {
    R[itr++] = **in_ptr;
  } while(in_iternext(in_iter));
  NpyIter_Deallocate(in_iter);
  /*** end of numpy input data ***/

  /***** end of input data construction *****/

  /*******************
  * construct output *
  *******************/
  data = (double*) malloc(Nao * Nao * sizeof(double));

  for(i=0;i<Nao;i++){
    for(j=i;j<Nao;j++){
      gcMatrix(data, R, Z, 
               center, exp, cef, ng, lm_xyz, 
               Nao, N, i, j);
      if(j!=i) data[i+j*Nao] = data[j+i*Nao];
    }
  }

  mat_dim[0] = Nao;
  mat_dim[1] = Nao;
  py_out = PyArray_SimpleNewFromData(2, mat_dim, NPY_DOUBLE, data);

  /*********************************
  * clean up and return the result *
  *********************************/
  free(exp);
  free(cef);
  free(ng);
  free(center);
  free(R);
  free(Z);

  Py_INCREF(py_out);
  return Py_BuildValue("O", py_out);

  /*  in case bad things happen */
  fail:
    Py_XDECREF(py_out);
    return NULL;
}

/*  define functions in module */
static PyMethodDef GCInt[] ={
  {"gcint", gcint, METH_VARARGS,
      "analytic Gaussian-Coulomb integral"},
  {NULL, NULL, 0, NULL}
};

/* module initialization */
PyMODINIT_FUNC initgcint(void){
  (void) Py_InitModule("gcint", GCInt);
  /* IMPORTANT: this must be called */
  import_array();
}
