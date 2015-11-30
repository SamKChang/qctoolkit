#include <Python.h>
#include <numpy/arrayobject.h>
#include <omp.h>
#include "gaussian.h"

/************************************************
*  main function for Gaussian-Coulomb integral  *
************************************************/
//gcMatrix(data, R, Z, center, exp, cef, ng, 
//         lm_xyz, Nao, N, i, j);
void gcMatrix(double *data,   //output result
              double *R,      //all the rest are input
              double *Z,
              double *center,
              double *exp,
              double *cef,
              int *ng,
              int *lm_xyz,
              int Nao,
              int N,
              int aoi,
              int aoj
             ){

  /* input related variables */
  int i, j, k, i0 = 0, j0 = 0, I;
  double ci[3], cj[3];
  int lmi[3], lmj[3];
  int ngi = ng[aoi], ngj = ng[aoj];
  double *expi = (double*) malloc(ngi * sizeof(double));
  double *expj = (double*) malloc(ngj * sizeof(double));
  double *cefi = (double*) malloc(ngi * sizeof(double));
  double *cefj = (double*) malloc(ngj * sizeof(double));

  /* gaussian integral variables */
  double cI[3], cij[3], P[3], Pi[3], Pj[3], PI[3], PI2;
  int lmij[3];
  int t, u, v;
  double p, p2, mu, ZI, x, HC_cef, Ni, Nj, Nij;
  double Hx, Hy, Hz, num, element_ij = 0;

  for(i=0;i<aoi;i++) i0 += ng[i];
  for(j=0;j<aoj;j++) j0 += ng[j];
  for(k=0;k<3;k++){
    lmi[k] = lm_xyz[k + 3*aoi];
    lmj[k] = lm_xyz[k + 3*aoj];
    lmij[k] = lmi[k] + lmj[k];
    ci[k] = center[k + 3*aoi];
    cj[k] = center[k + 3*aoj];
    cij[k] = ci[k] - cj[k];
  }
  for(i=0;i<ngi;i++){
    expi[i] = exp[i0+i];
    cefi[i] = cef[i0+i];
    for(j=0;j<ngj;j++){
      expj[j] = exp[j0+j];
      cefj[j] = cef[j0+j];
      p = expi[i] + expj[j];
      p2 = 2*p;
      mu = expi[i] * expj[j] / p;
      Ni = cefi[i] * Norm(expi[i], lmi);
      Nj = cefj[j] * Norm(expj[j], lmj);
      Nij = Ni * Nj;
      for(k=0;k<3;k++){
        P[k] = (expi[i]*ci[k] + expj[j]*cj[k])/p;
        Pi[k] = P[k] - ci[k];
        Pj[k] = P[k] - cj[k];
      }
      for(I=0;I<N;I++){
        PI2 = 0;
        for(k=0;k<3;k++){
          cI[k] = R[k+3*I];
          PI[k] = P[k] - cI[k];
          PI2 += PI[k] * PI[k];
        }
        ZI = Z[I];
        x = p*PI2;

        for(t=0;t<lmij[0]+1;t++){
          Hx = Hermite(lmi[0], lmj[0], t, p2, mu, 
                       cij[0], Pi[0], Pj[0]);
          for(u=0;u<lmij[1]+1;u++){
            Hy = Hermite(lmi[1], lmj[1], u, p2, mu, 
                         cij[1], Pi[1], Pj[1]);
            for(v=0;v<lmij[2]+1;v++){
              Hz = Hermite(lmi[2], lmj[2], v, p2, mu, 
                           cij[2], Pi[2], Pj[2]);
              HC_cef = HCcef(t, u, v, 0, PI, p2, x);
              num = ZI*(Nij*Hx*Hy*Hz*HC_cef);
              element_ij -= num*(2*M_PI/p);
/*
              if(num*(2*M_PI/p)*num > 1E-7){
                if(aoi==1 && aoj==1){
                  printf("% 10.6E a=%7.3E b=%7.3E: ", 
                         -num*(2*M_PI/p), expi[i], expj[j]);
                  printf("lm_a=[%d %d %d] lm_b=[%d %d %d]\n  ",
                         lmi[0],lmi[1],lmi[2],
                         lmj[0],lmj[1],lmj[2]);
                  printf("Nij:%6.2E Ni:%6.2E Nj:%6.2E ", 
                         Nij, Ni, Nj);
                  printf("Hx=%6.2E Hy=%6.2E Hz=%6.2E HC=%6.2E\n  ",
                         Hx, Hy, Hz, HC_cef);
                  printf("HCcef(%d, %d, %d, 0, PI, %3.1f %3.1f)=%f\n",
                         t, u, v, p2, x, HC_cef);
                  printf("   PI: %6.4f %6.4f %6.4f\n", 
                         PI[0], PI[1], PI[2]);
                }
              }
*/
            }
          }
        }
      }
    }
  }

  data[aoj+aoi*Nao] = element_ij;
  free(expi);
  free(expj);
  free(cefi);
  free(cefj);
}

/*********************
*  python interface  *
*********************/
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
  double *data, *overlap;
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
  overlap = (double*) malloc(Nao * Nao * sizeof(double));

  /* renormalization */
  renormalize(center, exp, cef, ng, lm_xyz, Nao);

  /* orthogonalize */
  // no effect at the moment, for debug purpose
  //orthogonalize(overlap, center, exp, cef, ng, lm_xyz, Nao);

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
  free(overlap);

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
