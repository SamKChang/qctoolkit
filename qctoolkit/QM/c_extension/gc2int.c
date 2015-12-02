#include <Python.h>
#include <numpy/arrayobject.h>
#include <omp.h>
#include "gaussian.h"

/**********************************
*  necessary Numpy functionality  *
**********************************/
// return 1D memory array pointer to allocated numpy array 
double* pyvector_to_Carrayptrs(PyArrayObject *arrayin){
  int n = arrayin->dimensions[0];
  return (double*) arrayin->data;
}

///************************************************
//*  main function for Gaussian-Coulomb integral  *
//************************************************/
////gc2Matrix(data, center, exp, cef, ng, 
////         lm_xyz, Nao, N, i, j, k, l);
//double gc2Matrix(double *center,
//               double *exp,
//               double *cef,
//               int *ng,
//               int *lm_xyz,
//               int Nao,
//               int aoi,
//               int aoj,
//               int aok,
//               int aol
//              ){
//
//  /* input related variables */
//  int s, i, j, k, l, i0 = 0, j0 = 0, k0 = 0, l0 = 0;
//  double ci[3], cj[3], ck[3], cl[3];
//  int lmi[3], lmj[3], lmk[3], lml[3];
//  int ngi = ng[aoi], ngj = ng[aoj], ngk = ng[aok], ngl = ng[aol];
//  double *expi = (double*) malloc(ngi * sizeof(double));
//  double *expj = (double*) malloc(ngj * sizeof(double));
//  double *expk = (double*) malloc(ngk * sizeof(double));
//  double *expl = (double*) malloc(ngl * sizeof(double));
//  double *cefi = (double*) malloc(ngi * sizeof(double));
//  double *cefj = (double*) malloc(ngj * sizeof(double));
//  double *cefk = (double*) malloc(ngk * sizeof(double));
//  double *cefl = (double*) malloc(ngl * sizeof(double));
//
//  /* gaussian integral variables */
//  double cI[3], cij[3], P[3], Pi[3], Pj[3], PI[3], PI2;
//  int lmij[3];
//  int t1, u1, v1, t2, u2, v2;
//  double p, p2, mu, x, Ni, Nj, Nij;
//  double Hx1, Hy1, Hz1, HC_cef1;
//  double Hx2, Hy2, Hz2, HC_cef2; 
//  double element_ijkl = 0;
//
////  for(i=0;i<aoi;i++) i0 += ng[i];
////  for(j=0;j<aoj;j++) j0 += ng[j];
////  for(k=0;k<3;k++){
////    lmi[k] = lm_xyz[k + 3*aoi];
////    lmj[k] = lm_xyz[k + 3*aoj];
////    lmij[k] = lmi[k] + lmj[k];
////    ci[k] = center[k + 3*aoi];
////    cj[k] = center[k + 3*aoj];
////    cij[k] = ci[k] - cj[k];
////  }
////  for(i=0;i<ngi;i++){
////    expi[i] = exp[i0+i];
////    cefi[i] = cef[i0+i];
////    for(j=0;j<ngj;j++){
////      expj[j] = exp[j0+j];
////      cefj[j] = cef[j0+j];
////      p = expi[i] + expj[j];
////      p2 = 2*p;
////      mu = expi[i] * expj[j] / p;
////      Ni = cefi[i] * Norm(expi[i], lmi);
////      Nj = cefj[j] * Norm(expj[j], lmj);
////      Nij = Ni * Nj;
////      for(k=0;k<3;k++){
////        P[k] = (expi[i]*ci[k] + expj[j]*cj[k])/p;
////        Pi[k] = P[k] - ci[k];
////        Pj[k] = P[k] - cj[k];
////      }
////      for(I=0;I<Nao;I++){
////        PI2 = 0;
////        for(k=0;k<3;k++){
////          cI[k] = 0;
////          PI[k] = P[k] - cI[k];
////          PI2 += PI[k] * PI[k];
////        }
////        x = p*PI2;
////
////        for(t=0;t<lmij[0]+1;t++){
////          Hx = Hermite(lmi[0], lmj[0], t, p2, mu, 
////                       cij[0], Pi[0], Pj[0]);
////          for(u=0;u<lmij[1]+1;u++){
////            Hy = Hermite(lmi[1], lmj[1], u, p2, mu, 
////                         cij[1], Pi[1], Pj[1]);
////            for(v=0;v<lmij[2]+1;v++){
////              Hz = Hermite(lmi[2], lmj[2], v, p2, mu, 
////                           cij[2], Pi[2], Pj[2]);
////              HC_cef = HCcef(t, u, v, 0, PI, p2, x);
////              num = (Nij*Hx*Hy*Hz*HC_cef);
////              element_ij -= num*(2*M_PI/p);
////            }
////          }
////        }
////      }
////    }
////  }
//
//  free(expi);
//  free(expj);
//  free(expk);
//  free(expl);
//  free(cefi);
//  free(cefj);
//  free(cefk);
//  free(cefl);
//
//  return element_ijkl;
//}

/*********************
*  python interface  *
*********************/
static PyObject* gc2int(PyObject* self, PyObject* args){

  /*** input variables ***/
  /* dictionary */
  PyObject *in_dict;
  PyObject *dictList;
  PyObject *dictList_fast;
  PyObject *pyStr;
  PyObject *item;
  /* numpy array */
  PyArrayObject *in_array1, *in_array2;
  NpyIter *in_iter;
  NpyIter_IterNextFunc *in_iternext;
  /* C data type for input */
  int Ngo, Nao;
  int N, itr;
  /* basis function variables */
  double **in_ptr; // address to numpy pointer
  double *center;  // gaussian center
  double *cef;     // contraction coefficients
  int *ng;         // number of gaussians per AO
  double *exp;     // gaussian exponents
  int *lm_xyz;     // angular momentum

  /* python output variables */
  PyObject *py_out;
  PyObject *py_out2;
  double *data, *overlap;
  int i, j, k, l;
  int mat_dim[4];

  /*  parse numpy array and two integers as argument */
  if (!PyArg_ParseTuple(args, "OO!O!",
                        &in_dict, 
                        &PyArray_Type, &in_array1,
                        &PyArray_Type, &in_array2
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

  /*** end of numpy input data ***/

  /***** end of input data construction *****/


  /*******************
  * construct output *
  *******************/
  for(i=0;i<4;i++) mat_dim[i] = Nao;
  py_out = (PyArrayObject*) PyArray_FromDims(4, mat_dim, 
                                             NPY_DOUBLE);
  data = pyvector_to_Carrayptrs(py_out);

  /* renormalization */
  renormalize(center, exp, cef, ng, lm_xyz, Nao);

  /* orthogonalize */
  // no effect at the moment, for debug purpose
  //py_out2 = (PyArrayObject*) PyArray_FromDims(2, mat_dim, NPY_DOUBLE);
  //overlap = pyvector_to_Carrayptrs(py_out2);
  //orthogonalize(overlap, center, exp, cef, ng, lm_xyz, Nao);

  for(i=0;i<Nao;i++){
    for(j=i;j<Nao;j++){
      for(k=0;k<Nao;k++){
        for(l=k;l<Nao;l++){
          data[l + k*Nao + j*Nao*Nao + i*Nao*Nao*Nao] = 
            gc2Matrix(center, exp, cef, ng, lm_xyz, 
                      Nao, i, j, k, l);
          if(k!=l){
            data[k + l*Nao + j*Nao*Nao + i*Nao*Nao*Nao] =
            data[l + k*Nao + j*Nao*Nao + i*Nao*Nao*Nao];
          }
        }
      }
      if(j!=i){
        data[Nao*Nao + i*Nao*Nao + j*Nao*Nao*Nao] = 
        data[Nao*Nao + j*Nao*Nao + i*Nao*Nao*Nao];
      }
    }
  }

  //py_out = PyArray_SimpleNewFromData(2, mat_dim, NPY_DOUBLE, data);

  /*********************************
  * clean up and return the result *
  *********************************/
  free(exp);
  free(cef);
  free(ng);
  free(center);
  free(overlap);

  Py_INCREF(py_out);
  return Py_BuildValue("O", py_out);

  /*  in case bad things happen */
  fail:
    Py_XDECREF(py_out);
    return NULL;
}

/*  define functions in module */
static PyMethodDef GC2Int[] ={
  {"gc2int", gc2int, METH_VARARGS,
      "analytic two-center Gaussian-Coulomb integral"},
  {NULL, NULL, 0, NULL}
};

/* module initialization */
PyMODINIT_FUNC initgc2int(void){
  (void) Py_InitModule("gc2int", GC2Int);
  /* IMPORTANT: this must be called */
  import_array();
}
