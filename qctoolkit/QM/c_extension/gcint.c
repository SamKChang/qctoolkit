#include <Python.h>
#include <numpy/arrayobject.h>
#include <omp.h>
#include <gsl/gsl_sf_gamma.h>

/* Hermite-Gaussian expansion coefficients */
double Hermite(int mi, int mj, int t, double p2, double mu,
               double Xij, double Xpi, double Xpj){
  double Hij, H0, H1, H2;
  if((t>mi+mj)||(t<0)||(mi<0)||(mj<0)){
    Hij = 0;
  }else if((mi>=mj)&&(mi>0)){
    H0 = Hermite(mi-1,mj,t-1,p2,mu,Xij,Xpi,Xpj)/p2;
    H1 = Hermite(mi-1,mj,t,p2,mu,Xij,Xpi,Xpj)*Xpi;
    H2 = Hermite(mi-1,mj,t+1,p2,mu,Xij,Xpi,Xpj)*(t+1);
    Hij = H0 + H1 + H2;
  }else if((mj>mi)&&(mj>0)){
    H0 = Hermite(mi,mj-1,t-1,p2,mu,Xij,Xpi,Xpj)/p2;
    H1 = Hermite(mi,mj-1,t,p2,mu,Xij,Xpi,Xpj)*Xpj;
    H2 = Hermite(mi,mj-1,t+1,p2,mu,Xij,Xpi,Xpj)*(t+1);
    Hij = H0 + H1 + H2;
  }else{
    Hij = exp(-mu*pow(Xij,2));
  }
  return Hij;
}

/* factorial function */
int fac(n){
  if (n<-1) return 0;
  else if (n<=0) return 1;
  else{
    int val = 1, k = n;
    while(k>0){
      val *= k;
      k = k-1;
    }
    return val;
  }
}
/* 2-factorial function */
int fac2(n){
  if (n<-1) return 0;
  else if (n<=0) return 1;
  else{
    int val = 1, k = n;
    while(k>0){
      val *= k;
      k = k-2;
    }
    return val;
  }
}

/* Gaussian normalization constant */
double Norm(float exp, int* lm){
  int nx = fac2(2*lm[0]-1);
  int ny = fac2(2*lm[1]-1);
  int nz = fac2(2*lm[2]-1);
  int n_xyz = lm[0] + lm[1] + lm[2];
  double N_xyz = nx * ny * nz;
  double num1 = pow(2*exp/M_PI, 1.5);
  double num2 = pow(4*exp, n_xyz);
  double out = pow(num1 * num2 / N_xyz, 0.5);
  return out;
}

/* Boys function */
double F(int n, double x){
  if(x>0){
    double g1, g2, factor;
    factor = 0.5 * pow(x, -n-0.5);
    g1 = gsl_sf_gamma(0.5+n);
    g2 = gsl_sf_gamma_inc_P(0.5+n, x);
    return factor * g1 * g2;
  }else{
    return 1.0/(1.0+2.0*n);
  }
}

/* Hermite-Coulomb integral coefficient */
double HCcef(int t, int u, int v, int n, 
             double* PI, double p2, double x){
  double cef = 0;
  if((t==0)&&(u==0)&&(v==0)){
    cef += F(n, x) * pow(-p2, n);
  }else if(((t>=u)||(t>=v))&&(t>0)){
    cef += HCcef(t-1,u,v,n+1,PI,p2,x) * PI[0];
    cef += HCcef(t-2,u,v,n+1,PI,p2,x) * (t-1);
  }else if((u>=v)&&(u>0)){
    cef += HCcef(t,u-1,v,n+1,PI,p2,x) * PI[1];
    cef += HCcef(t,u-2,v,n+1,PI,p2,x) * (u-1);
  }else if(v>0){
    cef += HCcef(t,u,v-1,n+1,PI,p2,x) * PI[2];
    cef += HCcef(t,u,v-2,n+1,PI,p2,x) * (v-1);
  }
  return cef;
}

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
  double p, p2, mu, q, ZI, x, HC_cef, Ni, Nj, Nij;
  double Hx, Hy, Hz, num, element_ij = 0;

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

  for(i=0;i<ngi;i++){
    for(j=0;j<ngj;j++){
      p = expi[i] + expj[j];
      p2 = 2*p;
      mu = expi[i] * expj[j] / p;
      Ni = cefi[i] * Norm(expi[i], lmi);
      Nj = cefj[j] * Norm(expj[j], lmj);
      Nij = Ni * Nj;
      for(k=0;k<3;k++){
        P[k] = (expi[k]*ci[k] + expj[k]*cj[k])/p;
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

  /* check for data consistancy, for debug purpose */
//  for(i=0;i<Nao;i++){
//    printf("%d\n", ng[i]);
//  }


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
