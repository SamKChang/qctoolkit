// C routine to calculate ESP for a point

#include <Python.h>
#include <numpy/arrayobject.h>
#include <xc.h>

/////////////////////////////////////////
// main function for writing CUBE file //
/////////////////////////////////////////
// or what ever computation need to be done
// only minor changes should be necessary
// to interface other type of calculation
static void libxc_c(double *rho, double *sigma, 
                    double *vrho, double *vsigma,
                    int N, int xcID){
  xc_func_type func;
  int i;
  if(xc_func_init(&func, xcID, XC_UNPOLARIZED) != 0){
    fprintf(stderr, "Functional '%d' not found\n", xcID);
    return 1;
  }

  switch(func.info->family){
    case XC_FAMILY_LDA:
      xc_lda_vxc(&func, N, rho, vrho);
      for(i=0;i<N;i++) vsigma[i] = 0;
      break;
    case XC_FAMILY_GGA:
      xc_gga_vxc(&func, N, rho, sigma, vrho, vsigma);
      break;
    case XC_FAMILY_HYB_GGA:
      xc_gga_vxc(&func, N, rho, sigma, vrho, vsigma);
      break;
  }
  xc_func_end(&func);
}

static void nplist_to_C_double_array(PyArrayObject *in_array, 
                              double *out_list,
                              int N)
{
  NpyIter *in_iter;
  NpyIter_IterNextFunc *in_iternext;
  double **in_ptr;
  int itr;

  in_iter = NpyIter_New(in_array, NPY_ITER_READONLY,
                        NPY_KEEPORDER, NPY_NO_CASTING, NULL);
  if (in_iter == NULL) goto fail;
  in_iternext = NpyIter_GetIterNext(in_iter, NULL);
  if (in_iternext == NULL){
    NpyIter_Deallocate(in_iter);
    goto fail;
  }
  /* interator pointer to actual numpy data */
  in_ptr = (double **) NpyIter_GetDataPtrArray(in_iter);
  itr=0;
  do {
    out_list[itr++] = **in_ptr;
  } while(in_iternext(in_iter));
  NpyIter_Deallocate(in_iter);

  fail:
    return NULL;
}

//////////////////////////////////////////////
// redirect PyArray data contant to C-array //
//////////////////////////////////////////////
static double *pyvector_to_Carrayptrs(PyArrayObject *arrayin){
  int n=arrayin->dimensions[0];
  /* pointer to arrayin data as double */
  return (double *) arrayin->data;
}

// !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*! //
///////////////////////////////////////
// !!! PYTHON INTERFACE FROM HERE!!! //
///////////////////////////////////////
// !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*! //

//////////////////////////////
// python callable function //
//////////////////////////////
static PyObject * libxc_vxc(PyObject * self, PyObject * args){
  PyArrayObject *py_rho;
  PyArrayObject *py_sigma;
  PyArrayObject *py_vrho_out;
  PyArrayObject *py_vsigma_out;
  double *rho;
  double *sigma;
  double *vrho_out;
  double *vsigma_out;
  int N;
  int dims[1];
  int xcID;

  // parse arguments check and/or error handling
  if(!PyArg_ParseTuple(args, "O!O!ii", 
                       &PyArray_Type, &py_rho,
                       &PyArray_Type, &py_sigma,
                       &N,
                       &xcID
                      )) return NULL;
  if((py_rho == NULL)||(py_sigma == NULL))
    return NULL;

  /*** access numpy array data ***/
  rho = (double *) malloc(N * sizeof(double));
  sigma = (double *) malloc(N * sizeof(double));
  nplist_to_C_double_array(py_rho, rho, N);
  nplist_to_C_double_array(py_sigma, sigma, N);
  dims[0] = N;

  py_vrho_out = (PyArrayObject*) PyArray_FromDims(1, dims, NPY_DOUBLE);
  py_vsigma_out=(PyArrayObject*) PyArray_FromDims(1, dims, NPY_DOUBLE);
  vrho_out = pyvector_to_Carrayptrs(py_vrho_out);
  vsigma_out = pyvector_to_Carrayptrs(py_vsigma_out);
  /*** end of list input data ***/

  libxc_c(rho, sigma, vrho_out, vsigma_out, N, xcID);

  free(rho);
  free(sigma);
  // freeing xcString gives error...

  return Py_BuildValue("OO", py_vrho_out, py_vsigma_out);
}

/////////////////////////////////////////
// register python module symbol table //
/////////////////////////////////////////
// PyMethodDef: struct of four field
//   string: method_name
//   PyCFunction: method_function
//   int: flag
//   string: documentation
static PyMethodDef libxc_vxcMethods[] = {
  {"libxc_vxc", libxc_vxc, METH_VARARGS, "libxc_vxc interface"},
  {NULL, NULL, 0, NULL} // sentinel?
};

//////////////////////////
//  method constructor  //
//////////////////////////
// the name MUST be init{name} otherwise python cannot find it
// depends on numpy C-API, import_array()
// and/or import_ufunc() are necessary
// otherwise the code return segfault
PyMODINIT_FUNC initlibxc_vxc(void){
  Py_InitModule("libxc_vxc", libxc_vxcMethods);
  import_array(); // necessary!
}
