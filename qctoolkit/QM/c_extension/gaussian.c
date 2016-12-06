#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_erf.h>
#include <math.h>
#include "gaussian.h"

void setZero(double* smallNumber){
  if(*smallNumber * *smallNumber < 1E-20)
    *smallNumber = 0.0;
}

void setD2Cef(int s, int lm, double exp, double* cef){
  if (s == 2) *cef = 4 * exp * exp;
  else if (s == 0) *cef = -2 * exp * (2 * lm + 1);
  else if (s == -2) *cef = lm * (lm - 1);
  else if (s == 1) *cef = -2 * exp;
  else if (s == -1) *cef = lm;
}

/* 2-factorial function */
// approximation necessary for large n
int fac2(int n){
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

/*******************************************
*  Hermite-Gaussian expansion coefficients *
*******************************************/
// used to evaluate Gaussian-Gaussian overlap integral
double Hermite(int mi, int mj, int t, double p2, double mu,
               double Xij, double Xpi, double Xpj){
  double Hij, H0, H1, H2;
  if((t>mi+mj)||(t<0)||(mi<0)||(mj<0)){
    Hij = 0;
  }else if((mi>=mj)&&(mi>0)){
    // annihilate higher order of first Gaussian
    H0 = Hermite(mi-1,mj,t-1,p2,mu,Xij,Xpi,Xpj)/p2;
    H1 = Hermite(mi-1,mj,t,p2,mu,Xij,Xpi,Xpj)*Xpi;
    H2 = Hermite(mi-1,mj,t+1,p2,mu,Xij,Xpi,Xpj)*(t+1);
    Hij = H0 + H1 + H2;
  }else if((mj>mi)&&(mj>0)){
    // annihilate higher order of second Gaussian
    H0 = Hermite(mi,mj-1,t-1,p2,mu,Xij,Xpi,Xpj)/p2;
    H1 = Hermite(mi,mj-1,t,p2,mu,Xij,Xpi,Xpj)*Xpj;
    H2 = Hermite(mi,mj-1,t+1,p2,mu,Xij,Xpi,Xpj)*(t+1);
    Hij = H0 + H1 + H2;
  }else{
    Hij = exp(-mu*pow(Xij,2));
  }
  return Hij;
}

/************************************
*  Gaussian normalization constant  *
************************************/
// each Gaussian orbital is assumed to be normalized
// the normalization constant is determined by
// its order and exponents
// therefore any single Gaussian AO has coefficient=1.000
// the normalization constant is assumed as common knowledge
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
// calling gsl gamma/error function library
double F(int n, double x){
// gamma function implementation
// directly calculate top n value
  if(x>0){
    // analytic expression of boys function from gamma functions
    // result from Mathematica
    double g1, g2, factor;
    g1 = gsl_sf_gamma(0.5+n);
    g2 = gsl_sf_gamma_inc_P(0.5+n, x);
    factor = 0.5 * pow(x, -n-0.5);
    return factor * g1 * g2;
  }else{
    return 1.0/(1.0+2.0*n);
  }

//// approximated Boy's funciton, 
//// very similar results as gamma function implimentation
//// test on e-e energy of H2O with aug-cc-pvdz basis 
//// the difference is less than 1E-6 Ha
//// but the compute time is considerably longer
//  int i;
//  double expx, gamma_n, term, crit, Fn = 0;
//
//  expx = exp(-x);
//  gamma_n = gsl_sf_gamma(n+0.5);
//  crit = 1.0E-12;
//
//  if(x<=0){
//    Fn = 1.0/(1.0+2.0*n);
//  }else if(x<=10){
//    i = 0;
//    while(i<100){
//      term = pow(x, i) / gsl_sf_gamma(n+i+1.5);
//      term *= gamma_n;
//      if(term/Fn < crit){
//        break;
//      }else if(i > 100){
//        printf("Boy's function not converged for x<10\n");
//      }
//      Fn += term;
//      i++;
//    }
//    Fn *= 0.5 * expx;
//  }else{
//    i = 0;
//    while(i<100){
//      term = pow(x, -i) /  gsl_sf_gamma(n-i+1.5);
//      term *= gamma_n;
//      if(term/Fn < crit){
//        break;
//      }else if(i > 100){
//        printf("Boy's function not converged for x>10\n");
//      }
//      Fn += term;
//      i++;
//    }
//    Fn *= -0.5 * expx;
//    Fn += gamma_n / (2*pow(x, n+0.5));
//  }
//  return Fn;
}

/* Hermite-Coulomb integral coefficient */
double HCcef(int t, int u, int v, int n, 
             double* PI, double p2, double x){
  double cef = 0;
  if((t==0)&&(u==0)&&(v==0)){
    cef += F(n, x) * pow(-p2, n);
  }else if(((t>=u)||(t>=v))&&(t>0)){
    // annihilate x-direction (t)
    cef += HCcef(t-1,u,v,n+1,PI,p2,x) * PI[0];
    cef += HCcef(t-2,u,v,n+1,PI,p2,x) * (t-1);
  }else if((u>=v)&&(u>0)){
    // annihilate y-direction (u)
    cef += HCcef(t,u-1,v,n+1,PI,p2,x) * PI[1];
    cef += HCcef(t,u-2,v,n+1,PI,p2,x) * (u-1);
  }else if(v>0){
    // annihilate z-direction (v)
    cef += HCcef(t,u,v-1,n+1,PI,p2,x) * PI[2];
    cef += HCcef(t,u,v-2,n+1,PI,p2,x) * (v-1);
  }
  return cef;
}

/********************************************
*  Overlap integral of two atomic orbitals  *
********************************************/
// aoverlap = aoOverlap(center, exp, cef, ng, lm_xyz, aoi, aoj
//                      flag, t1, u1, v1, t2, v2, u2);
// t1, u1, v1, t2, v2, u2 are used for kinetic energy integrals
// it must be set to zero for normalization overlap integral
double aoOverlap(double *center, double *exp, double *cef, 
                 int *ng,int *lm_xyz, int aoi, int aoj, int flag,
                 int t1, int u1, int v1, int t2, int u2, int v2){
  int ngi = ng[aoi], ngj = ng[aoj];
  int i, j, k, lmi[3], lmj[3], i0 = 0, j0 = 0;
  double p, p2, mu, cef_out;
  double ci[3], cj[3], cij[3];
  double P[3], Pi[3], Pj[3];
  double norm;
  double Hx, Hy, Hz, overlap = 0.0;
  double *expi, *expj, *cefi, *cefj;
  double factor;

  expi = (double*) malloc(ngi * sizeof(double));
  expj = (double*) malloc(ngj * sizeof(double));
  cefi = (double*) malloc(ngi * sizeof(double));
  cefj = (double*) malloc(ngj * sizeof(double));

	for(i=0;i<aoi;i++) i0 += ng[i];
  for(j=0;j<aoj;j++) j0 += ng[j];
  for(k=0;k<3;k++){
    lmi[k] = lm_xyz[k + aoi * 3];
    lmj[k] = lm_xyz[k + aoj * 3];
    ci[k] = center[k + aoi * 3];
    cj[k] = center[k + aoj * 3];
    cij[k] = ci[k] - cj[k];
  }

  for(i = 0; i < ngi; i++){
    expi[i] = exp[i0 + i];
    cefi[i] = cef[i0 + i];
    for(j = 0; j < ngj; j++){
      expj[j] = exp[j0 + j];
      cefj[j] = cef[j0 + j];
      cef_out = cefi[i] * cefj[j];
      norm = Norm(expi[i], lmi) * Norm(expj[j], lmj);
      p = expi[i] + expj[j];
      p2 = 2*p;
      mu = expi[i] * expj[j] / p;
      for(k = 0; k < 3; k++){
        P[k] = (expi[i] * ci[k] + expj[j] * cj[k]) / p;
        Pi[k] = P[k] - ci[k];
        Pj[k] = P[k] - cj[k];
      }
      Hx = Hermite(lmi[0] + t1, lmj[0] + t2, 0, p2, mu, 
                   cij[0], Pi[0], Pj[0]);
      Hy = Hermite(lmi[1] + u1, lmj[1] + u2, 0, p2, mu, 
                   cij[1], Pi[1], Pj[1]);
      Hz = Hermite(lmi[2] + v1, lmj[2] + v2, 0, p2, mu, 
                   cij[2], Pi[2], Pj[2]);
      factor = cef_out * norm * Hx*Hy*Hz * pow(M_PI/p, 1.5);

      if (flag != 0){
        double cx = 1.0, cy = 1.0, cz = 1.0;       
        if(flag == 1){
          setD2Cef(t2, lmj[0], expj[j], &cx);
        } else if(flag == 2) {
          setD2Cef(u2, lmj[1], expj[j], &cy);
        } else if(flag == 3) {
          setD2Cef(v2, lmj[2], expj[j], &cz);
        }
        factor *= cx * cy * cz;
      }
      overlap += factor;
    }
  }

  free(expi);
  free(expj);
  free(cefi);
  free(cefj);

  return overlap;
}

/************************************
*  Renormalization of AO expansion  *
************************************/
// renormalize(center, exp, cef, ng, lm_xyz, Nao);
void renormalize(double *center, double *exp, double *cef, 
                 int *ng, int *lm_xyz, int Nao){
  int ao, i, ngo=0;
  double overlap;

  for(ao=0;ao<Nao;ao++){
    overlap = aoOverlap(center, exp, cef, ng, lm_xyz, 
                        ao, ao, 0, 0, 0, 0, 0, 0, 0);
    for(i=ngo;i<ngo+ng[ao];i++){
      cef[i] /= sqrt(overlap);
    }
    ngo += ng[ao];
  }
}

/***********************************
*  Density normalization integral  *
***********************************/
// nint = densityIntegral(center, exp, cef, ng, lm_xyz, ao);
double densityIntegral(double *center, double *exp, double *cef, 
                       int *ng, int *lm_xyz, int ao){
  int i, s, lm[3], i0 = 0;
  double p, p2, mu, cef_out;
  double norm;
  double Hx, Hy, Hz, nint = 0;
  double *expi, *cefi;

  expi = (double*) malloc(ng[ao] * sizeof(double));
  cefi = (double*) malloc(ng[ao] * sizeof(double));

	for(i=0;i<ao;i++) i0 += ng[i];
  for(s=0;s<3;s++) lm[s] = lm_xyz[s+ao*3];

  for(i=0;i<ng[ao];i++){
    expi[i] = exp[i0+i];
    cefi[i] = cef[i0+i];
    cef_out = cefi[i];
    norm = Norm(expi[i], lm);
    p = expi[i];
    p2 = 2*p;
    mu = p;
    Hx = Hermite(lm[0], 0, 0, p2, mu, 0, 0, 0);
    Hy = Hermite(lm[1], 0, 0, p2, mu, 0, 0, 0);
    Hz = Hermite(lm[2], 0, 0, p2, mu, 0, 0, 0);
    nint += cef_out * norm * Hx*Hy*Hz * pow(M_PI/p, 1.5);
  }

  free(expi);
  free(cefi);

  return nint;
}

/*****************************************
*  Renormalization of density expansion  *
*****************************************/
// densityRenormalize(center, exp, cef, ng, lm_xyz, Nao);
void densityRenormalize(double *center, double *exp, double *cef, 
                        int *ng, int *lm_xyz, int Nao){
  int ao, i, ngo=0;
  double nint;

  for(ao=0;ao<Nao;ao++){
    nint = densityIntegral(center, exp, cef, ng, lm_xyz, ao);
    for(i=ngo;i<ngo+ng[ao];i++){
      cef[i] /= nint;
    }
    ngo += ng[ao];
  }
}

///****************************************
//*  Orthogonalization of overlap matrix  *
//****************************************/ 
//// calculate and diagonalize overlap matrix using LAPACK functions
//// at the moment not used
//// orthogonalize(overlap, center, exp, cef, ng, lm_xyz, Nao);
//void orthogonalize(double *overlap, double *center, 
//                   double *exp, double *cef, int *ng, 
//                   int *lm_xyz, int Nao){
//  int i, j;
//
//  printf("before diagonalization\n");
//  for(i=0;i<Nao;i++){
//    for(j=0;j<Nao;j++){
//      overlap[j+i*Nao] = aoOverlap(center, exp, cef, 
//                                   ng, lm_xyz, i, j,
//                                   0, 0, 0, 0, 0, 0, 0);
//      if(pow(overlap[j+i*Nao],2) > 0){
//        printf("S(%d, %d) = %10.8f\n", i+1, j+1, overlap[j+i*Nao]);
//      }
//    }
//  }
//
//  /* LAPACK variables */
//  int n = Nao, lda = Nao, info, lwork = -1;
//  double wkopt;
//  double *work;
//  double *w = (double *) malloc(Nao * Nao * sizeof(double));
//  dsyev_("Vectors", "Upper", &n, overlap, 
//         &lda, w, &wkopt, &lwork, &info);
//  lwork = (int)wkopt;
//  work = (double*)malloc(lwork*sizeof(double));
//  dsyev_("Vectors", "Upper", &n, overlap, 
//         &lda, w, work, &lwork, &info);
//  /* end of LAPACK routine */
//
////  printf("after diagonalization\n");
////  for(i=0;i<Nao;i++){
////    for(j=0;j<Nao;j++){
////      if(overlap[j+i*Nao] > 0){
////        printf("S(%d, %d) = %f\n", i, j, overlap[j+i*Nao]);
////      }
////    }
////  }
////  for(i=0;i<Nao;i++) printf(" %f", w[i]);
////  printf("\n");
//  free(w);
//  free(work);
//}

/************************************************
*  main function for Gaussian-Coulomb integral  *
************************************************/
//veMatrix(R, Z, center, exp, cef, ng, 
//         lm_xyz, Nao, N, i, j);
double veMatrix(double *R,      //all the rest are input
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
  double *expi, *expj, *cefi, *cefj;
  expi = (double*) malloc(ngi * sizeof(double));
  expj = (double*) malloc(ngj * sizeof(double));
  cefi = (double*) malloc(ngi * sizeof(double));
  cefj = (double*) malloc(ngj * sizeof(double));

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
            }
          }
        }
      }
    }
  }

  free(expi);
  free(expj);
  free(cefi);
  free(cefj);
  return element_ij;
}
/************************************************
*  main function for Gaussian-Coulomb integral  *
************************************************/
//eeKernel(R, Z, center, exp, cef, ng, 
//         lm_xyz, Nao, N, i, j);
void eeKernel(double *R,      //all the rest are input
              double *Z, 
              double *center,
              double *exp,
              double *cef,
              double *element_ij,
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
  double *expi, *expj, *cefi, *cefj;
  expi = (double*) malloc(ngi * sizeof(double));
  expj = (double*) malloc(ngj * sizeof(double));
  cefi = (double*) malloc(ngi * sizeof(double));
  cefj = (double*) malloc(ngj * sizeof(double));

  /* gaussian integral variables */
  double cI[3], cij[3], P[3], Pi[3], Pj[3], PI[3], PI2;
  int lmij[3];
  int t, u, v;
  double p, p2, mu, ZI, x, HC_cef, Ni, Nj, Nij;
  double Hx, Hy, Hz, num, element_ijI;

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
  for(I=0;I<N;I++){
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
              element_ij[I] += num*(2*M_PI/p);
            }
          }
        }
      }
    }
  }

  free(expi);
  free(expj);
  free(cefi);
  free(cefj);
}

/************************************************
*  main function for Gaussian-Coulomb integral  *
************************************************/
//vnMatrix(R, Z, center, exp, cef, ng, 
//         lm_xyz, Nao, N, i, j);
double vnMatrix(double *R,      //all the rest are input
                double *Z, 
                double *center,
                double *exp,
                double *cef,
                int *ng,
                int *lm_xyz,
                int Nao,
                int N,
                int ana
             ){  

  /* input related variables */
  int s, a, a0 = 0, I;
  double ca[3];
  int lma[3];
  int nga = ng[ana];
  double *expa, *cefa;
  expa = (double*) malloc(nga * sizeof(double));
  cefa = (double*) malloc(nga * sizeof(double));

  /* gaussian integral variables */
  double cI[3], P[3], PI[3], PI2;
  int t, u, v;
  double p, p2, mu, ZI, x, HC_cef, Na;
  double Hx, Hy, Hz, num, element_a = 0;

  for(s=0;s<ana;s++) a0 += ng[s];
  for(s=0;s<3;s++){
    lma[s] = lm_xyz[s + 3*ana];
    ca[s] = center[s + 3*ana];
  }
  for(a=0;a<nga;a++){
    expa[a] = exp[a0+a];
    cefa[a] = cef[a0+a];
    p = expa[a];
    p2 = 2*p;
    mu = p;
    Na = cefa[a] * Norm(expa[a], lma);
    for(s=0;s<3;s++){
      P[s] = ca[s];
    }
    for(I=0;I<N;I++){
      PI2 = 0;
      for(s=0;s<3;s++){
        cI[s] = R[s + 3*I];
        PI[s] = P[s] - cI[s];
        PI2 += PI[s] * PI[s];
      }
      ZI = Z[I];
      x = p*PI2;

      for(t=0;t<lma[0]+1;t++){
        Hx = Hermite(lma[0], 0, t, p2, mu,
                     0, 0, 0);
        for(u=0;u<lma[1]+1;u++){
          Hy = Hermite(lma[1], 0, u, p2, mu,
                       0, 0, 0);
          for(v=0;v<lma[2]+1;v++){
            Hz = Hermite(lma[2], 0, v, p2, mu,
                         0, 0, 0);
            HC_cef = HCcef(t, u, v, 0, PI, p2, x);
            num = ZI*(Na*Hx*Hy*Hz*HC_cef);
            element_a -= num*(2*M_PI/p);
          }
        }
      }
    }
  }

  free(expa);
  free(cefa);
  return element_a;
}

/************************************************
*  main function for Gaussian-Coulomb integral  *
************************************************/
// calculate 4D array for 2-center Gaussian-Coulomb
// electron-electron repulsion matrix from atomic Gaussian orbitals
// NOTE: Exchange interaction need to be considered 
//       to compute two-electron energy
// eeMatrix(center, exp, cef, ng, 
//          lm_xyz, Nao, N, i, j, k, l);
double eeMatrix(double *center,
                double *exp,
                double *cef,
                int *ng,
                int *lm_xyz,
                int Nao,
                int aoi,
                int aoj,
                int aok,
                int aol
              ){

  /* input related variables */
  int s, i, j, k, l, i0 = 0, j0 = 0, k0 = 0, l0 = 0;
  double ci[3], cj[3], ck[3], cl[3];
  int lmi[3], lmj[3], lmk[3], lml[3];
  int ngi = ng[aoi], ngj = ng[aoj], ngk = ng[aok], ngl = ng[aol];
  double *expi, *expj, *expk, *expl;
  double *cefi, *cefj, *cefk, *cefl;
  expi = (double*) malloc(ngi * sizeof(double));
  expj = (double*) malloc(ngj * sizeof(double));
  expk = (double*) malloc(ngk * sizeof(double));
  expl = (double*) malloc(ngl * sizeof(double));
  cefi = (double*) malloc(ngi * sizeof(double));
  cefj = (double*) malloc(ngj * sizeof(double));
  cefk = (double*) malloc(ngk * sizeof(double));
  cefl = (double*) malloc(ngl * sizeof(double));

  /* gaussian integral variables */
  double P[3], Pi[3], Pj[3];
  double Q[3], Qk[3], Ql[3];
  double PQ[3], PQ2, alpha;
  double cij[3], ckl[3];
  int lmij[3], lmkl[3];
  int t1, u1, v1, t2, u2, v2;
  double p, p2, mu, Ni, Nj, Nij;
  double q, q2, nu, Nk, Nl, Nkl;
  double x;
  double Hx1, Hy1, Hz1;
  double Hx2, Hy2, Hz2;
  double factor1, factor2, HC_cef; 
  double element_ijkl = 0;

  /* setup variables and save to allocated memory */
  for(s=0;s<aoi;s++) i0 += ng[s];
  for(s=0;s<aoj;s++) j0 += ng[s];
  for(s=0;s<aok;s++) k0 += ng[s];
  for(s=0;s<aol;s++) l0 += ng[s];
  for(s=0;s<3;s++){
    lmi[s] = lm_xyz[s + 3*aoi];
    lmj[s] = lm_xyz[s + 3*aoj];
    lmk[s] = lm_xyz[s + 3*aok];
    lml[s] = lm_xyz[s + 3*aol];
    lmij[s] = lmi[s] + lmj[s];
    lmkl[s] = lmk[s] + lml[s];
    ci[s] = center[s + 3*aoi];
    cj[s] = center[s + 3*aoj];
    ck[s] = center[s + 3*aok];
    cl[s] = center[s + 3*aol];
    cij[s] = ci[s] - cj[s];
    ckl[s] = ck[s] - cl[s];
  }
  /* loop for i_th AO */
  for(i=0;i<ngi;i++){
    expi[i] = exp[i0+i];
    cefi[i] = cef[i0+i];
    Ni = cefi[i] * Norm(expi[i], lmi);
    /* loop for j_th AO */
    for(j=0;j<ngj;j++){
      expj[j] = exp[j0+j];
      cefj[j] = cef[j0+j];
      Nj = cefj[j] * Norm(expj[j], lmj);
      /* variables for first spatial variale, r1 */
      Nij = Ni * Nj;
      p = expi[i] + expj[j];
      p2 = 2*p;
      mu = expi[i] * expj[j] / p;
      for(s=0;s<3;s++){
        P[s] = (expi[i]*ci[s] + expj[j]*cj[s]) / p;
        Pi[s] = P[s] - ci[s];
        Pj[s] = P[s] - cj[s];
      }
      /* loop for k_th AO */
      for(k=0;k<ngk;k++){
        expk[k] = exp[k0+k];
        cefk[k] = cef[k0+k];
        Nk = cefk[k] * Norm(expk[k], lmk);
        /* loop for l_th AO */
        for(l=0;l<ngl;l++){
          expl[l] = exp[l0+l];
          cefl[l] = cef[l0+l];
          Nl = cefl[l] * Norm(expl[l], lml);
          /* variables for second spatial variale, r2 */
          Nkl = Nk * Nl;
          q = expk[k] + expl[l];
          q2 = 2*q;
          nu = expk[k] * expl[l] / q;
          factor1 = 2*pow(M_PI, 5.0/2.0);
          factor1 /= p*q * sqrt(p+q);
          PQ2 = 0;
          for(s=0;s<3;s++){
            Q[s] = (expk[k]*ck[s] + expl[l]*cl[s]) / q;
            Qk[s] = Q[s] - ck[s];
            Ql[s] = Q[s] - cl[s];
            PQ[s] = P[s] - Q[s];
            setZero(&PQ[s]);
            PQ2 += PQ[s] * PQ[s];
          }
          alpha = p*q / (p+q);
          x = alpha * PQ2;
          /* Hermite coefficient for x1 */
          for(t1=0;t1<lmij[0]+1;t1++){
            Hx1 = Hermite(lmi[0], lmj[0], t1, p2, mu, 
                          cij[0], Pi[0], Pj[0]);
            /* Hermite coefficient for y1 */
            for(u1=0;u1<lmij[1]+1;u1++){
              Hy1 = Hermite(lmi[1], lmj[1], u1, p2, mu, 
                            cij[1], Pi[1], Pj[1]);
              /* Hermite coefficient for z1 */
              for(v1=0;v1<lmij[2]+1;v1++){
                Hz1 = Hermite(lmi[2], lmj[2], v1, p2, mu, 
                              cij[2], Pi[2], Pj[2]);
                /* Hermite coefficient for x2 */
                for(t2=0;t2<lmkl[0]+1;t2++){
                  Hx2 = Hermite(lmk[0], lml[0], t2, q2, nu, 
                                ckl[0], Qk[0], Ql[0]);
                  /* Hermite coefficient for y2 */
                  for(u2=0;u2<lmkl[1]+1;u2++){
                    Hy2 = Hermite(lmk[1], lml[1], u2, q2, nu, 
                                  ckl[1], Qk[1], Ql[1]);
                    /* Hermite coefficient for z2 */
                    for(v2=0;v2<lmkl[2]+1;v2++){
                      Hz2 = Hermite(lmk[2], lml[2], v2, q2, nu, 
                                    ckl[2], Qk[2], Ql[2]);
                      factor2 = Nij*Nkl*Hx1*Hy1*Hz1*Hx2*Hy2*Hz2;
                      factor2 *= pow(-1, t2+u2+v2);
                      /* 2-center integral */
                      HC_cef = HCcef(t1+t2, u1+u2, v1+v2, 
                                     0, PQ, 2*alpha, x);
                      element_ijkl += factor1*factor2*HC_cef;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  free(expi);
  free(expj);
  free(expk);
  free(expl);
  free(cefi);
  free(cefj);
  free(cefk);
  free(cefl);

  return element_ijkl;
}

/****************************************************************
*  main function for Gaussian-Coulomb density-fitting integral  *
****************************************************************/
// calculate 3D array for 2-center Gaussian-Coulomb
// electron-electron repulsion matrix from atomic Gaussian orbitals
// expanded by density basis and orbital basis
// It return 3D array, OUT[a, i, j] where a is the basis index
// for density gaussian expansion while i,j denote basis index
// for orbital gaussian expansion
// NOTE: Exchange interaction need to be considered 
//       to compute two-electron energy
// neMatrix(data, center, exp, cef, ng, 
//          lm_xyz, Nao, N, i, j, k, l);
double neMatrix(double *center,
                double *exp,
                double *cef,
                int *ng,
                int *lm_xyz,
                int Nao,
                double *fcenter,
                double *fexp,
                double *fcef,
                int *fng,
                int *flm_xyz,
                int fNao,
                int ana,
                int aoi,
                int aoj
              ){

  /* input related variables */
  int s, a, i, j, a0 = 0, i0 = 0, j0 = 0;
  double ca[3], ci[3], cj[3];
  int lma[3], lmi[3], lmj[3];
  int nga = fng[ana], ngi = ng[aoi], ngj = ng[aoj];
  double *expa, *expi, *expj, *cefa, *cefi, *cefj;
  expa = (double*) malloc(nga * sizeof(double));
  expi = (double*) malloc(ngi * sizeof(double));
  expj = (double*) malloc(ngj * sizeof(double));
  cefa = (double*) malloc(nga * sizeof(double));
  cefi = (double*) malloc(ngi * sizeof(double));
  cefj = (double*) malloc(ngj * sizeof(double));

  /* gaussian integral variables */
  double Q[3];
  double P[3], Pi[3], Pj[3];
  double PQ[3], PQ2, alpha;
  double cij[3];
  int lmij[3];
  int t1, u1, v1, t2, u2, v2;
  double q, q2, nu, Na;
  double p, p2, mu, Ni, Nj, Nij;
  double x;
  double Hx1, Hy1, Hz1;
  double Hx2, Hy2, Hz2;
  double factor1, factor2, HC_cef; 
  double element_aij = 0;

  /* setup variables and save to allocated memory */
  for(s=0;s<ana;s++) a0 += fng[s];
  for(s=0;s<aoi;s++) i0 += ng[s];
  for(s=0;s<aoj;s++) j0 += ng[s];
  for(s=0;s<3;s++){
    lma[s] = flm_xyz[s + 3*ana];
    lmi[s] = lm_xyz[s + 3*aoi];
    lmj[s] = lm_xyz[s + 3*aoj];
    lmij[s] = lmi[s] + lmj[s];
    ca[s] = fcenter[s + 3*ana];
    ci[s] = center[s + 3*aoi];
    cj[s] = center[s + 3*aoj];
    cij[s] = ci[s] - cj[s];
  }
  /* loop for a_th gaussian density */
  for(a=0;a<nga;a++){
    expa[a] = fexp[a0+a];
    cefa[a] = fcef[a0+a];
    for(s=0;s<3;s++) Q[s] = ca[s];
    Na = cefa[a] * Norm(expa[a], lma);
    q = expa[a];
    q2 = 2*q;
    nu = q;
    /* loop for i_th AO */
    for(i=0;i<ngi;i++){
      expi[i] = exp[i0+i];
      cefi[i] = cef[i0+i];
      Ni = cefi[i] * Norm(expi[i], lmi);
      /* loop for j_th AO */
      for(j=0;j<ngj;j++){
        expj[j] = exp[j0+j];
        cefj[j] = cef[j0+j];
        Nj = cefj[j] * Norm(expj[j], lmj);
        /* variables for first spatial variale, r1 */
        Nij = Ni * Nj;
        p = expi[i] + expj[j];
        p2 = 2*p;
        mu = expi[i] * expj[j] / p;
        PQ2 = 0;
        for(s=0;s<3;s++){
          P[s] = (expi[i]*ci[s] + expj[j]*cj[s]) / p;
          Pi[s] = P[s] - ci[s];
          Pj[s] = P[s] - cj[s];
          PQ[s] = P[s] - Q[s];
          setZero(&PQ[s]);
          PQ2 += PQ[s] * PQ[s];
        }
        alpha = p*q / (p+q);
        x = alpha * PQ2;
        factor1 = 2*pow(M_PI, 5.0/2.0);
        factor1 /= p*q * sqrt(p+q);
        /* Hermite coefficient for x2 */
        for(t2=0;t2<lma[0]+1;t2++){
          Hx2 = Hermite(lma[0], 0, t2, q2, nu, 
                        0, 0, 0);
          /* Hermite coefficient for y2 */
          for(u2=0;u2<lma[1]+1;u2++){
            Hy2 = Hermite(lma[1], 0, u2, q2, nu, 
                          0, 0, 0);
            /* Hermite coefficient for z2 */
            for(v2=0;v2<lma[2]+1;v2++){
              Hz2 = Hermite(lma[2], 0, v2, q2, nu, 
                            0, 0, 0);
              /* Hermite coefficient for x1 */
              for(t1=0;t1<lmij[0]+1;t1++){
                Hx1 = Hermite(lmi[0], lmj[0], t1, p2, mu, 
                              cij[0], Pi[0], Pj[0]);
                /* Hermite coefficient for y1 */
                for(u1=0;u1<lmij[1]+1;u1++){
                  Hy1 = Hermite(lmi[1], lmj[1], u1, p2, mu, 
                                cij[1], Pi[1], Pj[1]);
                  /* Hermite coefficient for z1 */
                  for(v1=0;v1<lmij[2]+1;v1++){
                    Hz1 = Hermite(lmi[2], lmj[2], v1, p2, mu, 
                                  cij[2], Pi[2], Pj[2]);
                    factor2 = Nij*Na*Hx1*Hy1*Hz1*Hx2*Hy2*Hz2;
                    factor2 *= pow(-1, t2+u2+v2);
                    /* 2-center integral */
                    HC_cef = HCcef(t1+t2, u1+u2, v1+v2, 
                                   0, PQ, 2*alpha, x);
                    element_aij += factor1*factor2*HC_cef;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  free(expa);
  free(expi);
  free(expj);
  free(cefa);
  free(cefi);
  free(cefj);

  return element_aij;
}

/****************************************************************
*  main function for Gaussian-Coulomb density-fitting integral  *
****************************************************************/
// calculate 2D array for 2-center Gaussian-Coulomb
// electron-electron repulsion matrix from atomic Gaussian density
// It return 2D array, OUT[a, b] where a, b are the basis indices
// for density gaussian expansion
// NOTE: Exchange interaction need to be considered 
//       to compute two-electron energy
// neMatrix(data, center, exp, cef, ng, 
//          lm_xyz, Nao, N, i, j, k, l);
double nnMatrix(double *center,
                double *exp,
                double *cef,
                int *ng,
                int *lm_xyz,
                int Nao,
                int ana,
                int anb
              ){

  /* input related variables */
  int s, a, b, a0 = 0, b0 = 0;
  double ca[3], cb[3];
  int lma[3], lmb[3];
  int nga = ng[ana], ngb = ng[anb];
  double *expa, *expb, *cefa, *cefb;
  expa = (double*) malloc(nga * sizeof(double));
  expb = (double*) malloc(ngb * sizeof(double));
  cefa = (double*) malloc(nga * sizeof(double));
  cefb = (double*) malloc(ngb * sizeof(double));

  /* gaussian integral variables */
  double Q[3];
  double P[3];
  double PQ[3], PQ2, alpha;
  int t1, u1, v1, t2, u2, v2;
  double q, q2, nu, Na;
  double p, p2, mu, Nb;
  double x;
  double Hx1, Hy1, Hz1;
  double Hx2, Hy2, Hz2;
  double factor1, factor2, HC_cef; 
  double element_ab = 0;

  /* setup variables and save to allocated memory */
  for(s=0;s<ana;s++) a0 += ng[s];
  for(s=0;s<anb;s++) b0 += ng[s];
  for(s=0;s<3;s++){
    lma[s] = lm_xyz[s + 3*ana];
    lmb[s] = lm_xyz[s + 3*anb];
    ca[s] = center[s + 3*ana];
    cb[s] = center[s + 3*anb];
  }
  /* loop for a_th gaussian density */
  for(a=0;a<nga;a++){
    expa[a] = exp[a0+a];
    cefa[a] = cef[a0+a];
    for(s=0;s<3;s++) P[s] = ca[s];
    Na = cefa[a] * Norm(expa[a], lma);
    p = expa[a];
    p2 = 2*p;
    mu = p;
    for(b=0;b<ngb;b++){
      expb[b] = exp[b0+b];
      cefb[b] = cef[b0+b];
      for(s=0;s<3;s++) Q[s] = cb[s];
      Nb = cefb[b] * Norm(expb[b], lmb);
      q = expb[b];
      q2 = 2*q;
      nu = q;
      PQ2 = 0;
      for(s=0;s<3;s++){
        PQ[s] = P[s] - Q[s];
        setZero(&PQ[s]);
        PQ2 += PQ[s] * PQ[s];
      }
      alpha = p*q / (p+q);
      x = alpha * PQ2;
      factor1 = 2*pow(M_PI, 5.0/2.0);
      factor1 /= p*q * sqrt(p+q);
      /* Hermite coefficient for x1 */
      for(t1=0;t1<lma[0]+1;t1++){
        Hx1 = Hermite(lma[0], 0, t1, p2, mu,
                      0, 0, 0);
        /* Hermite coefficient for y1 */
        for(u1=0;u1<lma[1]+1;u1++){
          Hy1 = Hermite(lma[1], 0, u1, p2, mu,
                        0, 0, 0);
          /* Hermite coefficient for z1 */
          for(v1=0;v1<lma[2]+1;v1++){
            Hz1 = Hermite(lma[2], 0, v1, p2, mu,
                          0, 0, 0);
            /* Hermite coefficient for x2 */
            for(t2=0;t2<lmb[0]+1;t2++){
              Hx2 = Hermite(lmb[0], 0, t2, q2, nu,
                            0, 0, 0);
              /* Hermite coefficient for y2 */
              for(u2=0;u2<lmb[1]+1;u2++){
                Hy2 = Hermite(lmb[1], 0, u2, q2, nu,
                              0, 0, 0);
                /* Hermite coefficient for z2 */
                for(v2=0;v2<lmb[2]+1;v2++){
                  Hz2 = Hermite(lmb[2], 0, v2, q2, nu,
                                0, 0, 0);
                  factor2 = Na*Nb*Hx1*Hy1*Hz1*Hx2*Hy2*Hz2;
                  factor2 *= pow(-1, t2+u2+v2);
                  /* 2-center integral */
                  HC_cef = HCcef(t1+t2, u1+u2, v1+v2, 
                                 0, PQ, 2*alpha, x);
                  element_ab += factor1*factor2*HC_cef;
                }
              }
            }
          }
        }
      }
    }
  }
  free(expa);
  free(expb);
  free(cefa);
  free(cefb);

  return element_ab;
}

/******************************************************
*  main function for Gaussian 2nd-derivativeintegral  *
******************************************************/
//keMatrix(R, Z, center, exp, cef, ng, 
//         lm_xyz, Nao, N, i, j);
double keMatrix(double *center,
                double *exp,
                double *cef,
                int *ng,
                int *lm_xyz,
                int aoi,
                int aoj){

  double element_ij = 0;
  element_ij += aoOverlap(center, exp, cef, ng, lm_xyz, aoi, aoj, 
                          1, 0, 0, 0, 2, 0, 0);
  element_ij += aoOverlap(center, exp, cef, ng, lm_xyz, aoi, aoj, 
                          1, 0, 0, 0, 0, 0, 0);
  element_ij += aoOverlap(center, exp, cef, ng, lm_xyz, aoi, aoj, 
                          1, 0, 0, 0, -2, 0, 0);
  element_ij += aoOverlap(center, exp, cef, ng, lm_xyz, aoi, aoj, 
                          2, 0, 0, 0, 0, 2, 0);
  element_ij += aoOverlap(center, exp, cef, ng, lm_xyz, aoi, aoj, 
                          2, 0, 0, 0, 0, 0, 0);
  element_ij += aoOverlap(center, exp, cef, ng, lm_xyz, aoi, aoj, 
                          2, 0, 0, 0, 0, -2, 0);
  element_ij += aoOverlap(center, exp, cef, ng, lm_xyz, aoi, aoj, 
                          3, 0, 0, 0, 0, 0, 2);
  element_ij += aoOverlap(center, exp, cef, ng, lm_xyz, aoi, aoj, 
                          3, 0, 0, 0, 0, 0, 0);
  element_ij += aoOverlap(center, exp, cef, ng, lm_xyz, aoi, aoj, 
                          3, 0, 0, 0, 0, 0, -2);
  return -0.5*element_ij;
}

/*******************************************************
*  main function for Gaussian 1st-derivative integral  *
*******************************************************/
//knMatrix(R, Z, center, exp, cef, ng, 
//         lm_xyz, Nao, N, i, j);
double knMatrix(double *R,      //all the rest are input
                double *Z, 
                double *center,
                double *exp,
                double *cef,
                int *ng,
                int *lm_xyz,
                int Nao,
                int N,
                int ana
             ){  

  /* input related variables */
  int s, a, a0 = 0, I;
  double ca[3];
  int lma[3];
  int nga = ng[ana];
  double *expa, *cefa;
  expa = (double*) malloc(nga * sizeof(double));
  cefa = (double*) malloc(nga * sizeof(double));

  /* gaussian integral variables */
  double cI[3], P[3], PI[3], PI2;
  int t, u, v;
  double p, p2, mu, ZI, x, HC_cef, Na;
  double Hx, Hy, Hz, num, element_a = 0;

  for(s=0;s<ana;s++) a0 += ng[s];
  for(s=0;s<3;s++){
    lma[s] = lm_xyz[s + 3*ana];
    ca[s] = center[s + 3*ana];
  }
  for(a=0;a<nga;a++){
    expa[a] = exp[a0+a];
    cefa[a] = cef[a0+a];
    p = expa[a];
    p2 = 2*p;
    mu = p;
    Na = cefa[a] * Norm(expa[a], lma);
    for(s=0;s<3;s++){
      P[s] = ca[s];
    }
    for(I=0;I<N;I++){
      PI2 = 0;
      for(s=0;s<3;s++){
        cI[s] = R[s + 3*I];
        PI[s] = P[s] - cI[s];
        PI2 += PI[s] * PI[s];
      }
      ZI = Z[I];
      x = p*PI2;

      for(t=0;t<lma[0]+1;t++){
        Hx = Hermite(lma[0], 0, t, p2, mu,
                     0, 0, 0);
        for(u=0;u<lma[1]+1;u++){
          Hy = Hermite(lma[1], 0, u, p2, mu,
                       0, 0, 0);
          for(v=0;v<lma[2]+1;v++){
            Hz = Hermite(lma[2], 0, v, p2, mu,
                         0, 0, 0);
            HC_cef = HCcef(t, u, v, 0, PI, p2, x);
            num = ZI*(Na*Hx*Hy*Hz*HC_cef);
            element_a -= num*(2*M_PI/p);
          }
        }
      }
    }
  }

  free(expa);
  free(cefa);
  return element_a;
}
