#include <omp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_erf.h>
#include <math.h>
#include "gaussian.h"

// double checked with NWChem overlap matrix 
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

// double checked with NWChem overlap matrix 
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
  //printf("%10.6E\n", pow(num1 * num2 / N_xyz, 0.5));
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

//  approximated Boy's funciton, 
//  very similar results as gamma function implimentation
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

// double checked with NWChem overlap matrix 
/********************************************
*  Overlap integral of two atomic orbitals  *
********************************************/
// aoverlap = aoOverlap(center, exp, cef, ng, lm_xyz, ao, ao);
double aoOverlap(double *center, double *exp, double *cef, int *ng,
                 int *lm_xyz, int aoi, int aoj){
  int ngi = ng[aoi], ngj = ng[aoj];
  int i, j, k, lmi[3], lmj[3], lmij[3], t, u, v, i0=0, j0=0;
  double p, p2, mu, cef_out;
  double ci[3], cj[3], cij[3];
  double P[3], Pi[3], Pj[3];
  double norm;
  double Hx, Hy, Hz, overlap = 0;

  double *expi = (double*) malloc(ngi * sizeof(double));
  double *expj = (double*) malloc(ngj * sizeof(double));
  double *cefi = (double*) malloc(ngi * sizeof(double));
  double *cefj = (double*) malloc(ngj * sizeof(double));

	for(i=0;i<aoi;i++) i0 += ng[i];
  for(j=0;j<aoj;j++) j0 += ng[j];
  for(k=0;k<3;k++){
    lmi[k] = lm_xyz[k+aoi*3];
    lmj[k] = lm_xyz[k+aoj*3];
    ci[k] = center[k+aoi*3];
    cj[k] = center[k+aoj*3];
    cij[k] = ci[k] - cj[k];
  }

  for(i=0;i<ngi;i++){
    expi[i] = exp[i0+i];
    cefi[i] = cef[i0+i];
    for(j=0;j<ngj;j++){
      expj[j] = exp[j0+j];
      cefj[j] = cef[j0+j];
      cef_out = cefi[i]*cefj[j];
      norm = Norm(expi[i], lmi)*Norm(expj[j], lmj);
      p = expi[i] + expj[j];
      p2 = 2*p;
      mu = expi[i]*expj[j] / p;
      for(k=0;k<3;k++){
        P[k] = (expi[i]*ci[k] + expj[j]*cj[k])/p;
        Pi[k] = P[k] - ci[k];
        Pj[k] = P[k] - cj[k];
      }
      for(t=0;t<lmij[0]+1;t++){
        Hx = Hermite(lmi[0], lmj[0], t, p2, mu, 
                     cij[0], Pi[0], Pj[0]);
        for(u=0;u<lmij[1]+1;u++){
          Hy = Hermite(lmi[1], lmj[1], u, p2, mu, 
                       cij[1], Pi[1], Pj[1]);
          for(v=0;v<lmij[2]+1;v++){
            Hz = Hermite(lmi[2], lmj[2], v, p2, mu, 
                         cij[2], Pi[2], Pj[2]);
            overlap += cef_out * norm * Hx*Hy*Hz
                       * pow(M_PI/p, 1.5);
          }
        }
      }
    }
  }

  free(expi);
  free(expj);
  free(cefi);
  free(cefj);

  return overlap;
}

// double checked with NWChem overlap matrix
/************************************
*  Renormalization of AO expansion  *
************************************/
// renormalize(center, exp, cef, ng, lm_xyz, Nao);
void renormalize(double *center, double *exp, double *cef, 
                 int *ng, int *lm_xyz, int Nao){
  int ao, i, ngo=0;
  double overlap, gocef;

  for(ao=0;ao<Nao;ao++){
    overlap = aoOverlap(center, exp, cef, ng, lm_xyz, ao, ao);
    for(i=ngo;i<ngo+ng[ao];i++){
      cef[i] /= sqrt(overlap);
    }
    ngo += ng[ao];
  }
}

/****************************************
*  Orthogonalization of overlap matrix  *
****************************************/ 
// calculate and diagonalize overlap matrix using LAPACK functions
// at the moment not used
// orthogonalize(overlap, center, exp, cef, ng, lm_xyz, Nao);
void orthogonalize(double *overlap, double *center, 
                   double *exp, double *cef, int *ng, 
                   int *lm_xyz, int Nao){
  int i, j;

  printf("before diagonalization\n");
  for(i=0;i<Nao;i++){
    for(j=0;j<Nao;j++){
      overlap[j+i*Nao] = aoOverlap(center, exp, cef, 
                                   ng, lm_xyz, i, j);
      if(pow(overlap[j+i*Nao],2) > 0){
        printf("S(%d, %d) = %10.8f\n", i+1, j+1, overlap[j+i*Nao]);
      }
    }
  }

  /* LAPACK variables */
  int n = Nao, lda = Nao, info, lwork = -1;
  double wkopt;
  double *work;
  double *w = (double *) malloc(Nao * Nao * sizeof(double));
  dsyev_("Vectors", "Upper", &n, overlap, 
         &lda, w, &wkopt, &lwork, &info);
  lwork = (int)wkopt;
  work = (double*)malloc(lwork*sizeof(double));
  dsyev_("Vectors", "Upper", &n, overlap, 
         &lda, w, work, &lwork, &info);
  /* end of LAPACK routine */

//  printf("after diagonalization\n");
//  for(i=0;i<Nao;i++){
//    for(j=0;j<Nao;j++){
//      if(overlap[j+i*Nao] > 0){
//        printf("S(%d, %d) = %f\n", i, j, overlap[j+i*Nao]);
//      }
//    }
//  }
//  for(i=0;i<Nao;i++) printf(" %f", w[i]);
//  printf("\n");
  free(w);
  free(work);
}

/************************************************
*  main function for Gaussian-Coulomb integral  *
************************************************/
//gc2Matrix(data, center, exp, cef, ng, 
//         lm_xyz, Nao, N, i, j, k, l);
double gc2Matrix(double *center,
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
  double *expi = (double*) malloc(ngi * sizeof(double));
  double *expj = (double*) malloc(ngj * sizeof(double));
  double *expk = (double*) malloc(ngk * sizeof(double));
  double *expl = (double*) malloc(ngl * sizeof(double));
  double *cefi = (double*) malloc(ngi * sizeof(double));
  double *cefj = (double*) malloc(ngj * sizeof(double));
  double *cefk = (double*) malloc(ngk * sizeof(double));
  double *cefl = (double*) malloc(ngl * sizeof(double));

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
  for(i=0;i<ngi;i++){
    expi[i] = exp[i0+i];
    cefi[i] = cef[i0+i];
    Ni = cefi[i] * Norm(expi[i], lmi);
    for(j=0;j<ngj;j++){
      expj[j] = exp[j0+j];
      cefj[j] = cef[j0+j];
      Nj = cefj[j] * Norm(expj[j], lmj);
      Nij = Ni * Nj;
      p = expi[i] + expj[j];
      p2 = 2*p;
      mu = expi[i] * expj[j] / p;
      for(s=0;s<3;s++){
        P[s] = (expi[i]*ci[s] + expj[j]*cj[s]) / p;
        Pi[s] = P[s] - ci[s];
        Pj[s] = P[s] - cj[s];
      }
      for(k=0;k<ngk;k++){
        expk[k] = exp[k0+k];
        cefk[k] = cef[k0+k];
        Nk = cefk[k] * Norm(expk[k], lmk);
        for(l=0;l<ngl;l++){
          expl[l] = exp[l0+l];
          cefl[l] = cef[l0+l];
          Nl = cefl[l] * Norm(expl[l], lml);
          Nkl = Nk * Nl;
          q = expk[k] + expl[l];
          q2 = 2*q;
          nu = expk[k] * expl[l] / p;
          factor1 = 2*pow(M_PI, 5.0/2.0);
          factor1 /= p*q * sqrt(p+q);
          PQ2 = 0;
          for(s=0;s<3;s++){
            Q[s] = (expk[k]*ck[s] + expl[l]*cl[s]) / q;
            Qk[s] = Q[s] - ck[s];
            Ql[s] = Q[s] - cl[s];
            PQ[s] = P[s] - Q[s];
            PQ2 += PQ[s] * PQ[s];
          }
          alpha = p*q / (p+q);
          x = alpha * PQ2;
          for(t1=0;t1<lmij[0]+1;t1++){
            Hx1 = Hermite(lmi[0], lmj[0], t1, p2, mu, 
                          cij[0], Pi[0], Pj[0]);
            for(u1=0;u1<lmij[1]+1;u1++){
              Hy1 = Hermite(lmi[1], lmj[1], u1, p2, mu, 
                            cij[1], Pi[1], Pj[1]);
              for(v1=0;v1<lmij[2]+1;v1++){
                Hz1 = Hermite(lmi[2], lmj[2], v1, p2, mu, 
                              cij[2], Pi[2], Pj[2]);
                for(t2=0;t2<lmkl[0]+1;t2++){
                  Hx2 = Hermite(lmk[0], lml[0], t2, q2, nu, 
                                ckl[0], Qk[0], Ql[0]);
                  for(u2=0;u2<lmkl[1]+1;u2++){
                    Hy2 = Hermite(lmk[1], lml[1], u2, q2, nu, 
                                  ckl[1], Qk[1], Ql[1]);
                    for(v2=0;v2<lmkl[2]+1;v2++){
                      Hz2 = Hermite(lmk[2], lml[2], v2, q2, nu, 
                                    ckl[2], Qk[2], Ql[2]);
                      factor2 = Nij*Nkl*Hx1*Hy1*Hz1*Hx2*Hy2*Hz2;
                      double tmp = Nij*Nkl*Hx1*Hy1*Hz1*Hx2*Hy2*Hz2;
                      factor2 *= pow(-1, t2+u2+v2);
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


/* 2-factorial function */
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
