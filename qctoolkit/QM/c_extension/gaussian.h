// set small numbers to zero for numerical stability
void setZero(double*);

// leading coefficient of gaussian deriative
void setD2Cef(int, int, double, double*);

// Hermite-Gaussian expansion coefficients
double Hermite(int, int, int, double, double,
               double, double, double);

// Hermite-Coulomb integral coefficient
double HCcef(int, int, int, int, 
             double*, double, double);

// 2-factorial function
int fac2(int);

// Gaussian normalization constant
double Norm(float, int*);

// overlap = aoOverlap(center, exp, cef, ng, lm_xyz, ao, ao);
double aoOverlap(double*, double*, double*, int*, int*, 
                 int, int, int, int, int, int, int, int, int);

// nint = densityIntegral(center, exp, cef, ng, lm_xyz, ao, ao);
double densityIntegral(double*, double*, double*, int*,
                       int*, int);

// renormalize(center, exp, cef, ng, lm_xyz, Nao);
void renormalize(double*, double*, double*, 
                 int*, int*, int);

// densityRenormalize(center, exp, cef, ng, lm_xyz, Nao);
void densityRenormalize(double*, double*, double*, 
                        int*, int*, int);

// orthogonalize(overlap, center, exp, cef, ng, lm_xyz, Nao);
void orthogonalize(double*, double*, 
                   double*, double*, int*, 
                   int*, int);

// Boys function
double F(int, double);

// one-center electorn-nucleus integral
double veMatrix(double*, double*, double*, double*, double*,
                int*, int*, int, int, int, int);

void eeKernel(double*, double*, double*, double*, double*, double*,
              int*, int*, int, int, int, int);

// 4D two-center electorn-repulsion integral
double eeMatrix(double*, double*, double*, int*, int*,
                int, int, int, int, int);

// 3D two-center electorn-repulsion integral
double neMatrix(double*, double*, double*, int*, int*, int, 
                double*, double*, double*, int*, int*, int, 
                int, int, int);

// 2D two-center electorn-repulsion integral
double nnMatrix(double*, double*, double*, int*, int*,
                int, int, int);

// one-center electorn-nucleus integral
double vnMatrix(double*, double*, double*, double*, double*,
                int*, int*, int, int, int);

// electron wavefunction second order derivative integral
double keMatrix(double*, double*, double*,
                int*, int*, int, int);

// electron density first order derivative integral
double knMatrix(double*, double*, double*, double*, double*,
                int*, int*, int, int, int);
