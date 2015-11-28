/* Hermite-Gaussian expansion coefficients */
double Hermite(int, int, int, double, double,
               double, double, double);

/* Hermite-Coulomb integral coefficient */
double HCcef(int, int, int, int, 
             double*, double, double);

/* 2-factorial function */
int fac2(int);

/* Gaussian normalization constant */
double Norm(float, int*);

//overlap = aoOverlap(center, exp, cef, ng, lm_xyz, ao, ao);
double aoOverlap(double*, double*, double*, int*,
                 int*, int, int);

// renormalize(center, exp, cef, ng, lm_xyz, Nao);
void renormalize(double*, double*, double*, 
                 int*, int*, int);

// orthogonalize(overlap, center, exp, cef, ng, lm_xyz, Nao);
void orthogonalize(double*, double*, 
                   double*, double*, int*, 
                   int*, int);

/* Boys function */
double F(int, double);
