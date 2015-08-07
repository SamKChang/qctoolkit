class Descriptor{
    int size;
    double* vector;
    char* type;
  public:
    Descriptor(int size, char name[30]);
    double* getVector();
};

class CoulombMatrix{
    int base_dim;
    double *matrix;
  public:
    CoulombMatrix(int size); //constructor
    double* getVector();
};
