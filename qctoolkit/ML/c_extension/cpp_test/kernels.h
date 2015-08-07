// interface of evaluate functionality
// inconsistent implementation of evaluate gives error
// for the performance/multithreading purpose,
// python implementation is NOT ideal
// kernels.h, kernels.cpp and pyinterface.cpp must be modified
// when new kernel is added.
class Kernel{
  public:
    virtual double evaluate(double* vec1, double* vec2, 
                            int size) = 0;
};

class Gaussian : public Kernel{
    double sigma;
  public:
    Gaussian(double input);
    double evaluate(double* vec1, double* vec2, int size);
};
