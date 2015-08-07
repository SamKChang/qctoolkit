#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

typedef void* kernelpt;

EXTERNC kernelpt kernel_create(char*, double *input);
EXTERNC void kernel_free(char*, kernelpt);
EXTERNC double kernel_evaluate(char*, kernelpt, 
                               double*, double*,int);

#undef EXTERNC
