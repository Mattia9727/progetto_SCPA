#ifndef _ELLPACKPRODUCTH_
#define _ELLPACKPRODUCTH_

#include "../../matrices/format/headers/ellpack.h"

#define NUM_THREADS 40


void ellpack_product(ellpack_matrix* mat, matrix* vector, matrix* result);
double omp_ellpack_product(ellpack_matrix mat, matrix vector, matrix* result);

double optimized_ellpack_product(ellpack_matrix mat, matrix vector, matrix* result, int nThreads);

#endif