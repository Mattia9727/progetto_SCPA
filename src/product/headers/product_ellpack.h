#ifndef _ELLPACKPRODUCTH_
#define _ELLPACKPRODUCTH_

#include "../../matrices/format/headers/ellpack.h"

#define NUM_THREADS 40


void ellpack_product(ellpack_matrix* mat, matrix* vector, matrix* result);
double omp_ellpack_product(ellpack_matrix mat, matrix vector, matrix* result);

double optimized_omp_ellpack_product(ellpack_matrix mat, matrix vector, matrix* result);

ellpack_matrix* move_ellpack_matrix_on_gpu(ellpack_matrix* mat);
void free_ellpack_matrix_on_gpu(double* AS, int* JA);
matrix* move_matrix_on_gpu(matrix* vector);
void optimized_cuda_ellpack_product(ellpack_matrix* mat, matrix vector, matrix* result);
double optimized_cuda_ellpack_product_in(ellpack_matrix host_mat, matrix vector, matrix* result);

#endif