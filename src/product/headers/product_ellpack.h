#ifndef _ELLPACKPRODUCTH_
#define _ELLPACKPRODUCTH_

#include "../../matrices/format/headers/ellpack.h"

#define NUM_THREADS 40

//product_ellpack_omp.c
void ellpack_product(ellpack_matrix* mat, matrix* vector, matrix* result);
double omp_ellpack_product(ellpack_matrix mat, matrix vector, matrix* result);
double optimized_omp_ellpack_product(ellpack_matrix mat, matrix vector, matrix* result, int nThreads);

//product_ellpack_gpu.c
double* move_AS_on_gpu(ellpack_matrix* mat);
long* move_JA_on_gpu(ellpack_matrix* mat);

double* h_move_AS_on_gpu(h_ellpack_matrix* mat);
long* h_move_JA_on_gpu(h_ellpack_matrix* mat);
long* h_move_maxnz_on_gpu(h_ellpack_matrix* mat);
long* h_move_hackOffsets_on_gpu(h_ellpack_matrix* mat);

double* move_matrix_coeff_on_gpu(matrix* vector);

void optimized_cuda_ellpack_product(ellpack_matrix* mat, matrix vector, matrix* result);
void optimized_cuda_h_ellpack_product(int m, int n, long* maxnz, double* AS, long* JA, long* hackOffsets, long hackSize, long numMatrix, double* coeff, double* myRes);

double optimized_cuda_ellpack_product_in(ellpack_matrix host_mat, matrix vector, matrix* result);

void print_gpu_matrix_d(double* matrix, int m, int n);
void print_gpu_matrix_i(int* matrix, int m, int n);
double optimized_cuda_ellpack_product_in(ellpack_matrix host_mat, matrix vector, matrix* result);
double optimized_cuda_h_ellpack_product_in(h_ellpack_matrix host_mat, matrix vector, matrix* result);



#endif