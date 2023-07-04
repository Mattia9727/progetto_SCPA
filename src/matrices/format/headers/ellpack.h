#ifndef _ELLPACKH_
#define _ELLPACKH_

#include "coo.h"

#define ROWS 10
#define COLS 10

typedef struct {
    int m;
    int n;
    int maxnz;
    int* JA;
    double* AS;
} ellpack_matrix;

typedef struct {
    long m;
    long n;
    long hackSize;
    long numMatrix;
    long matDim;
    long* hackOffsets;
    long* maxnz;
    long* JA;
    double* AS;
} h_ellpack_matrix;

typedef struct {
    int m;
    int n;
    int hackSize;
    int numMatrix;
    int matDim;
    int* hackOffsets;
    int* maxnz;
    int* JA;
    double* AS;
} h_ellpack_matrix_bis;

void print_ellpack_matrix(ellpack_matrix matrix);
ellpack_matrix convert_to_ellpack(sparse_matrix* matrix);
ellpack_matrix convert_coo_to_ellpack(coo_matrix mat);
void free_ellpack_matrix(ellpack_matrix* matrix);

void print_h_ellpack_matrix(h_ellpack_matrix matrix);
h_ellpack_matrix convert_coo_to_h_ellpack(coo_matrix mat);
h_ellpack_matrix_bis convert_coo_to_h_ellpack_bis(coo_matrix mat);
void free_h_ellpack_matrix(h_ellpack_matrix* matrix);
void free_h_ellpack_matrix(h_ellpack_matrix_bis* matrix);

void fprint_ellpack_matrix(ellpack_matrix matrix);
void fprint_h_ellpack_matrix(h_ellpack_matrix matrix);
void fprint_h_ellpack_matrix_bis(h_ellpack_matrix_bis matrix);

#endif