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
    int m;
    int n;
    int hackSize;
    int numMatrix;
    int* hackOffsets;
    int* maxnz;
    int* JA;
    double* AS;
} h_ellpack_matrix;

void print_ellpack_matrix(ellpack_matrix matrix);
ellpack_matrix convert_to_ellpack(sparse_matrix* matrix);
ellpack_matrix convert_coo_to_ellpack(coo_matrix mat);
void free_ellpack_matrix(ellpack_matrix* matrix);

void print_h_ellpack_matrix(h_ellpack_matrix matrix);
h_ellpack_matrix convert_coo_to_h_ellpack(coo_matrix mat);
void free_h_ellpack_matrix(h_ellpack_matrix* matrix);

void fprint_ellpack_matrix(ellpack_matrix matrix);
void fprint_h_ellpack_matrix(h_ellpack_matrix matrix);

#endif