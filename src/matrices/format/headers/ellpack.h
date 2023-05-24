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

void print_ellpack_matrix(ellpack_matrix matrix);
ellpack_matrix convert_to_ellpack(sparse_matrix* matrix);
ellpack_matrix convert_coo_to_ellpack(coo_matrix mat);
void free_ellpack_matrix(ellpack_matrix* matrix);

#endif