#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "headers/ellpack.h"

void print_ellpack_matrix(ellpack_matrix matrix){
    printf("M = %d, N = %d, MAXNZ = %d \n\n",matrix.m, matrix.n, matrix.maxnz);
    printf("AS (Matrice di nonzeri): \n");
    for(int i=0; i < matrix.m; i++) {
        printf("[");
        for (int j = 0; j < matrix.maxnz; j++) {
            printf("%f ",matrix.AS[i * matrix.maxnz + j]);
        }
        printf("]\n");
    }
    printf("\nJA (Indici di nonzeri): \n");
    for(int i=0; i < matrix.m; i++) {
        printf("[");
        for (int j = 0; j < matrix.maxnz; j++) {
            printf("%d ",matrix.JA[i * matrix.maxnz + j]);
        }
        printf("]\n");
    }
    printf("\n\n");
}

ellpack_matrix convert_to_ellpack(sparse_matrix* matrix){
    ellpack_matrix new_matrix;
    new_matrix.m = matrix->m;
    new_matrix.n = matrix->n;
    new_matrix.maxnz = 0;
    for(int i=0; i < matrix->m; i++){
        int maxnz=0;
        for(int j=0; j < matrix->n; j++){
            if (matrix->coeff[i][j] != 0.0) maxnz++;
        }
        if (new_matrix.maxnz < maxnz) new_matrix.maxnz = maxnz;
    }
    new_matrix.AS = (double*)calloc(new_matrix.maxnz * new_matrix.m, sizeof(double));
    new_matrix.JA = (int*)calloc(new_matrix.maxnz * new_matrix.m,sizeof(int));
    for(int i=0; i < matrix->m; i++){
        int k=0;
        for(int j=0; j < new_matrix.n; j++){
            if (matrix->coeff[i][j] != 0) {
                new_matrix.AS[i*new_matrix.maxnz+k] = matrix->coeff[i][j];
                new_matrix.JA[i*new_matrix.maxnz+k] = j;
                k++;
            }
        }
    }
    return new_matrix;
}

ellpack_matrix convert_coo_to_ellpack(coo_matrix mat){
    ellpack_matrix converted_matrix;
    converted_matrix.m = mat.m;
    converted_matrix.n = mat.n;
    converted_matrix.maxnz = 0;
    int* row_arr = (int*)calloc(mat.m,sizeof(int));
    for (int i=0; i<mat.nz; i++){
        row_arr[mat.rows[i]]++;
    }
    for (int i=0; i<mat.m; i++){
        if (converted_matrix.maxnz < row_arr[i]) converted_matrix.maxnz = row_arr[i];
    }
    free(row_arr);
    unsigned long mat_dim = (unsigned long)converted_matrix.m*(unsigned long)converted_matrix.maxnz;
    converted_matrix.AS = (double*)calloc(mat_dim,sizeof(double));
    converted_matrix.JA = (int*)calloc(mat_dim,sizeof(int));
    int row;
    int* col_arr = (int*)calloc(converted_matrix.m,sizeof(int));
    for (int i=0; i<mat.nz; i++){
        row = mat.rows[i];
        converted_matrix.AS[(unsigned long)row*converted_matrix.maxnz +col_arr[row]] = (double)mat.values[i];
        converted_matrix.JA[(unsigned long)row*converted_matrix.maxnz +col_arr[row]] = mat.cols[i];
        col_arr[row]++;
    }
    free(col_arr);

    return converted_matrix;
}

void free_ellpack_matrix(ellpack_matrix* matrix){
    free(matrix->AS);
    free(matrix->JA);
}
