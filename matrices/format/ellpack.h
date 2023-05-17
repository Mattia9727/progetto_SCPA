#ifndef _ELLPACKH_
#define _ELLPACKH_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "../matrix_generator.h"
#include "../matrix_retriever.h"

#define ROWS 10
#define COLS 10

typedef struct {
    int m;
    int n;
    int maxnz;
    int** JA;
    double** AS;
} ellpack_matrix;

void print_ellpack_matrix(ellpack_matrix matrix){
    printf("M = %d, N = %d, MAXNZ = %d \n\n",matrix.m, matrix.n, matrix.maxnz);
    printf("AS (Matrice di nonzeri): \n");
    for(int i=0; i < matrix.m; i++) {
        printf("[");
        for (int j = 0; j < matrix.maxnz; j++) {
            printf("%f ",matrix.AS[i][j]);
        }
        printf("]\n");
    }
    printf("\nJA (Indici di nonzeri): \n");
    for(int i=0; i < matrix.m; i++) {
        printf("[");
        for (int j = 0; j < matrix.maxnz; j++) {
            printf("%d ",matrix.JA[i][j]);
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
    new_matrix.AS = (double**)malloc(sizeof(double) * new_matrix.maxnz * new_matrix.m);
    new_matrix.JA = (int**)malloc(sizeof(int) * new_matrix.maxnz * new_matrix.m);
    for(int i=0; i < matrix->m; i++){
        int k=0;
        new_matrix.AS[i] = (double *) calloc(new_matrix.maxnz, sizeof(double));
        new_matrix.JA[i] = (int *) calloc(new_matrix.maxnz, sizeof(int));
        for(int j=0; j < new_matrix.n; j++){
            if (matrix->coeff[i][j] != 0) {
                new_matrix.AS[i][k] = matrix->coeff[i][j];
                new_matrix.JA[i][k] = j;
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
    int* row_arr;
    row_arr = (int*)calloc(mat.m,sizeof(int));
    for (int i=0; i<mat.nz; i++){
        row_arr[mat.rows[i]]++;
    }
    for (int i=0; i<mat.m; i++){
        if (converted_matrix.maxnz < row_arr[i]) converted_matrix.maxnz = row_arr[i];
    }
    free(row_arr);
    converted_matrix.AS = (double**)calloc(converted_matrix.m * converted_matrix.maxnz,sizeof(double*));
    converted_matrix.JA = (int**)calloc(converted_matrix.m * converted_matrix.maxnz,sizeof(int*));

    for (int i=0; i<mat.m; i++){
        converted_matrix.AS[i] = (double *) calloc(converted_matrix.maxnz, sizeof(double));
        converted_matrix.JA[i] = (int *) calloc(converted_matrix.maxnz, sizeof(int));
    }

    int* col_arr;
    col_arr = (int*)calloc(converted_matrix.m,sizeof(int));
    for (int i=0; i<mat.nz; i++){
        converted_matrix.AS[mat.rows[i]][col_arr[mat.rows[i]]] = (double)mat.values[i];
        converted_matrix.JA[mat.rows[i]][col_arr[mat.rows[i]]] = mat.cols[i];
        col_arr[mat.rows[i]]++;
    }
    free(col_arr);


    //print_ellpack_matrix(converted_matrix);

    return converted_matrix;
}

void free_ellpack_matrix(ellpack_matrix* matrix){
    for (int i=0; i<matrix->maxnz; i++){
        free(matrix->AS[i]);
        free(matrix->JA[i]);
    }
    free(matrix->AS);
    free(matrix->JA);
}

#endif