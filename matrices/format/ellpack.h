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
    float** AS;
} ellpack_matrix;

void PrintELLPACKMatrix(ellpack_matrix matrix){
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

ellpack_matrix ConvertToELLPACK(sparse_matrix* matrix){
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
    new_matrix.AS = (float**)malloc(sizeof(float) * new_matrix.maxnz * new_matrix.m);
    new_matrix.JA = (int**)malloc(sizeof(int) * new_matrix.maxnz * new_matrix.m);
    for(int i=0; i < matrix->m; i++){
        int k=0;
        new_matrix.AS[i] = (float *) calloc(new_matrix.maxnz, sizeof(float));
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

ellpack_matrix ConvertCOOToELLPACK(coo_matrix mat){
    ellpack_matrix converted_matrix;
    converted_matrix.m = mat.m;
    converted_matrix.n = mat.n;
    converted_matrix.maxnz = 0;
    int* row_arr;
    row_arr = malloc(sizeof(int)*mat.m);
    for (int i=0; i<mat.m; i++) {row_arr[i]=0;}
    for (int i=0; i<mat.nz; i++){
        row_arr[mat.rows[i]]++;
    }
    for (int i=0; i<mat.m; i++){
        if (converted_matrix.maxnz < row_arr[mat.rows[i]]) converted_matrix.maxnz = row_arr[mat.rows[i]];
    }
    free(row_arr);
    converted_matrix.AS = (float**)calloc(converted_matrix.m * converted_matrix.maxnz,sizeof(float));
    converted_matrix.JA = (int**)calloc(converted_matrix.m * converted_matrix.maxnz,sizeof(int));

    for (int i=0; i<mat.m; i++){
        converted_matrix.AS[i] = (float *) calloc(converted_matrix.maxnz, sizeof(float));
        converted_matrix.JA[i] = (int *) calloc(converted_matrix.maxnz, sizeof(int));
    }

    int* col_arr;
    col_arr = malloc(sizeof(int)*converted_matrix.m);
    for (int i=0; i<converted_matrix.m; i++) {col_arr[i]=0;}
    for (int i=0; i<mat.nz; i++){
        converted_matrix.AS[mat.rows[i]][col_arr[mat.rows[i]]] = (float)mat.values[i];
        converted_matrix.JA[mat.rows[i]][col_arr[mat.rows[i]]] = mat.cols[i];
        col_arr[mat.rows[i]]++;
    }
    free(col_arr);


    //PrintELLPACKMatrix(converted_matrix);
    //PrintELLPACKMatrix(converted_matrix);

    return converted_matrix;
}

void ELLPACKProduct(ellpack_matrix* mat, matrix* vector, matrix* result){

    //int op = 0;
    for (int i = 0; i < result->m; i++) {
        result->coeff[i] = (float *) calloc(result->n, sizeof(float));
        for (int k = 0; k < result->n; k++) {
            float t = 0;
            for (int j = 0; j < mat->maxnz; j++) {
                t = t + mat->AS[i][j]*vector->coeff[mat->JA[i][j]][k];
                //op++;
            }
            result->coeff[i][k] = t;
        }
    }
    //printf("\n\nOperazioni effettuate: %d\n\n",op); //Esegue matrix->M * vector->n * MAXNZ operazioni :)
}

void OmpELLPACKProduct(ellpack_matrix* mat, matrix* vector, matrix* result){

    //int op = 0;

    float t;

    int i;

    #pragma omp parallel for shared(result, mat, vector) private(t,i)
    for (i = 0; i < result->m; i++) {
        result->coeff[i] = (float *) calloc(result->n, sizeof(float));
        for (int k = 0; k < result->n; k++) {
            t = 0;
            for (int j = 0; j < mat->maxnz; j++) {
                t = t + mat->AS[i][j]*vector->coeff[mat->JA[i][j]][k];
                //op++;
            }
            result->coeff[i][k] = t;
        }
    }
    //printf("\n\nOperazioni effettuate: %d\n\n",op); //Esegue matrix->M * vector->n * MAXNZ operazioni :)
}

void OptimizedELLPACKProduct(ellpack_matrix* mat, matrix* vector, matrix* result){

    int op = 0;
    for (int i = 0; i < result->m; i++) {
        result->coeff[i] = (float *) calloc(result->n, sizeof(float));
        for (int k = 0; k < result->n; k++) {
            float t = 0;
            int prev_JA = -1;
            for (int j = 0; j < mat->maxnz; j++) {
                if (prev_JA>=mat->JA[i][j]) break;
                prev_JA = mat->JA[i][j];
                t = t + mat->AS[i][j]*vector->coeff[prev_JA][k];
                op++;
            }
            result->coeff[i][k] = t;
        }
    }
    printf("\n\nOperazioni effettuate: %d\n\n",op); //Esegue matrix->M * vector->n * MAXNZ operazioni :)
}

#endif