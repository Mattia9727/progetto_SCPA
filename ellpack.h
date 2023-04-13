#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "matrix_generator.h"

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

multivector ELLPACKProduct(ellpack_matrix* matrix, multivector* vector){
    if (matrix->n != vector->m){
        printf("Impossibile eseguire prodotto matrice multivettore con queste strutture dati.");
        exit(0);
    }
    int m = matrix->m;
    int maxnzr = matrix->maxnz;
    multivector result;
    result.m = matrix->m;
    result.n = vector->n;
    result.coeff = (float **) malloc(sizeof(float *) * result.m);
    printf("MAXNZ: %d",matrix->maxnz);
    fflush(stdout);
    int op = 0;
    for (int i = 0; i < result.m; i++) {
        result.coeff[i] = (float *) calloc(result.n, sizeof(float));
        for (int k = 0; k < result.n; k++) {
            float t = 0;
            for (int j = 0; j < matrix->maxnz; j++) {
                t = t + matrix->AS[i][j]*vector->coeff[matrix->JA[i][j]][k];
                op++;
            }
            result.coeff[i][k] = t;
        }
    }
    printf("\n\nOperazioni effettuate: %d\n\n",op); //Esegue matrix.M * vector.n * MAXNZ operazioni :)
    return result;
}

multivector OmpELLPACKProduct(ellpack_matrix* matrix, multivector* vector){
    if (matrix->n != vector->m){
        printf("Impossibile eseguire prodotto matrice multivettore con queste strutture dati.");
        exit(0);
    }
    int m = matrix->m;
    int maxnzr = matrix->maxnz;
    multivector result;
    result.m = matrix->m;
    result.n = vector->n;
    result.coeff = (float **) malloc(sizeof(float *) * result.m);
    printf("MAXNZ: %d",matrix->maxnz);
    fflush(stdout);
    int op = 0;
    for (int i = 0; i < result.m; i++) {
        result.coeff[i] = (float *) calloc(result.n, sizeof(float));
        for (int k = 0; k < result.n; k++) {
            float t = 0;
            #pragma omp parallel for
            for (int j = 0; j < matrix->maxnz; j++) {
                int tid = omp_get_thread_num();
                printf("Sto parallelizzando con thread %d",tid);
                t = t + matrix->AS[i][j]*vector->coeff[matrix->JA[i][j]][k];
                op++;
            }
            result.coeff[i][k] = t;
        }
    }
    printf("\n\nOperazioni effettuate: %d\n\n",op); //Esegue matrix.M * vector.n * MAXNZ operazioni :)
    return result;
}

multivector OptimizedELLPACKProduct(ellpack_matrix* matrix, multivector* vector){
    if (matrix->n != vector->m){
        printf("Impossibile eseguire prodotto matrice multivettore con queste strutture dati.");
        exit(0);
    }
    int m = matrix->m;
    int maxnzr = matrix->maxnz;
    multivector result;
    result.m = matrix->m;
    result.n = vector->n;
    result.coeff = (float **) malloc(sizeof(float *) * result.m);
    printf("MAXNZ: %d",matrix->maxnz);
    fflush(stdout);
    int op = 0;
    for (int i = 0; i < result.m; i++) {
        result.coeff[i] = (float *) calloc(result.n, sizeof(float));
        for (int k = 0; k < result.n; k++) {
            float t = 0;
            int prev_JA = -1;
            for (int j = 0; j < matrix->maxnz; j++) {
                if (prev_JA>=matrix->JA[i][j]) break;
                prev_JA = matrix->JA[i][j];
                t = t + matrix->AS[i][j]*vector->coeff[prev_JA][k];
                op++;
            }
            result.coeff[i][k] = t;
        }
    }
    printf("\n\nOperazioni effettuate: %d\n\n",op); //Esegue matrix.M * vector.n * MAXNZ operazioni :)
    return result;
}