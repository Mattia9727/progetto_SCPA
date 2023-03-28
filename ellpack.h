#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "matrix_generator.h"

#define ROWS 10
#define COLS 10

typedef struct {
    int M;
    int N;
    int MAXNZ;
    int** JA;
    float** NZ;
} ellpack_matrix;

void PrintELLPACKMatrix(ellpack_matrix matrix){
    printf("M = %d, N = %d, MAXNZ = %d \n\n",matrix.M, matrix.N, matrix.MAXNZ);
    printf("NZ (Matrice di nonzeri): \n");
    for(int i=0; i < matrix.M; i++) {
        printf("[");
        for (int j = 0; j < matrix.MAXNZ; j++) {
            printf("%f ",matrix.NZ[i][j]);
        }
        printf("]\n");
    }
    printf("\nJA (Indici di nonzeri): \n");
    for(int i=0; i < matrix.M; i++) {
        printf("[");
        for (int j = 0; j < matrix.MAXNZ; j++) {
            printf("%d ",matrix.JA[i][j]);
        }
        printf("]\n");
    }
}

ellpack_matrix ConvertToELLPACK(sparse_matrix* matrix){
    ellpack_matrix new_matrix;
    new_matrix.M = matrix->m;
    new_matrix.N = matrix->n;
    new_matrix.MAXNZ = 0;
    for(int i=0; i < matrix->m; i++){
        int maxnz=0;
        for(int j=0; j < matrix->n; j++){
            if (matrix->coeff[i][j] != 0.0) maxnz++;
        }
        if (new_matrix.MAXNZ < maxnz) new_matrix.MAXNZ = maxnz;
    }
    new_matrix.NZ = malloc(sizeof(float) * new_matrix.MAXNZ * new_matrix.M);
    new_matrix.JA = malloc(sizeof(int) * new_matrix.MAXNZ * new_matrix.M);
    for(int i=0; i < matrix->m; i++){
        int k=0;
        new_matrix.NZ[i] = (int *) calloc(new_matrix.MAXNZ, sizeof(float));
        new_matrix.JA[i] = (int *) calloc(new_matrix.MAXNZ, sizeof(int));
        for(int j=0; j < new_matrix.N; j++){
            if (matrix->coeff[i][j] != 0) {
                new_matrix.NZ[i][k] = matrix->coeff[i][j];
                new_matrix.JA[i][k] = j;
                k++;
            }
        }
    }
    return new_matrix;
}