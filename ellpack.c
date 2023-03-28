//
// Created by mattia971 on 3/28/23.
//
#include "stdio.h"
#include "stdlib.h"
#include <time.h>

#define ROWS 10
#define COLS 10

//PROVVISORIO
typedef struct {
    int m;
    int n;
    int nz;
    float** coeff;
}sparse_matrix;

sparse_matrix GenerateMatrix() {
    sparse_matrix new_matrix;
    new_matrix.m = ROWS;
    new_matrix.n = COLS;
    new_matrix.nz = 0;

    new_matrix.coeff = (float **)malloc(sizeof(float*)*ROWS);
    int max_nz = 5;         // numero massimo di elementi non nulli per riga
    int row_nz;             // numero di elementi non nulli nella riga corrente
    int i, j, k;

    // Allocazione della memoria per il doppio puntatore
    new_matrix.coeff = (float **) malloc(sizeof(float *) * ROWS);
    if (new_matrix.coeff == NULL) {
        printf("Errore di allocazione della memoria.\n");
        exit(0);
    }

    // Generazione casuale dei valori della matrice
    srand(time(NULL));
    for (i = 0; i < ROWS; i++) {
        row_nz = rand() % (max_nz + 1);     // numero casuale di elementi non nulli nella riga
        new_matrix.coeff[i] = (float *) calloc(COLS, sizeof(float));     // allocazione della memoria per la riga i-esima
        if (new_matrix.coeff[i] == NULL) {
            printf("Errore di allocazione della memoria.\n");
            exit(0);
        }
        for (j = 0; j < row_nz; j++) {
            k = rand() % COLS;      // indice casuale di colonna per il valore non nullo
            new_matrix.coeff[i][k] = rand() % 10 + 1;     // valore casuale tra 1 e 10
        }
    }

    // Stampa della matrice
    for (i = 0; i < ROWS; i++) {
        for (j = 0; j < COLS; j++) {
            printf("%f ", new_matrix.coeff[i][j]);
        }
        printf("\n");
    }

    for(int i=0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            if (new_matrix.coeff[i][j] !=0) new_matrix.nz++;
        }
    }

    return new_matrix;
}


//PROVVISORIO

typedef struct {
    int M;
    int N;
    int MAXNZ;
    int** JA;
    float** NZ;
} ellpack_matrix;

void PrintELLPACKMatrix(ellpack_matrix matrix){
    printf("M = %d, N = %d, MAXNZ = %d ",matrix.M, matrix.N, matrix.MAXNZ);
    for(int i=0; i < matrix.M; i++) {
        printf("[");
        for (int j = 0; j < matrix.MAXNZ; j++) {
            printf("%f ",matrix.NZ[i][j]);
        }
        printf("]\n");
    }
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

int main(int argc, char** argv){
    sparse_matrix matrix = GenerateMatrix();
    printf("OK\n\n");
    /*for(int i=0; i<ROWS; i++){
        for(int j=0; j<COLS; j++){
            printf("%f ",matrix.coeff[i][j]);
        }
        printf("\n");
    }*/
    printf("%d",matrix.m);
    printf("%d",matrix.n);
    fflush(stdout);
    ellpack_matrix ellpackMatrix = ConvertToELLPACK(&matrix);
    PrintELLPACKMatrix(ellpackMatrix);
    return 1;
}
