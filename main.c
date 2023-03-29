#include <malloc.h>
#include <stdlib.h>
#include "product.h"

void main(){
    sparse_matrix mat;
    mat.m = 2;
    mat.n = 3;
    mat.nz = 3;
    float exMat[][3] = {{1.0,1.0,0.0},{0.0,0.0,1.0}};
    mat.coeff = (float**)malloc(sizeof(float*)*mat.m);
    for(int i = 0; i < mat.m; i++){
        mat.coeff[i] = &(exMat[i][0]);
    }

    csr_matrix converted_matrix = convertToCsr(mat);
    printf("AS\n");
    for(int i = 0; i < mat.nz; i++){
        printf("%f",converted_matrix.as[i]);
        printf(" ");
    }
    printf("\nJS\n");
    for(int i = 0; i < mat.nz; i++){
        printf("%d",converted_matrix.js[i]);
        printf(" ");
    }
    printf("\nIRP\n");
    for(int i = 0; i <= mat.m; i++){
        printf("%d",converted_matrix.irp[i]);
        printf(" ");
    }
    printf("\n");

    multivector vector;
    vector.m = 3;
    vector.n = 2;
    float exVec[][2] = {{1.0,1.0},{1.0,1.0},{1.0,1.0}};
    vector.coeff = (float**)malloc(sizeof(float*)*vector.m);
    for(int i = 0; i < vector.m; i++){
        vector.coeff[i] = &(exVec[i][0]);
    } 
    printf("sto per calcolare\n");
    float** result = calcola_prodotto_seriale(converted_matrix,vector);

    // Stampa della matrice
    for (int i = 0; i < mat.m; i++) {
        for (int j = 0; j < vector.n; j++) {
            printf("%f ", result[i][j]);
        }
        printf("\n");
    }

}