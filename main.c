#include "csr.h"
#include <malloc.h>
#include <stdlib.h>

void main(){
    sparse_matrix mat;
    mat.m = 4;
    mat.n = 4;
    mat.nz = 7;
    float exMat[][4] = {{11.0,12.0,0.0,0.0},{0.0,22.0,23.0,0.0},{0.0,0.0,33.0,0.0},{0.0,0.0,43.0,44.0}};
    mat.coeff = (float**)malloc(sizeof(float*)*4);
    for(int i = 0; i < 4; i++){
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


}