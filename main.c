#include "csr.h"
#include <malloc.h>
#include <stdlib.h>

void main(){
    sparse_matrix mat;
    mat.m = 3;
    mat.n = 6;
    mat.nz = 7;
    float exMat[][6] = {{2.0,0.0,0.0,5.0,1.0,0.0},{0.0,2.0,3.0,0.0,0.0,0.0},{0.0,0.0,2.0,1.0,0.0,0.0}};
    mat.coeff = (float**)malloc(sizeof(float*)*3);
    for(int i = 0; i < 3; i++){
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