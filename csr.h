#include "matrix_generator.h"
#include <malloc.h>
#include <stdlib.h>

#ifndef _CRSH_
#define _CSRH_ 

typedef struct{
    int m;              //Numero righe matrice
    int n;              //Numero colonne matrice
    int nz;             //Numero non zeri
    int*        irp;    //Vettore dei puntatori all'inizio di ciascuna riga
    float*      as;     //Vettore dei coefficienti
    int*        ja;     //Vettore degli indici di colonna
} csr_matrix;

void stampaMatriceCsr(csr_matrix mat){
    printf("AS\n");
    for(int i = 0; i < mat.nz; i++){
        printf("%f",mat.as[i]);
        printf(" ");
    }
    printf("\nJS\n");
    for(int i = 0; i < mat.nz; i++){
        printf("%d",mat.ja[i]);
        printf(" ");
    }
    printf("\nIRP\n");
    for(int i = 0; i <= mat.m; i++){
        printf("%d",mat.irp[i]);
        printf(" ");
    }
    printf("\n");
}

csr_matrix convertToCsr(sparse_matrix matrix){
    //Nuova matrice
    csr_matrix convertedMatrix;
    //Alloco i vettori
    convertedMatrix.m = matrix.m;
    convertedMatrix.n = matrix.n;
    convertedMatrix.nz = matrix.nz;
    convertedMatrix.as = (float*)malloc(matrix.nz*sizeof(float));
    if(convertedMatrix.as == NULL){
        exit(1);
    }
    convertedMatrix.ja = (int*)malloc(matrix.nz*sizeof(int));
    if(convertedMatrix.ja == NULL){
        exit(1);
    }
    convertedMatrix.irp = (int*)malloc((matrix.m +1)*sizeof(int));
    if(convertedMatrix.irp == NULL){
        exit(1);
    }
    int nz = 0;
    //Eseguo la conversione
    for(int i = 0; i < convertedMatrix.m; i++){
        convertedMatrix.irp[i] = nz;
        for(int j = 0; j < convertedMatrix.n; j++){
            float elem = matrix.coeff[i][j];
            if(elem != 0){
                convertedMatrix.as[nz] = elem;
                convertedMatrix.ja[nz] = j; 
                nz++;
            }
        }
    }
    convertedMatrix.irp[convertedMatrix.m] = nz;

    return convertedMatrix;
    
}

#endif