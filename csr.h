#include "matrix_generator.h"
#include "test_matrices/matrix_retriever.h"
#include <malloc.h>
#include <stdlib.h>

#ifndef _CRSH_
#define _CSRH_ 

typedef struct{
    int m;              //Numero righe matrice
    int n;              //Numero colonne matrice
    int nz;             //Numero non zeri
    int*        irp;    //Vettore dei puntatori all'inizio di ciascuna riga
    double*      as;    //Vettore dei coefficienti
    int*        ja;     //Vettore degli indici di colonna
} csr_matrix;

void stampaMatriceCsr(csr_matrix mat){
    printf("AS\n");
    for(int i = 0; i < mat.nz; i++){
        printf("%f",mat.as[i]);
        printf(" ");
    }
    printf("\nJA\n");
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

csr_matrix convertToCsr(sparse_matrix mat){
    //Nuova matrice
    csr_matrix convertedMatrix;
    //Alloco i vettori
    convertedMatrix.m = mat.m;
    convertedMatrix.n = mat.n;
    convertedMatrix.nz = mat.nz;
    convertedMatrix.as = (double*)malloc(mat.nz*sizeof(double));
    if(convertedMatrix.as == NULL){
        exit(1);
    }
    convertedMatrix.ja = (int*)malloc(mat.nz*sizeof(int));
    if(convertedMatrix.ja == NULL){
        exit(1);
    }
    convertedMatrix.irp = (int*)malloc((mat.m +1)*sizeof(int));
    if(convertedMatrix.irp == NULL){
        exit(1);
    }
    
    int nz = 0;
    //Eseguo la conversione
    for(int i = 0; i < convertedMatrix.m; i++){
        convertedMatrix.irp[i] = nz;
        for(int j = 0; j < convertedMatrix.n; j++){
            float elem = mat.coeff[i][j];
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

csr_matrix convertToCsrFromCoo(coo_matrix mat){

    int cumsum = 0;
    int temp, row, dest, last = 0;
    csr_matrix convertedMatrix;

    // Preparo la nuova struttura
    convertedMatrix.m = mat.m;
    convertedMatrix.n = mat.n;
    convertedMatrix.nz = mat.nz;
    convertedMatrix.as = (double*)malloc(mat.nz*sizeof(double));
    if(convertedMatrix.as == NULL){
        exit(1);
    }
    convertedMatrix.ja = (int*)malloc(mat.nz*sizeof(int));
    if(convertedMatrix.ja == NULL){
        exit(1);
    }
    convertedMatrix.irp = (int*)malloc((mat.m +1)*sizeof(int));
    if(convertedMatrix.irp == NULL){
        exit(1);
    }
    
    //Converto
    for(int i = 0; i < mat.nz; i++){
        convertedMatrix.irp[mat.rows[i]] = 1;
    }
    
    for(int i = 0; i < mat.m; i++){
        temp = convertedMatrix.irp[i];
        convertedMatrix.irp[i] = cumsum;
        cumsum += temp;
    }
    convertedMatrix.irp[mat.m] = mat.nz;
    
    for(int i = 0; i < mat.nz; i++){
        row = mat.rows[i];
        dest = convertedMatrix.irp[row];
        convertedMatrix.as[dest] = mat.values[i];
        convertedMatrix.ja[dest] = mat.cols[i];

        convertedMatrix.irp[row]++;
    }

    for(int i = 0; i <= mat.m; i++){
        temp = convertedMatrix.irp[i];
        convertedMatrix.irp[i] = last;
        last = temp;
    }
    return convertedMatrix;
}

#endif