#ifndef _CRSH_
#define _CSRH_ 

#include "../matrix_generator.h"
#include "coo.h"
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


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
            double elem = mat.coeff[i][j];
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
    //calcola il numero di elementi non zero per ogni riga di A
    memset(convertedMatrix.irp, 0, (mat.m+1) * sizeof(int));

    for (int n = 0; n < mat.nz; n++) {
        convertedMatrix.irp[mat.rows[n]]++;
    }

    //calcola la somma cumulativa degli elementi non zero per riga per ottenere Bp[]
    for (int i = 0, cumsum = 0; i < mat.m; i++) {
        int temp = convertedMatrix.irp[i];
        convertedMatrix.irp[i] = cumsum;
        cumsum += temp;
    }
    convertedMatrix.irp[mat.m] = mat.nz;

    for (int n = 0; n < mat.nz; n++) {
        int row = mat.rows[n];
        int dest = convertedMatrix.irp[row];

        convertedMatrix.ja[dest] = mat.cols[n];
        convertedMatrix.as[dest] = mat.values[n];

        convertedMatrix.irp[row]++;
    }

     for (int i = 0, last = 0; i <= mat.m; i++) {
        int temp = convertedMatrix.irp[i];
        convertedMatrix.irp[i] = last;
        last = temp;
    }
    return convertedMatrix;
}

void free_csr_matrix(csr_matrix* matrix){
    free(matrix->irp);
    free(matrix->as);
    free(matrix->ja);
}

#endif