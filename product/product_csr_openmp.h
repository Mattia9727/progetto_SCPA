
#ifndef _PRODUCTCSROPENMPH_
#define _PRODUCTCSROPENMPH_ 

#include <stdio.h>
#include <omp.h>
#include "product.h"


void calcola_prodotto_per_righe_csr_openmp(csr_matrix csrMatrix, matrix multivector,matrix* result){
    
    if(csrMatrix.n != multivector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }

    float partialSum;

    int i;

    #pragma omp parallel for schedule(static,16) shared(result, csrMatrix, multivector) private(partialSum,i)
    for(i = 0; i < csrMatrix.m; i++){
        for(int k = 0; k < multivector.n; k++){
            partialSum = 0;
            for(int j = csrMatrix.irp[i]; j < csrMatrix.irp[i+1]; j++){
                partialSum += csrMatrix.as[j]*multivector.coeff[csrMatrix.ja[j]][k];
            }
            result->coeff[i][k] = partialSum;
        }
    }    
}  

void calcola_prodotto_per_righe_csr_openmp_ottimizzato(csr_matrix csrMatrix, matrix multivector_trasposto,matrix* result){
    if(csrMatrix.n != multivector_trasposto.n){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }

    float partialSum;

    int i;
    
    #pragma omp parallel for shared(result, csrMatrix, multivector_trasposto) private(partialSum,i)
    for(i = 0; i < csrMatrix.m; i++){
        for(int k = 0; k < multivector_trasposto.m; k++){
            partialSum = 0;
            for(int j = csrMatrix.irp[i]; j < csrMatrix.irp[i+1]; j++){
                partialSum += csrMatrix.as[j]*multivector_trasposto.coeff[k][csrMatrix.ja[j]];
            }
            result->coeff[i][k] = partialSum;
        }
    }

}

/*
void calcola_prodotto_per_righe_csr_openmp_ottimizzato(csr_matrix csrMatrix, matrix multivector_trasposto,matrix* result){
    if(csrMatrix.n != multivector_trasposto.n){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }

    float partialSum;

    int i,k;
    
    #pragma omp parallel for shared(result, csrMatrix, multivector_trasposto) private(partialSum,i,k)
    for(i = 0; i < csrMatrix.m; i++){
        for(k = 0; k < multivector_trasposto.m; k++){
            partialSum = 0;
            for(int j = csrMatrix.irp[i]; j < csrMatrix.irp[i+1]; j++){
                partialSum += csrMatrix.as[j]*multivector_trasposto.coeff[k][csrMatrix.ja[j]];
            }
            result->coeff[i][k] = partialSum;
        }
    }

}
*/

#endif