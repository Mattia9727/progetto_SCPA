#include "csr.h"

#ifndef _PRODUCTH_
#define _PRODUCTH_

matrix prepara_risultato(int m, int n){
    matrix result;
    result.m = m;
    result.n = n;

    float** coeff= (float**)malloc(sizeof(float*)*m);
    if(coeff == NULL){
        printf("Errore malloc\n");
        exit(1);
    }
    for(int i = 0; i < m; i++){
        coeff[i] = (float*)malloc(sizeof(float)*n);
    }
    result.coeff = coeff;
    return result;
}

void calcola_prodotto_seriale(csr_matrix csrMatrix, matrix vector, matrix* result){
    
    if(csrMatrix.n != vector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }
    
    float t;
    //ciclo nelle righe della matrice
    for(int i = 0; i < csrMatrix.m; i++){
        //ciclo nelle colonne del vettore
        for(int k = 0; k < vector.n; k++){
            t = 0;
            // ciclo per somma elementi riga matrice - elementi colonna vettore
            for(int j = csrMatrix.irp[i]; j < csrMatrix.irp[i+1]; j++){
                t += csrMatrix.as[j]*vector.coeff[csrMatrix.ja[j]][k];
            }
            result->coeff[i][k] = t;
        }
    }
}

#endif