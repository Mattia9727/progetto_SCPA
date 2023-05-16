#ifndef _PRODUCTH_
#define _PRODUCTH_

#include "../matrices/format/csr.h"
#include "../matrices/matrix_generator.h"

#define ALPHA 10000

matrix prepara_risultato(int m, int n){
    matrix result;
    result.m = m;
    result.n = n;

    double** coeff= (double**)malloc(sizeof(double*)*m);
    if(coeff == NULL){
        printf("Errore malloc\n");
        exit(1);
    }
    for(int i = 0; i < m; i++){
        coeff[i] = (double*)malloc(sizeof(double)*n);
    }
    result.coeff = coeff;
    return result;
}

void calcola_prodotto_seriale(csr_matrix csrMatrix, matrix vector, matrix* result){
    
    if(csrMatrix.n != vector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }
    
    double t;
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

void checkResult(matrix m1, matrix m2){
    for(int i = 0; i < m1.m; i++){
        for(int j = 0; j < m1.n; j++){
            if((int)(m1.coeff[i][j]*ALPHA) != (int)(m2.coeff[i][j]*ALPHA))
                exit(-1);
        }
    }
}

#endif