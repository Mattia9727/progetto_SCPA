#include "csr.h"

float** calcola_prodotto_seriale(csr_matrix csrMatrix, multivector vector){
    
    if(csrMatrix.n != vector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }
    //Converti la matrice in formato csr
    //csr_matrix csrMatrix = convertToCsr(mat);
    float** result = (float**)malloc(sizeof(float*)*csrMatrix.m);
    if(result == NULL){
        printf("Errore malloc\n");
        exit(1);
    }
    for(int i = 0; i < csrMatrix.m; i++){
        result[i] = (float*)malloc(sizeof(float)*vector.n);
    }
    
    float t;
    for(int i = 0; i < csrMatrix.m; i++){
        for(int k = 0; k < vector.n; k++){
            t = 0;
            for(int j = csrMatrix.irp[i]; j < csrMatrix.irp[i+1]; j++){
                t += csrMatrix.as[j]*vector.coeff[csrMatrix.js[j]][k];
            }
            result[i][k] = t;
        }
    }
    return result;
}