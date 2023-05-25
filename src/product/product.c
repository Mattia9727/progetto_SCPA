#include <time.h>
#include "headers/product.h"

void prepara_risultato(int m, int n, matrix* result){
    result->m = m;
    result->n = n;
    double * coeff= (double*)calloc((unsigned long)m*n,sizeof(double));
    if(coeff == NULL){
        printf("Errore malloc\n");
        exit(1);
    }
    result->coeff = coeff;
}

void free_matrix(matrix* result){
    free(result->coeff);
}

double calcola_prodotto_seriale(csr_matrix csrMatrix, matrix vector, matrix* result){
    
    if(csrMatrix.n != vector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }
    
    double t;
    int irp_1, irp_2;
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    //ciclo nelle righe della matrice
    for(int i = 0; i < csrMatrix.m; i++){
        irp_1 = csrMatrix.irp[i];
        irp_2 = csrMatrix.irp[i+1];
        //ciclo nelle colonne del vettore
        
        for(int k = 0; k < vector.n; k++){
            
            t = 0;
            // ciclo per somma elementi riga matrice - elementi colonna vettore
            for(int j = irp_1; j < irp_2; j++){
                t += csrMatrix.as[j]*vector.coeff[csrMatrix.ja[j]*vector.n + k];
            }
            result->coeff[i * vector.n + k] = t;
            
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    return (double)( end.tv_sec - start.tv_sec )+ ( end.tv_nsec - start.tv_nsec )/ (double)1000000000L;

}

void check_result(matrix m1, matrix m2){
    for(int i = 0; i < m1.m; i++){
        for(int j = 0; j < m1.n; j++){
            if((int)(m1.coeff[i*m1.n + j]*ALPHA) != (int)(m2.coeff[i*m2.n +j]*ALPHA))
                exit(-1);
        }
    }
}
