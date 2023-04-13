#include <malloc.h>
#include <stdlib.h>
#include "product_csr_openmp.h"
#include <time.h>


void main(){

    coo_matrix mat;
    csr_matrix converted_matrix;
    matrix multivector, result;
    clock_t begin, end;
    int col_multivector[7] = {3,4,8,12,16,32,64}; 

    mat = getMatrix("test_matrices/matrix_files/cage4.mtx");
    converted_matrix = convertToCsrFromCoo(mat);
    //stampaMatriceCsr(converted_matrix);

    for(int i = 0; i < 7; i++){
        printf("COLONNE MULTIVETTORE %d\n",col_multivector[i]);
        multivector = GenerateMultivector(mat.n,col_multivector[i]); 
        begin = clock();
        result = prepara_risultato(converted_matrix.m, multivector.n);
        calcola_prodotto_seriale(converted_matrix,multivector, &result);
        end = clock();
        printf("\tserial %f\n",(double)(end - begin) / CLOCKS_PER_SEC);
        // Stampa della matrice
        //stampaMatrice(converted_matrix.m, vector.n,result);
        begin = clock();
        calcola_prodotto_per_righe_csr_openmp(converted_matrix,multivector, &result);
        end = clock();
        printf("\tparall %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
        // Stampa della matrice
        //stampaMatrice(converted_matrix.m, vector.n,result);
 
    }
       
}