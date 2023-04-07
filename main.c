#include <malloc.h>
#include <stdlib.h>
#include "product_csr_openmp.h"
#include <time.h>

void main(){

    sparse_matrix mat = GenerateSparseMatrix(20,20,5);

    csr_matrix converted_matrix = convertToCsr(mat);
    
    //stampaMatriceCsr(converted_matrix);

    matrix multivector = GenerateMultivector(20,10); 
    clock_t begin = clock();
    matrix result = prepara_risultato(converted_matrix.m, multivector.n);
    calcola_prodotto_seriale(converted_matrix,multivector, &result);
    clock_t end = clock();
    printf("TEMPO SERIALE %f\n",(double)(end - begin) / CLOCKS_PER_SEC);
    // Stampa della matrice
    //stampaMatrice(converted_matrix.m, vector.n,result);
    begin = clock();
    calcola_prodotto_per_righe_csr_openmp(converted_matrix,multivector, &result);
    end = clock();
    printf("TEMPO PARALLELO %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
    // Stampa della matrice
    //stampaMatrice(converted_matrix.m, vector.n,result);

}