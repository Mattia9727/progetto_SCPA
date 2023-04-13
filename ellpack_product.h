#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "ellpack.h"

multivector calcola_prodotto_seriale_ellpack(ellpack_matrix ellpackMatrix, multivector vector, bool opt){

    if(ellpackMatrix.n != vector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }

    multivector result;
    if (opt == true) {result = OptimizedELLPACKProduct(&ellpackMatrix, &vector);}
    else {result = ELLPACKProduct(&ellpackMatrix, &vector);}

    PrintMatrix(result);
    return result;
}

multivector calcola_prodotto_omp_ellpack(ellpack_matrix ellpackMatrix, multivector vector, bool opt){

    if(ellpackMatrix.n != vector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }

    multivector result;
    result = OmpELLPACKProduct(&ellpackMatrix, &vector);

    PrintMatrix(result);
    return result;
}