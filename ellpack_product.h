#include <stdbool.h>
#include "ellpack.h"

matrix calcola_prodotto_seriale_ellpack(ellpack_matrix ellpackMatrix, matrix vector, bool opt){

    if(ellpackMatrix.n != vector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }

    matrix result;
    if (opt == true) {result = OptimizedELLPACKProduct(&ellpackMatrix, &vector);}
    else {result = ELLPACKProduct(&ellpackMatrix, &vector);}

    stampaMatrice(result);
    return result;
}

/*
matrix calcola_prodotto_omp_ellpack(ellpack_matrix ellpackMatrix, matrix vector, bool opt){

    if(ellpackMatrix.n != vector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }

    matrix result;
    result = OmpELLPACKProduct(&ellpackMatrix, &vector);

    stampaMatrice(result);
    return result;
}*/
