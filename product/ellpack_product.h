#ifndef _ELLPACKPRODUCTH_
#define _ELLPACKPRODUCTH_

#include <stdbool.h>
#include "../matrices/format/ellpack.h"

matrix calcola_prodotto_seriale_ellpack(ellpack_matrix ellpackMatrix, matrix vector, bool opt){

    if(ellpackMatrix.n != vector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }

    matrix result;
    result.m = ellpackMatrix.m;
    result.n = vector.n;
    result.coeff = (float **) malloc(sizeof(float *) * result.m);
    if (opt == true) {
        clock_t begin = clock();
        OptimizedELLPACKProduct(&ellpackMatrix, &vector, &result);
        clock_t end = clock();
        printf("\tserial_opt time: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
    }
    else {
        clock_t begin = clock();
        ELLPACKProduct(&ellpackMatrix, &vector, &result);
        clock_t end = clock();
        printf("\tserial_noopt time: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
    }

    //stampaMatrice(result);
    return result;
}

matrix calcola_prodotto_omp_ellpack(ellpack_matrix ellpackMatrix, matrix vector, bool opt){

    if(ellpackMatrix.n != vector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }

    matrix result;
    result.m = ellpackMatrix.m;
    result.n = vector.n;
    result.coeff = (float **) malloc(sizeof(float *) * result.m);

    clock_t begin = clock();
    OmpELLPACKProduct(&ellpackMatrix, &vector, &result);
    clock_t end = clock();
    printf("\tomp time: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
    //stampaMatrice(result);
    return result;
}

#endif