#ifndef _ELLPACKPRODUCTH_
#define _ELLPACKPRODUCTH_

#include <stdbool.h>
#include "../matrices/format/ellpack.h"


void ellpack_product(ellpack_matrix* mat, matrix* vector, matrix* result){

    //int op = 0;
    for (int i = 0; i < result->m; i++) {
        result->coeff[i] = (double *) calloc(result->n, sizeof(double));
        for (int k = 0; k < result->n; k++) {
            double t = 0;
            for (int j = 0; j < mat->maxnz; j++) {
                t = t + mat->AS[i][j]*vector->coeff[mat->JA[i][j]][k];
                //op++;
            }
            result->coeff[i][k] = t;
        }
    }
    //printf("\n\nOperazioni effettuate: %d\n\n",op); //Esegue matrix->M * vector->n * MAXNZ operazioni :)
}

void omp_ellpack_product(ellpack_matrix mat, matrix vector, matrix* result){

    //int op = 0;

    double t;

    int i;

    #pragma omp parallel for schedule(static, 1) shared(result, mat, vector) private(t,i)
    for (i = 0; i < result->m; i++) {
        for (int k = 0; k < result->n; k++) {
            t = 0;
            for (int j = 0; j < mat.maxnz; j++) {
                t = t + mat.AS[i][j]*vector.coeff[mat.JA[i][j]][k];
                //op++;
            }
            result->coeff[i][k] = t;
        }
    }
    //printf("\n\nOperazioni effettuate: %d\n\n",op); //Esegue matrix->M * vector->n * MAXNZ operazioni :)
}

void optimized_ellpack_product(ellpack_matrix* mat, matrix* vector, matrix* result){

    int op = 0;
    for (int i = 0; i < result->m; i++) {
        result->coeff[i] = (double *) calloc(result->n, sizeof(double));
        for (int k = 0; k < result->n; k++) {
            double t = 0;
            int prev_JA = -1;
            for (int j = 0; j < mat->maxnz; j++) {
                if (prev_JA>=mat->JA[i][j]) break;
                prev_JA = mat->JA[i][j];
                t = t + mat->AS[i][j]*vector->coeff[prev_JA][k];
                op++;
            }
            result->coeff[i][k] = t;
        }
    }
    printf("\n\nOperazioni effettuate: %d\n\n",op); //Esegue matrix->M * vector->n * MAXNZ operazioni :)
}

matrix calcola_prodotto_seriale_ellpack(ellpack_matrix ellpackMatrix, matrix vector, bool opt){

    if(ellpackMatrix.n != vector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }

    matrix result;
    result.m = ellpackMatrix.m;
    result.n = vector.n;
    result.coeff = (double **) malloc(sizeof(double *) * result.m);
    if (opt == true) {
        clock_t begin = clock();
        optimized_ellpack_product(&ellpackMatrix, &vector, &result);
        clock_t end = clock();
        printf("\tserial_opt time: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
    }
    else {
        clock_t begin = clock();
        ellpack_product(&ellpackMatrix, &vector, &result);
        clock_t end = clock();
        printf("\tserial_noopt time: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
    }

    //stampa_matrice(result);
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
    result.coeff = (double **) malloc(sizeof(double *) * result.m);

    clock_t begin = clock();
    omp_ellpack_product(ellpackMatrix, vector, &result);
    clock_t end = clock();
    printf("\tomp time: %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
    //stampa_matrice(result);
    return result;
}

#endif