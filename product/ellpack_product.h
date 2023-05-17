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
    double t;

    int i;

    #pragma omp parallel for schedule(static) shared(result, mat, vector) private(t,i)
    for (i = 0; i < result->m; i++) {
        for (int k = 0; k < result->n; k++) {
            t = 0;
            for (int j = 0; j < mat.maxnz; j++) {
                t = t + mat.AS[i][j]*vector.coeff[mat.JA[i][j]][k];
            }
            result->coeff[i][k] = t;
        }
    }
}

void optimized_ellpack_product(ellpack_matrix mat, matrix vector, matrix* result){
    double t = 0;
    int i;
    int prev_JA = -1;
    int curr_JA;
    #pragma omp parallel for schedule(static) shared(result, mat, vector) private(t,i,prev_JA, curr_JA)
    for (i = 0; i < result->m; i++) {
        for (int k = 0; k < result->n; k++) {
            prev_JA = -1;
            t = 0;
            for (int j = 0; j < mat.maxnz; j++) {
                curr_JA = mat.JA[i][j];
                if (prev_JA<curr_JA){
                    prev_JA = curr_JA;
                    t = t + mat.AS[i][j]*vector.coeff[curr_JA][k];
                }else j = mat.maxnz;

            }
            result->coeff[i][k] = t;
        }
    }
}


#endif