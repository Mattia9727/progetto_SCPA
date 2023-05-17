#ifndef _COOH_
#define _COOH_ 

#include <stdlib.h>
#include "../matrix_generator.h"

typedef struct{
    int m;                  //Numero righe matrice
    int n;                  //Numero colonne matrice
    int nz;                 //Numero non zeri
    int*        cols;       //Vettore degli indici di colonna
    int*        rows;       //Vettore degli indici di riga
    double*     values;     //Vettore dei valori
} coo_matrix;

void stampa_matrice_coo(coo_matrix mat){
    for(int i = 0; i < mat.nz; i++){
            printf("%d %d %lf\n",mat.rows[i], mat.cols[i],mat.values[i]);
    }
}

matrix convert_to_simple_matrix(coo_matrix mat){
    matrix new_matrix;
    new_matrix.m = mat.m;
    new_matrix.n = mat.n;
    new_matrix.coeff = (double **)malloc(sizeof(float*)*mat.m);
    for(int i = 0; i < mat.m; i++){
        new_matrix.coeff[i] = (double *) calloc(mat.n, sizeof(double));     // allocazione della memoria per la riga i-esima
        if (new_matrix.coeff[i] == NULL) {
            printf("Errore di allocazione della memoria.\n");
            exit(0);
        }
    }
    if (new_matrix.coeff == NULL) {
        printf("Errore di allocazione della memoria.\n");
        exit(0);
    }
    for(int i = 0; i < mat.nz; i++){
        new_matrix.coeff[mat.rows[i]][mat.cols[i]] = (double)(mat.values[i]);
    }
    return new_matrix;
}

#endif