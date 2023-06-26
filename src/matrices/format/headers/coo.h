#ifndef _COOH_
#define _COOH_ 

#include "../../headers/matrix_generator.h"

typedef struct{
    int m;                  //Numero righe matrice
    int n;                  //Numero colonne matrice
    int nz;                 //Numero non zeri
    int*        cols;       //Vettore degli indici di colonna
    int*        rows;       //Vettore degli indici di riga
    double*     values;     //Vettore dei valori
} coo_matrix;

void stampa_matrice_coo(coo_matrix mat);

matrix convert_to_simple_matrix(coo_matrix mat);

coo_matrix get_matrix(char* matrixFileName);


void free_coo_matrix(coo_matrix* matrix);

#endif