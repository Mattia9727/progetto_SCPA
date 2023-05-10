#ifndef _COOH_
#define _COOH_ 

#include <stdlib.h>

typedef struct{
    int m;                  //Numero righe matrice
    int n;                  //Numero colonne matrice
    int nz;                 //Numero non zeri
    int*        cols;       //Vettore degli indici di colonna
    int*        rows;       //Vettore degli indici di riga
    double*     values;     //Vettore dei valori
} coo_matrix;

void stampaMatriceCoo(coo_matrix mat){
    for(int i = 0; i < mat.nz; i++){
            printf("%d %d %lf\n",mat.rows[i], mat.cols[i],mat.values[i]);
    }
}

#endif