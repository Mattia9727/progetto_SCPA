#ifndef _CRSH_
#define _CSRH_ 

#include "coo.h"
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


typedef struct{
    int m;              //Numero righe matrice
    int n;              //Numero colonne matrice
    int nz;             //Numero non zeri
    int*        irp;    //Vettore dei puntatori all'inizio di ciascuna riga
    double*      as;    //Vettore dei coefficienti
    int*        ja;     //Vettore degli indici di colonna
} csr_matrix;

void stampa_matrice_csr(csr_matrix mat);
csr_matrix convert_to_csr(sparse_matrix mat);
csr_matrix convert_to_csr_from_coo(coo_matrix mat);
void free_csr_matrix(csr_matrix* matrix);

#endif