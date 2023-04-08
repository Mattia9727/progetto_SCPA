#include <stdio.h>
#include <stdlib.h>
#include "mmio.c"

typedef struct{
    int m;                  //Numero righe matrice
    int n;                  //Numero colonne matrice
    int nz;                 //Numero non zeri
    int*        cols;       //Vettore degli indici di colonna
    int*        rows;       //Vettore degli indici di riga
    double*     values;     //Vettore dei valori
} coo_matrix;

coo_matrix getMatrix(char* matrixFileName){
    MM_typecode matcode;
    FILE *f;
    coo_matrix mat;
    
    //Apro il file
    f = fopen(matrixFileName, "r");
    if(f == NULL) {
        printf("Could not open the file.\n");
        exit(1);
    }
    // Leggo il banner
    if (mm_read_banner(f, &matcode) != 0){
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }
    // Prendo i parametri
    if ((mm_read_mtx_crd_size(f, &(mat.m), &(mat.n), &(mat.nz))) !=0)
        exit(1);

    printf("%d %d %d\n",mat.m, mat.n, mat.nz);

    /* alloco la memoria per la matrice */

    mat.rows = (int *) malloc(mat.nz * sizeof(int));
    mat.cols = (int *) malloc(mat.nz * sizeof(int));
    mat.values = (double *) malloc(mat.nz * sizeof(double));

    for (int i=0; i<mat.nz; i++)
    {
        if(mm_is_pattern(matcode) == 1){
            fscanf(f, "%d %d\n", &(mat.rows[i]), &(mat.cols[i]));
            mat.values[i] = 1.0;
        }else{
            fscanf(f, "%d %d %lg\n", &(mat.rows[i]), &(mat.cols[i]), &(mat.values[i]));
        }
        mat.cols[i]--;  /* adjust from 1-based to 0-based */
        mat.rows[i]--;
    }
    fclose(f);

    if(mm_is_symmetric(matcode) == 1){
        printf("E' SIMMETRICA\n");
    }
 
    return mat;
}



