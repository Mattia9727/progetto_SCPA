#ifndef _MATRIXRETRIEVERH_
#define _MATRIXRETRIEVERH_

#include <stdio.h>
#include <stdlib.h>
#include "mmio.c"
#include "format/coo.h"


coo_matrix getMatrix(char* matrixFileName){
    MM_typecode matcode;
    FILE *f;
    coo_matrix mat, total_mat;
    int is_pattern, diag_nz = 0;
    
    
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

    /* alloco la memoria per la matrice */

    mat.rows = (int *) malloc(mat.nz * sizeof(int));
    mat.cols = (int *) malloc(mat.nz * sizeof(int));
    mat.values = (double *) malloc(mat.nz * sizeof(double));

    is_pattern = mm_is_pattern(matcode);


    for (int i=0; i<mat.nz; i++)
    {
        if(is_pattern){
            fscanf(f, "%d %d\n", &(mat.rows[i]), &(mat.cols[i]));
            mat.values[i] = 1.0;
        }else{
            fscanf(f, "%d %d %lg\n", &(mat.rows[i]), &(mat.cols[i]), &(mat.values[i]));
        }
        mat.cols[i]--;  /* adjust from 1-based to 0-based */
        mat.rows[i]--;
        if(mat.cols[i] == mat.rows[i])
            diag_nz++;
    }
    fclose(f);

    if(mm_is_symmetric(matcode) == 1){
        total_mat.m = mat.m;
        total_mat.n = mat.n;
        total_mat.nz = (mat.nz - diag_nz) * 2 + diag_nz;
        total_mat.cols = (int*)malloc(total_mat.nz * sizeof(int));
        total_mat.rows = (int*)malloc(total_mat.nz * sizeof(int));
        total_mat.values = (double*)malloc(total_mat.nz * sizeof(double));
        int j = 0;
        for(int i = 0; i < mat.nz; i++){
            total_mat.cols[j] = mat.cols[i];
            total_mat.rows[j] = mat.rows[i];
            total_mat.values[j] = mat.values[i];
            if(mat.cols[i] != mat.rows[i]){
                j++;
                total_mat.cols[j] = mat.rows[i];
                total_mat.rows[j] = mat.cols[i];
                total_mat.values[j] = mat.values[i];
            }
            j++;
        }
        return total_mat;

    }

    return mat;
}

#endif



