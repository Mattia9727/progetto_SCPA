#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include "headers/ellpack.h"

#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})

#define min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})


void print_ellpack_matrix(ellpack_matrix matrix){
    printf("M = %d, N = %d, MAXNZ = %d \n\n",matrix.m, matrix.n, matrix.maxnz);
    printf("AS (Matrice di nonzeri): \n");
    for(int i=0; i < matrix.m; i++) {
        printf("[");
        for (int j = 0; j < matrix.maxnz; j++) {
            printf("%f ",matrix.AS[i * matrix.maxnz + j]);
        }
        printf("]\n");
    }
    printf("\nJA (Indici di nonzeri): \n");
    for(int i=0; i < matrix.m; i++) {
        printf("[");
        for (int j = 0; j < matrix.maxnz; j++) {
            printf("%d ",matrix.JA[i * matrix.maxnz + j]);
        }
        printf("]\n");
    }
    printf("\n\n");
}

ellpack_matrix convert_to_ellpack(sparse_matrix* matrix){
    ellpack_matrix new_matrix;
    new_matrix.m = matrix->m;
    new_matrix.n = matrix->n;
    new_matrix.maxnz = 0;
    for(int i=0; i < matrix->m; i++){
        int maxnz=0;
        for(int j=0; j < matrix->n; j++){
            if (matrix->coeff[i][j] != 0.0) maxnz++;
        }
        if (new_matrix.maxnz < maxnz) new_matrix.maxnz = maxnz;
    }
    new_matrix.AS = (double*)calloc(new_matrix.maxnz * new_matrix.m, sizeof(double));
    new_matrix.JA = (int*)calloc(new_matrix.maxnz * new_matrix.m,sizeof(int));
    for(int i=0; i < matrix->m; i++){
    int k=0;
        for(int j=0; j < new_matrix.n; j++){
            if (matrix->coeff[i][j] != 0) {
                new_matrix.AS[i*new_matrix.maxnz+k] = matrix->coeff[i][j];
                new_matrix.JA[i*new_matrix.maxnz+k] = j;
                k++;
            }
        }
    }
    return new_matrix;
}

ellpack_matrix convert_coo_to_ellpack(coo_matrix mat){
    ellpack_matrix converted_matrix;
    converted_matrix.m = mat.m;
    converted_matrix.n = mat.n;
    converted_matrix.maxnz = 0;
    int* row_arr = (int*)calloc(mat.m,sizeof(int));
    for (int i=0; i<mat.nz; i++){
        row_arr[mat.rows[i]]++;
    }
    for (int i=0; i<mat.m; i++){
        if (converted_matrix.maxnz < row_arr[i]) converted_matrix.maxnz = row_arr[i];
    }
    free(row_arr);
    unsigned long mat_dim = (unsigned long)converted_matrix.m*(unsigned long)converted_matrix.maxnz;
    converted_matrix.AS = (double*)calloc(mat_dim,sizeof(double));
    converted_matrix.JA = (int*)calloc(mat_dim,sizeof(int));
    int row;
    int* col_arr = (int*)calloc(converted_matrix.m,sizeof(int));
    for (int i=0; i<mat.nz; i++){
        row = mat.rows[i];
        converted_matrix.AS[(unsigned long)row*converted_matrix.maxnz +col_arr[row]] = (double)mat.values[i];
        converted_matrix.JA[(unsigned long)row*converted_matrix.maxnz +col_arr[row]] = mat.cols[i];
        col_arr[row]++;
    }
    free(col_arr);

    return converted_matrix;
}

void free_ellpack_matrix(ellpack_matrix* matrix){
    free(matrix->AS);
    free(matrix->JA);
}

void print_h_ellpack_matrix(h_ellpack_matrix matrix){
    printf("M = %d, N = %d, NumMatrix = %d, HackSize = %d, Matdim = %d\n\n",matrix.m, matrix.n, matrix.numMatrix, matrix.hackSize, matrix.matDim);

    printf("\nMaxnz: \n");
    for (long i=0; i<matrix.numMatrix; i++){
        printf("%d ",matrix.maxnz[i]);
    }

    printf("\nHackOffsets: \n");
    for (long i=0; i<=matrix.numMatrix; i++){
        printf("%d ",matrix.hackOffsets[i]);
    }

    //PROVA
    //for (int i=0; i<matrix.m*matrix.maxnz[0]; i++){
    //    printf("%f \n",matrix.AS[i]);
    //}
    //è TUTTO ALLOCATO BENE, NON FUNZIONA LA PRINT
    
    for (long i=0; i<matrix.numMatrix; i++){
        printf("\nAS %d (Matrice di nonzeri): \n",i);
        for(long j=0; j < matrix.hackSize; j++) {
            printf("[");
            for (long k = 0; k < matrix.maxnz[i]; k++) {
                printf("%f ",matrix.AS[matrix.hackOffsets[i] + k + j*matrix.maxnz[i]]);
            }
            printf("]\n");
        }
        printf("\nJA %d (Indici di nonzeri): \n",i);
        for(long j=0; j < matrix.hackSize; j++) {
            printf("[");
            for (long k = 0; k < matrix.maxnz[i]; k++) {
                printf("%d ",matrix.JA[matrix.hackOffsets[i] + k + j*matrix.maxnz[i]]);
            }
            printf("]\n");
        }
        printf("\n\n");
    }
    
}

h_ellpack_matrix convert_coo_to_h_ellpack(coo_matrix mat){
    h_ellpack_matrix converted_matrix;
    converted_matrix.m = mat.m;
    converted_matrix.n = mat.n;
    converted_matrix.hackSize = 32;
    converted_matrix.numMatrix = converted_matrix.m/converted_matrix.hackSize + 1;
    converted_matrix.hackOffsets = (long*)calloc(converted_matrix.numMatrix+1,sizeof(long));
    converted_matrix.maxnz = (long*)calloc(converted_matrix.numMatrix,sizeof(long));

    long* row_arr = (long*)calloc(mat.m,sizeof(long));
    for (long i=0; i<mat.nz; i++){
        //printf("rowarr i=%ld\n",i);
        row_arr[mat.rows[i]]++;
    }
    // CALCOLO DI TUTTI I MAXNZ PER OGNI COPPIA AS E JA DI SOTTOMATRICI
    for (long i=0; i<converted_matrix.numMatrix; i++){
        //printf("maxnz i=%ld\n",i);
        long max_nz=0;
        for (long j=0; j<converted_matrix.hackSize; j++){
            
            if(i * converted_matrix.hackSize + j < converted_matrix.m) max_nz = max(max_nz,row_arr[i * converted_matrix.hackSize + j]);
            else break;
            //printf("rowarr index:%d, nz=%d, maxnz=%d\n",i * converted_matrix.hackSize + j,row_arr[mat.rows[i * converted_matrix.hackSize + j]],max_nz);
        }
        //printf("\n");
        converted_matrix.maxnz[i]= max_nz;
    }

    //for (int i=0; i<converted_matrix.numMatrix; i++){
    //    printf("maxnz: %d\n",converted_matrix.maxnz[i]);
    //}

    // CALCOLO DEGLI OFFSET TRA UNA COPPIA AS E JA DI SOTTOMATRICI E LA SUCCESSIVA
    for (long i=1; i<converted_matrix.numMatrix+1; i++){
        converted_matrix.hackOffsets[i] = converted_matrix.hackOffsets[i-1] + converted_matrix.maxnz[i-1] * converted_matrix.hackSize;
        
    }

    //for (int i=0; i<converted_matrix.numMatrix+1; i++){
    //    printf("hackOffsets: %d\n",converted_matrix.hackOffsets[i]);
    //}
    free(row_arr);


    // ALLOCAZIONE AS E JA
    long mat_dim = 0;
    for (long i=0; i<converted_matrix.numMatrix; i++){
        //printf("hackSize=%d, maxnz=%d, prod=%d\n",converted_matrix.hackSize, converted_matrix.maxnz[i],converted_matrix.hackSize * converted_matrix.maxnz[i]);
        mat_dim += converted_matrix.hackSize * converted_matrix.maxnz[i];
    }
    printf("matdim=%ld, numMatrix=%d\n",mat_dim,converted_matrix.numMatrix);
    converted_matrix.matDim = mat_dim;
    //printf("MATDIM: %d\n",mat_dim);
    converted_matrix.AS = (double*)calloc(mat_dim,sizeof(double));
    converted_matrix.JA = (long*)calloc(mat_dim,sizeof(long));

    // COSTRUZIONE AS E JA
    long row;
    long* col_arr = (long*)calloc(converted_matrix.m,sizeof(long));
    //FILE *debug_col = fopen("../measurements/results/debug_col.csv", "w+");
    //for (long i=0; i<converted_matrix.m; i++){
    //    fprintf(debug_col,"row = %ld, col= %ld\n",i,col_arr[i]);
    //}
    //FILE *debug = fopen("../measurements/results/debug.csv", "w+");
    //printf("nznznz: %d\n\n",mat.nz);
    //for (int i=0; i<mat.nz; i++){
    //    fprintf(debug,"%d: col=%d, row =%d, val=%f \n",i,mat.cols[i],mat.rows[i],mat.values[i]);
    //}        
    //FILE *debug = fopen("../measurements/results/debug_H.csv", "w+");
    long counter=0;

    for (long i=0; i<mat.nz; i++){
        //printf("i=%ld\n",i);
        long row = mat.rows[i];
        long submatrix = row/converted_matrix.hackSize;
        long subrow = row%converted_matrix.hackSize;
        converted_matrix.AS[converted_matrix.hackOffsets[row/converted_matrix.hackSize] + subrow*converted_matrix.maxnz[submatrix] + (long)col_arr[row]] = (double)mat.values[i];
        converted_matrix.JA[converted_matrix.hackOffsets[row/converted_matrix.hackSize] + subrow*converted_matrix.maxnz[submatrix] + (long)col_arr[row]] = mat.cols[i];
        col_arr[row]++;
    }

    /*for (long i=0; i<converted_matrix.numMatrix; i++){        //Per prendere hackOffsets -> focus su matrici
        printf("i=%ld\n",i);
        long hackOffsets = converted_matrix.hackOffsets[i];
        long maxnz = converted_matrix.maxnz[i];
        //printf("for (int j=%d; j<%d; j++)\n",converted_matrix.hackOffsets[i],min(converted_matrix.hackOffsets[i+1],mat.nz));
        for (long j=converted_matrix.hackOffsets[i]; j<min(converted_matrix.hackOffsets[i+1],mat.nz); j++){
            row = mat.rows[counter];
            //fprintf(debug,"row = %ld, converted_matrix.AS[%ld+%ld] = (double)mat.values[%ld]\n",row,row*maxnz,col_arr[row],counter);
            //fflush(debug);
            //printf("row = %ld, converted_matrix.AS[%ld+%ld] = (double)mat.values[%ld]\n",row,row*maxnz,col_arr[row],counter);
            converted_matrix.AS[row*maxnz + (long)col_arr[row]] = (double)mat.values[counter];
            converted_matrix.JA[row*maxnz + (long)col_arr[row]] = mat.cols[counter];
            col_arr[row]++;
            counter++;
        }
    }*/
    free(col_arr);
    //for (int i=0; i<mat_dim; i++){
    //    printf("%d ",converted_matrix.AS[i]);
    //}

    return converted_matrix;
}

h_ellpack_matrix_bis convert_coo_to_h_ellpack_bis(coo_matrix mat){
    h_ellpack_matrix_bis converted_matrix;
    converted_matrix.m = mat.m;
    converted_matrix.n = mat.n;
    converted_matrix.hackSize = 32;
    converted_matrix.numMatrix = converted_matrix.m/converted_matrix.hackSize + 1;
    if(converted_matrix.m%converted_matrix.hackSize == 0) converted_matrix.numMatrix--;
    converted_matrix.hackOffsets = (int*)calloc(converted_matrix.numMatrix+1,sizeof(int));
    converted_matrix.maxnz = (int*)calloc(converted_matrix.numMatrix,sizeof(int));

    int* row_arr = (int*)calloc(mat.m,sizeof(int));
    for (int i=0; i<mat.nz; i++){
        //printf("rowarr i=%ld\n",i);
        row_arr[mat.rows[i]]++;
    }
    // CALCOLO DI TUTTI I MAXNZ PER OGNI COPPIA AS E JA DI SOTTOMATRICI
    for (int i=0; i<converted_matrix.numMatrix; i++){
        //printf("maxnz i=%ld\n",i);
        int max_nz=0;
        for (int j=0; j<converted_matrix.hackSize; j++){
            
            if(i * converted_matrix.hackSize + j < converted_matrix.m) max_nz = max(max_nz,row_arr[i * converted_matrix.hackSize + j]);
            else break;
            //printf("rowarr index:%d, nz=%d, maxnz=%d\n",i * converted_matrix.hackSize + j,row_arr[mat.rows[i * converted_matrix.hackSize + j]],max_nz);
        }
        //printf("\n");
        converted_matrix.maxnz[i]= max_nz;
    }

    //for (int i=0; i<converted_matrix.numMatrix; i++){
    //    printf("maxnz: %d\n",converted_matrix.maxnz[i]);
    //}

    // CALCOLO DEGLI OFFSET TRA UNA COPPIA AS E JA DI SOTTOMATRICI E LA SUCCESSIVA
    for (int i=1; i<converted_matrix.numMatrix+1; i++){
        converted_matrix.hackOffsets[i] = converted_matrix.hackOffsets[i-1] + converted_matrix.maxnz[i-1] * converted_matrix.hackSize;
        
    }

    //for (int i=0; i<converted_matrix.numMatrix+1; i++){
    //    printf("hackOffsets: %d\n",converted_matrix.hackOffsets[i]);
    //}
    free(row_arr);


    // ALLOCAZIONE AS E JA
    int mat_dim = 0;
    for (int i=0; i<converted_matrix.numMatrix; i++){
        //printf("hackSize=%d, maxnz=%d, prod=%d\n",converted_matrix.hackSize, converted_matrix.maxnz[i],converted_matrix.hackSize * converted_matrix.maxnz[i]);
        mat_dim += converted_matrix.hackSize * converted_matrix.maxnz[i];
    }
    printf("matdim=%ld, numMatrix=%d\n",mat_dim,converted_matrix.numMatrix);
    converted_matrix.matDim = mat_dim;
    //printf("MATDIM: %d\n",mat_dim);
    converted_matrix.AS = (double*)calloc(mat_dim,sizeof(double));
    converted_matrix.JA = (int*)calloc(mat_dim,sizeof(int));

    // COSTRUZIONE AS E JA
    int row;
    int* col_arr = (int*)calloc(converted_matrix.m,sizeof(int));
    //FILE *debug_col = fopen("../measurements/results/debug_col.csv", "w+");
    //for (long i=0; i<converted_matrix.m; i++){
    //    fprintf(debug_col,"row = %ld, col= %ld\n",i,col_arr[i]);
    //}
    //FILE *debug = fopen("../measurements/results/debug.csv", "w+");
    //printf("nznznz: %d\n\n",mat.nz);
    //for (int i=0; i<mat.nz; i++){
    //    fprintf(debug,"%d: col=%d, row =%d, val=%f \n",i,mat.cols[i],mat.rows[i],mat.values[i]);
    //}        
    //FILE *debug = fopen("../measurements/results/debug_H.csv", "w+");
    int counter=0;

    for (int i=0; i<mat.nz; i++){
        //printf("i=%ld\n",i);
        int row = mat.rows[i];
        int submatrix = row/converted_matrix.hackSize;
        int subrow = row%converted_matrix.hackSize;
        converted_matrix.AS[converted_matrix.hackOffsets[row/converted_matrix.hackSize] + subrow*converted_matrix.maxnz[submatrix] + col_arr[row]] = (double)mat.values[i];
        converted_matrix.JA[converted_matrix.hackOffsets[row/converted_matrix.hackSize] + subrow*converted_matrix.maxnz[submatrix] + col_arr[row]] = mat.cols[i];
        col_arr[row]++;
    }

    
    free(col_arr);


    return converted_matrix;
}

void free_h_ellpack_matrix(h_ellpack_matrix* matrix){
    free(matrix->AS);
    free(matrix->JA);
    free(matrix->maxnz);
    free(matrix->hackOffsets);
}

void free_h_ellpack_matrix(h_ellpack_matrix_bis* matrix){
    free(matrix->AS);
    free(matrix->JA);
    free(matrix->maxnz);
    free(matrix->hackOffsets);
}


void fprint_ellpack_matrix(ellpack_matrix matrix){

    FILE *ellp = fopen("../measurements/results/ellp.csv", "w+");

    fprintf(ellp,"M = %d, N = %d, MAXNZ = %d \n\n",matrix.m, matrix.n, matrix.maxnz);
    fprintf(ellp,"AS (Matrice di nonzeri): \n");
    for(int i=0; i < matrix.m; i++) {
        fprintf(ellp,"[");
        for (int j = 0; j < matrix.maxnz; j++) {
            fprintf(ellp,"%f ",matrix.AS[i * matrix.maxnz + j]);
        }
        fprintf(ellp,"]\n");
    }
    fprintf(ellp,"\nJA (Indici di nonzeri): \n");
    for(int i=0; i < matrix.m; i++) {
        fprintf(ellp,"[");
        for (int j = 0; j < matrix.maxnz; j++) {
            fprintf(ellp,"%d ",matrix.JA[i * matrix.maxnz + j]);
        }
        fprintf(ellp,"]\n");
    }
    fprintf(ellp,"\n\n");
    fclose(ellp);
}

void fprint_h_ellpack_matrix(h_ellpack_matrix matrix){

    FILE *hellp = fopen("../measurements/results/hellp.csv", "w+");

    fprintf(hellp,"M = %ld, N = %ld, NumMatrix = %ld, HackSize = %ld\n\n",matrix.m, matrix.n, matrix.numMatrix, matrix.hackSize);

    fprintf(hellp,"\nMaxnz: \n");
    for (int i=0; i<matrix.numMatrix; i++){
        fprintf(hellp,"%ld ",matrix.maxnz[i]);
    }

    fprintf(hellp,"\nHackOffsets: \n");
    for (int i=0; i<matrix.numMatrix; i++){
        fprintf(hellp,"%ld ",matrix.hackOffsets[i]);
    }

    //PROVA
    //for (int i=0; i<matrix.m*matrix.maxnz[0]; i++){
    //    printf("%f \n",matrix.AS[i]);
    //}
    //è TUTTO ALLOCATO BENE, NON FUNZIONA LA PRINT
    for (int i=0; i<matrix.numMatrix; i++){
        fprintf(hellp,"\nAS %ld (Matrice di nonzeri): \n",i);
        for(int j=0; j < matrix.hackSize; j++) {
            fprintf(hellp,"[");
            for (int k = 0; k < matrix.maxnz[i]; k++) {
                fprintf(hellp,"%f ",matrix.AS[matrix.hackOffsets[i] + k + j*matrix.maxnz[i]]);
            }
            fprintf(hellp,"]\n");
        }
        fprintf(hellp,"\nJA %ld (Indici di nonzeri): \n",i);
        for(int j=0; j < matrix.hackSize; j++) {
            fprintf(hellp,"[");
            for (int k = 0; k < matrix.maxnz[i]; k++) {
                fprintf(hellp,"%ld ",matrix.JA[matrix.hackOffsets[i] + k + j*matrix.maxnz[i]]);
            }
            fprintf(hellp,"]\n");
        }
        fprintf(hellp,"\n\n");
    }
    fclose(hellp);
}

void fprint_h_ellpack_matrix_bis(h_ellpack_matrix_bis matrix){

    FILE *hellp = fopen("../measurements/results/hellp.csv", "a");

    fprintf(hellp,"M = %d, N = %d, NumMatrix = %d, HackSize = %d\n\n",matrix.m, matrix.n, matrix.numMatrix, matrix.hackSize);
    int count = 0;
    for (int i=0; i<matrix.numMatrix; i++){
        if(matrix.maxnz[i] > 32 && matrix.maxnz[i] < 96) count++;
    }
    fprintf(hellp,"%d ",count);

    /*
    fprintf(hellp,"\nMaxnz: \n");
    for (int i=0; i<matrix.numMatrix; i++){
        fprintf(hellp,"%d ",matrix.maxnz[i]);
    }
    
    fprintf(hellp,"\nHackOffsets: \n");
    for (int i=0; i<matrix.numMatrix; i++){
        fprintf(hellp,"%d ",matrix.hackOffsets[i]);
    }
    */
    //PROVA
    //for (int i=0; i<matrix.m*matrix.maxnz[0]; i++){
    //    printf("%f \n",matrix.AS[i]);
    //}
    //è TUTTO ALLOCATO BENE, NON FUNZIONA LA PRINT
    /*
    for (int i=0; i<matrix.numMatrix; i++){
        fprintf(hellp,"\nAS %d (Matrice di nonzeri): \n",i);
        for(int j=0; j < matrix.hackSize; j++) {
            fprintf(hellp,"[");
            for (int k = 0; k < matrix.maxnz[i]; k++) {
                fprintf(hellp,"%f ",matrix.AS[matrix.hackOffsets[i] + k + j*matrix.maxnz[i]]);
            }
            fprintf(hellp,"]\n");
        }
        fprintf(hellp,"\nJA %d (Indici di nonzeri): \n",i);
        for(int j=0; j < matrix.hackSize; j++) {
            fprintf(hellp,"[");
            for (int k = 0; k < matrix.maxnz[i]; k++) {
                fprintf(hellp,"%d ",matrix.JA[matrix.hackOffsets[i] + k + j*matrix.maxnz[i]]);
            }
            fprintf(hellp,"]\n");
        }
        fprintf(hellp,"\n\n");
    }*/
    fclose(hellp);
}