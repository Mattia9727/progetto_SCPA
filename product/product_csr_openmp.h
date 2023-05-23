
#ifndef _PRODUCTCSROPENMPH_
#define _PRODUCTCSROPENMPH_ 

#include <stdio.h>
#include <omp.h>
#include "product.h"
#include <time.h>

#define NUM_THREADS 40


double calcola_prodotto_per_righe_csr_openmp(csr_matrix csrMatrix, matrix multivector,matrix* result){
    /*
    if(csrMatrix.n != multivector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }
    */
    double partialSum;

    int i,j,k, irp_1,irp_2;
    int m = csrMatrix.m,n=multivector.n;
    struct timespec start, end;
    int chunkSize = ((int)((csrMatrix.m/(float)NUM_THREADS)/16))*16;
    if(chunkSize < 16) chunkSize = 16;
    clock_gettime(CLOCK_MONOTONIC, &start);
    #pragma omp parallel for schedule(static,chunkSize) firstprivate(n,m) private(partialSum,i,j,k, irp_1, irp_2)
    for(i = 0; i < m; i++){
        irp_1 = csrMatrix.irp[i];
        irp_2 = csrMatrix.irp[i+1];
        
        for(k = 0; k < n; k++){
            
            partialSum = 0;
            #pragma omp reduction(+ : partialSum)
            for(j = irp_1; j < irp_2; j++){
                partialSum += csrMatrix.as[j]*multivector.coeff[csrMatrix.ja[j]* n + k];
            }
            result->coeff[i * n + k] = partialSum;
            
        }
    }    
    clock_gettime(CLOCK_MONOTONIC, &end);
    return (double)( end.tv_sec - start.tv_sec )+ ( end.tv_nsec - start.tv_nsec )/ (double)1000000000L;
} 

double calcola_prodotto_per_righe_csr_openmp_bis(csr_matrix csrMatrix, matrix multivector,matrix* result){
    /*
    if(csrMatrix.n != multivector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }
    */
    int i,j,k, irp_1,irp_2,col;
    double elem;
    int m = csrMatrix.m,n=multivector.n;
    struct timespec start, end;
    int chunkSize = ((int)((csrMatrix.m/(float)NUM_THREADS)/16))*16;
    if(chunkSize < 16) chunkSize = 16;
    clock_gettime(CLOCK_MONOTONIC, &start);
    #pragma omp parallel for schedule(static,chunkSize) firstprivate(n,m) private(elem, col,i,j,k, irp_1, irp_2)
    for(i = 0; i < m; i++){
        irp_1 = csrMatrix.irp[i];
        irp_2 = csrMatrix.irp[i+1];
        for(j = irp_1; j < irp_2; j++){
            elem = csrMatrix.as[j];
            col = csrMatrix.ja[j];
            for(int k = 0; k < n; k++){                
                result->coeff[i * n + k] += elem*multivector.coeff[col * n + k];
            }
        }
    }    
    clock_gettime(CLOCK_MONOTONIC, &end);
    return (double)( end.tv_sec - start.tv_sec )+ ( end.tv_nsec - start.tv_nsec )/ (double)1000000000L;
} 
/*
double calcola_prodotto_per_righe_csr_openmp_bis(csr_matrix csrMatrix, matrix multivector,matrix* result){
    
    if(csrMatrix.n != multivector.m){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }
    
    double elem;
    int col;
    int i;

    int chunckSize = ((int)((csrMatrix.m/(float)NUM_THREADS)/64))*64;
    if(chunckSize < 64) chunckSize = 64;
    clock_t begin, end;
    begin = clock();
    #pragma omp parallel for schedule(static,chunckSize) shared(result, multivector, csrMatrix) private(i,col,elem)
    for(i = 0; i < csrMatrix.m; i++){
        for(int j = csrMatrix.irp[i]; j < csrMatrix.irp[i+1]; j++){
            elem = csrMatrix.as[j];
            col = csrMatrix.ja[j];
            for(int k = 0; k < multivector.n; k++){                
                result->coeff[i][k] += elem*multivector.coeff[col][k];
            }
        }
    }    
    end = clock();
    return (double)(end - begin) / CLOCKS_PER_SEC;
}  
*/
double calcola_prodotto_per_righe_csr_openmp_trasposto(csr_matrix csrMatrix, matrix multivector_trasposto,matrix* result){
    if(csrMatrix.n != multivector_trasposto.n){
        printf("Prodotto non calcolabile tra la matrice e il multivettore inserito\n");
        exit(1);
    }

    double partialSum;

    int i,j,k, irp_1,irp_2;
    int mat_m = csrMatrix.m,vec_n=multivector_trasposto.n, vec_m = multivector_trasposto.m;
    struct timespec start, end;
    
    int chunkSize = ((int)((csrMatrix.m/(float)NUM_THREADS)/16))*16;
    if(chunkSize < 16) chunkSize = 16;
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    #pragma omp parallel for schedule(static,chunkSize) firstprivate(vec_n,vec_m, mat_m) private(partialSum,i,j,k, irp_1, irp_2)
    for(i = 0; i < mat_m; i++){
        irp_1 = csrMatrix.irp[i];
        irp_2 = csrMatrix.irp[i+1];
        for(k = 0; k < vec_m; k++){
            partialSum = 0;
            #pragma omp reduction(+ : partialSum)
            for(j = irp_1; j < irp_2; j++){
                partialSum += csrMatrix.as[j]*multivector_trasposto.coeff[k*vec_n + csrMatrix.ja[j]];
            }
            result->coeff[i * vec_m + k] = partialSum;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    return (double)( end.tv_sec - start.tv_sec )+ ( end.tv_nsec - start.tv_nsec )/ (double)1000000000L;
}
#endif