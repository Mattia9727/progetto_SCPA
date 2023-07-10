#include <stdbool.h>
#include "headers/product_ellpack.h"
#include <time.h>

void ellpack_product(ellpack_matrix* mat, matrix* vector, matrix* result){
    for (int i = 0; i < result->m; i++) {
        for (int k = 0; k < result->n; k++) {
            double t = 0;
            for (int j = 0; j < mat->maxnz; j++) {
                t = t + mat->AS[i*mat->maxnz+j]*vector->coeff[(mat->JA[i*mat->maxnz+j])*vector->n + k];
            }
            result->coeff[i * result->n + k] = t;
        }
    }
}

double omp_ellpack_product(ellpack_matrix mat, matrix vector, matrix* result, int nThreads){
    double t;

    int i,j,k;
    unsigned long maxnz= mat.maxnz, m = result->m, n= result->n;
    struct timespec start, end;
    int chunkSize = ((int)((mat.m/(float)nThreads)/16))*16;
    if(chunkSize < 16) chunkSize = 16;
    clock_gettime(CLOCK_MONOTONIC, &start);
    #pragma omp parallel for schedule(static,chunkSize) shared(result, mat, vector) firstprivate(m,n,maxnz) private(t,i,j,k)
    for (i = 0; i < m; i++) {
        for (k = 0; k < n; k++) {
            t = 0;
            #pragma omp reduction(+ : t)
            for (j = 0; j < maxnz; j++) {
                t = t + mat.AS[i*maxnz+j]*vector.coeff[(mat.JA[i*maxnz+j])*n + k];
            }
            result->coeff[i * n + k] = t;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    return (double)( end.tv_sec - start.tv_sec )+ ( end.tv_nsec - start.tv_nsec )/ (double)1000000000L;
}

double optimized_omp_ellpack_product(ellpack_matrix mat, matrix vector, matrix* result, int nThreads){
    double t = 0;
    int i,j,k;
    unsigned long maxnz= mat.maxnz, m = result->m, n= result->n;
    int prev_JA = -1;
    int curr_JA;
    struct timespec start, end;
    int chunkSize = ((int)((mat.m/(float)nThreads)/16))*16;
    if(chunkSize < 16) chunkSize = 16;
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    #pragma omp parallel for schedule(static,chunkSize) num_threads(nThreads) shared(result, mat, vector) firstprivate(m,n,maxnz) private(t,i,j,k,prev_JA, curr_JA)
    for (i = 0; i < m; i++) {
        for (int k = 0; k < n; k++) {
            prev_JA = -1;
            t = 0;
            for (int j = 0; j < maxnz; j++) {
                curr_JA = mat.JA[i*maxnz+j];
                if (prev_JA<curr_JA){
                    prev_JA = curr_JA;
                    t = t + mat.AS[i*maxnz+j]*vector.coeff[curr_JA*n+k];
                }else break;
            
            }
            result->coeff[i*n+k] = t;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    return (double)( end.tv_sec - start.tv_sec )+ ( end.tv_nsec - start.tv_nsec )/ (double)1000000000L;
}
